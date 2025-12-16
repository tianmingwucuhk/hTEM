*********************************************************************
**********************Day6_project_QC&integration********************
*********************************************************************

library(Seurat)
library(tidyverse)
library(ggplot2)
library(scutilsR)
library(gridExtra)
library(harmony)
library(ggplot2)

samples <- list.dirs("data/10x_matrix/", full.names = F, recursive = F)
seu.list <- pbapply::pblapply(samples, function(sn) {
  counts <- Read10X(file.path("data/10x_matrix/", sn))
  sn <- gsub("_", "-", sn) # 注意"_"在`CreateSeuratObject()`里有特殊的意义
  colnames(counts) <- paste(sn, colnames(counts), sep = "_")
  seu <- CreateSeuratObject(counts = counts,min.cells = 3, min.features = 200)
  return(seu)
})

D6 <- base::Reduce(f = merge, x = seu.list)

rm(seu.list)
gc()

table(D6$orig.ident)
#####qc###########
D6 <- PercentageFeatureSet(D6, pattern = "^MT-", col.name = "percent.mt")
D6 <- PercentageFeatureSet(D6, pattern = "^RPS", col.name = "percent.rps")
D6 <- PercentageFeatureSet(D6, pattern = "^RPL", col.name = "percent.rpl")

VlnPlot(D6, features = c("nFeature_RNA"), group.by = "orig.ident", pt.size = 0) + NoLegend() +
  geom_hline(yintercept = c(1000, 5000), linetype = "dashed", color = "blue")

VlnPlot(D6, features = c("nCount_RNA"), group.by = "orig.ident", pt.size = 0, log = T) + NoLegend() + 
  geom_hline(yintercept = c(3000, 50000), linetype = "dashed", color = "blue") + NoLegend()

VlnPlot(D6, features = c("percent.mt"), group.by = "orig.ident", pt.size = 0, log = T) +
  geom_hline(yintercept = c(10), linetype = "dashed", color = "blue")


D6 <- subset(D6, subset = nFeature_RNA > 200 & percent.mt < 5)

D6 <- NormalizeData(D6)
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
D6 <- CellCycleScoring(D6, s.features = s.genes, g2m.features = g2m.genes)


D6 <- scutilsR::MarkDoublets(D6, split.by = "orig.ident")
D6 <- subset(D6,subset = DF.classifications=='Singlet')

QuickCluster <- function(object) {
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object, nfeatures = 2000)
  object <- ScaleData(object)
  object <- RunPCA(object)
  object <- FindNeighbors(object, reduction = "pca", dims = 1:30)
  object <- FindClusters(object)
  return(object)
}
seu.list <- SplitObject(D6, split.by = "orig.ident")
seu.list <- lapply(seu.list, QuickCluster)
clusters <- lapply(seu.list, function(xx) xx$seurat_clusters) %>% base::Reduce(c, .)
D6$quick_clusters <- clusters[rownames(D6@meta.data)]
rm(seu.list)
gc()

D6 <- scutilsR::RemoveAmbientRNAs(D6, split.by = "orig.ident", cluster.name = "quick_clusters")
D6
head(D6@meta.data)
DefaultAssay(D6) <- 'decontX'

D6 <- PercentageFeatureSet(D6, pattern = "^MT-", col.name = "percent.mt")
D6 <- PercentageFeatureSet(D6, pattern = "^RPS", col.name = "percent.rps")
D6 <- PercentageFeatureSet(D6, pattern = "^RPL", col.name = "percent.rpl")
D6 <- subset(D6, subset = nFeature_decontX > 200 & percent.mt < 5)
qs::qsave(D6, file = 'output/D6_downstream.qs')


######integration##########
##
D6_perturbation <- qs::qread(file = 'output/D6_downstream.qs')
v3 <- qs::qread(file = 'data/01_NTSM_integrated_decontX.qs')
v3_d6<-subset(v3, subset = orig.ident == 'NTSM-V4-D6')

#fix cellid
Newcellnames <- gsub('NTSM-V4-D6','v3-d6',colnames(v3_d6))
v3_d6 <- RenameCells(v3_d6, new.names = Newcellnames)
v3_d6$orig.ident <- as.factor(v3_d6$orig.ident)
v3_d6$orig.ident <- factor(v3_d6$orig.ident,labels = 'v3_d6')

D6_perturbation[['RNA']] <- NULL
data <- list(v3_d6, D6_perturbation)
D6_integration <- merge(x=data[[1]],y=data[-1])


##
D6_integration <- NormalizeData(D6_integration)
D6_integration <- CellCycleScoring(D6_integration, s.features = s.genes, g2m.features = g2m.genes)
D6_integration <- FindVariableFeatures(D6_integration, selection.method = "vst", nfeatures = 2000)
D6_integration <- ScaleData(D6_integration, vars.to.regress = c('S.Score','G2M.Score'))
D6_integration <- RunPCA(D6_integration, features = VariableFeatures(object = D6_integration),reduction.name = "pca")


#harmony
#integrate with harmony
library(harmony)
D6_integration <- RunHarmony(D6_integration, group.by.vars=c('orig.ident'),reduction.use='pca',dim.use=1:50)
D6_integration <- RunUMAP(D6_integration,reduction = 'harmony',dims = 1:30,n.neighbors = 200,reduction.name = 'umap.harmony')

#bbknn
library(bbknnR)
D6_integration=RunBBKNN(D6_integration,reduction = 'pca',run_UMAP = TRUE, UMAP_name="umap.bbknn",run_TSNE= FALSE,
                        batch_key = 'orig.ident')

p1 =DimPlot(D6_integration,reduction = 'umap.harmony',group.by = 'orig.ident')
p2 =DimPlot(D6_integration,reduction = 'umap.bbknn', group.by = 'orig.ident')
p1+p2

FeaturePlot(D6_integration, features = c('PAX3'), pt.size = 0.2,reduction = "umap.harmony",split.by = 'orig.ident',
            ncol = 9, order = TRUE)

DimPlot(D6_integration,group.by = 'detailed_celltype',split.by = 'orig.ident',reduction = 'umap.bbknn')
ggsave(filename = 'D6_integration_detailed_celltype.tiff',width = 13,height = 5,dpi = 300,path = 'output/Figures/')
qs::qsave(D6_integration, file = 'output/D6_integration.qs')
D6_integration <- qs::qread(file = 'output/D6_integration.qs')

resolutions <- c(1.2,1.4,1.6,1.8,2.0,2.4,2.8,3.0)

D6_integration <- FindClusters(D6_integration, resolution = resolutions,graph.name = 'bbknn')
DimPlot(D6_integration,reduction = 'umap.bbknn',group.by =  paste0("bbknn_res.", resolutions), ncol = 3, label = T) & NoLegend()
ggsave(filename = 'D6_integration_bbknn_clusters.tiff',width = 14,height =11,dpi = 300,path = 'output/Figures/')


*********************************************************************
**********************Day7_project_QC&integration********************
*********************************************************************

library(Seurat)
library(tidyverse)
library(scutilsR)

samples <- list.dirs("data/10x_matrix/", full.names = F, recursive = F)
seu.list <- pbapply::pblapply(samples, function(sn) {
  counts <- Read10X(file.path("data/10x_matrix/", sn))
  sn <- gsub("_", "-", sn) 
  colnames(counts) <- paste(sn, colnames(counts), sep = "_")
  seu <- CreateSeuratObject(counts = counts,min.cells = 3, min.features = 200)
  return(seu)
})

D7 <- base::Reduce(f = merge, x = seu.list)

rm(seu.list)
gc()

#### QC ####
table(D7$orig.ident)

D7 <- PercentageFeatureSet(D7, pattern = "^MT-", col.name = "percent.mt")
D7 <- PercentageFeatureSet(D7, pattern = "^RPS", col.name = "percent.rps")
D7 <- PercentageFeatureSet(D7, pattern = "^RPL", col.name = "percent.rpl")

VlnPlot(D7, features = c("nFeature_RNA"), group.by = "orig.ident", pt.size = 0) + NoLegend() +
  geom_hline(yintercept = c(1000, 5000), linetype = "dashed", color = "blue")

VlnPlot(D7, features = c("nCount_RNA"), group.by = "orig.ident", pt.size = 0, log = T) + NoLegend() + 
  geom_hline(yintercept = c(3000, 50000), linetype = "dashed", color = "blue") + NoLegend()

VlnPlot(D7, features = c("percent.mt"), group.by = "orig.ident", pt.size = 0, log = T) +
  geom_hline(yintercept = c(10), linetype = "dashed", color = "blue")


D7 <- subset(D7, subset = nFeature_RNA > 200 & percent.mt < 5)

D7 <- NormalizeData(D7)
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
D7 <- CellCycleScoring(D7, s.features = s.genes, g2m.features = g2m.genes)

D7 <- scutilsR::MarkDoublets(D7, split.by = "orig.ident")

D7 <- subset(D7,subset = DF.classifications=='Singlet')

QuickCluster <- function(object) {
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object, nfeatures = 2000)
  object <- ScaleData(object)
  object <- RunPCA(object)
  object <- FindNeighbors(object, reduction = "pca", dims = 1:30)
  object <- FindClusters(object)
  return(object)
}
seu.list <- SplitObject(D7, split.by = "orig.ident")
seu.list <- lapply(seu.list, QuickCluster)
clusters <- lapply(seu.list, function(xx) xx$seurat_clusters) %>% base::Reduce(c, .)
D7$quick_clusters <- clusters[rownames(D7@meta.data)]
rm(seu.list)
gc()

D7 <- scutilsR::RemoveAmbientRNAs(D7, split.by = "orig.ident", cluster.name = "quick_clusters")
D7
head(D7@meta.data)
DefaultAssay(D7) <- 'decontX'

D7 <- PercentageFeatureSet(D7, pattern = "^MT-", col.name = "percent.mt")
D7 <- PercentageFeatureSet(D7, pattern = "^RPS", col.name = "percent.rps")
D7 <- PercentageFeatureSet(D7, pattern = "^RPL", col.name = "percent.rpl")
D7 <- subset(D7, subset = nFeature_decontX > 200 & percent.mt < 5)
qs::qsave(D7, file = 'output/D7_downstream.qs')



##Integration
V1_V2 <- qs::qread(file = 'output/D7_downstream.qs')
v3 <- qs::qread(file = '../../NTSM_project/NTSM_V4_TA/output/02_output/01_NTSM_integrated_decontX.qs')
V3_D7<-subset(v3, subset = orig.ident == 'NTSM-V4-D7')

Newcellnames <- gsub('NTSM-V4-D7','v3-D7',colnames(V3_D7))
V3_D7 <- RenameCells(V3_D7, new.names = Newcellnames)
V3_D7$orig.ident <- as.factor(V3_D7$orig.ident)
V3_D7$orig.ident <- factor(V3_D7$orig.ident,labels = 'v3-D7')

V1_V2[['RNA']] <- NULL
data <- list(V1_V2, V3_D7)
D7_integration <- merge(x=data[[1]],y=data[-1])
D7_integration$bbknn_res.3= NULL

##
D7_integration <- NormalizeData(D7_integration)
D7_integration <- CellCycleScoring(D7_integration, s.features = s.genes, g2m.features = g2m.genes)
D7_integration <- FindVariableFeatures(D7_integration, selection.method = "vst", nfeatures = 2000)
D7_integration <- ScaleData(D7_integration, vars.to.regress = c('S.Score','G2M.Score'))
D7_integration <- RunPCA(D7_integration, features = VariableFeatures(object = D7_integration),reduction.name = "pca")


#harmony
#integrate with harmony
library(harmony)
D7_integration <- RunHarmony(D7_integration, group.by.vars=c('orig.ident'),reduction.use='pca',dim.use=1:50)
D7_integration <- RunUMAP(D7_integration,reduction = 'harmony',dims = 1:30,n.neighbors = 200,reduction.name = 'umap.harmony')

#bbknn
library(bbknnR)
D7_integration=RunBBKNN(D7_integration,reduction = 'pca',run_UMAP = TRUE, UMAP_name="umap.bbknn",run_TSNE= FALSE,
                        batch_key = 'orig.ident')

p1 =DimPlot(D7_integration,reduction = 'umap.harmony',group.by = 'orig.ident')
p2 =DimPlot(D7_integration,reduction = 'umap.bbknn', group.by = 'orig.ident')
p1+p2

FeaturePlot(D7_integration, features = c('PAX3'), pt.size = 0.2,reduction = "umap.harmony",split.by = 'orig.ident',
            ncol = 9, order = TRUE)

DimPlot(D7_integration,group.by = 'detailed_celltype',split.by = 'orig.ident',reduction = 'umap.bbknn')
ggsave(filename = 'D7_integration_detailed_celltype.tiff',width = 13,height = 5,dpi = 300,path = 'output/Figures/')

resolutions <- c(1.2,1.4,1.6,1.8,2.0,2.4,2.8,3.0)

D7_integration <- FindClusters(D7_integration, resolution = resolutions,graph.name = 'bbknn')
DimPlot(D7_integration,reduction = 'umap.bbknn',group.by =  paste0("bbknn_res.", resolutions), ncol = 3, label = T) & NoLegend()
ggsave(filename = 'D7_integration_bbknn_clusters.tiff',width = 14,height =11,dpi = 300,path = 'output/Figures/')

qs::qsave(D7_integration, file = 'output/D7_integration.qs')
D7_integration <- qs::qread(file = 'output/D7_integration.qs')

gene_list<-unique(GOI)

pdf("output/Figures/D7_integrated_bbknn.pdf", width = 10, height = 5)  # 调整大小以适应每页的内容

for (gene in gene_list) {
  p <- FeaturePlot(
    D7_integration, 
    features = c(gene), 
    pt.size = 0.2, 
    reduction = "umap.bbknn", 
    split.by = 'orig.ident', 
    order = TRUE
  )
  print(p)
}

dev.off()

*********************************************************************
*****************************Day6&7_project_plot*********************
*********************************************************************

********************************************
******************feature_plot**************
********************************************
library(bbknnR)
library(tidyverse)
library(cowplot)
library(RColorBrewer)
D6 <- qs::qread(file = 'output/03_Figure6_Task/D6_annotation.qs')
s <- subset(D6, subset = group %in% c('NTC','VANGL2-KO'))
c <- subset (s,subset = celltype%in% c('PSM','SM','Caud.NP','NT'))


c=RunBBKNN(c,reduction = 'pca',run_UMAP = TRUE, UMAP_name="umap.subset",run_TSNE= FALSE,
           batch_key = 'orig.ident')

D6_Cell <- setNames(
  c("#FF1677", "#FAA31C", "#CE471C", "#D8EA60", "#D58A81", "#00FB0D",
    "#0D8AFF", "#0D8AFF", "#0D8AFF", "#16FEC8", "#16FEC8", "#5335C1",
    "#D194EE", "#D194EE", "#2A40FE", "#2A40FE", "#B6CAFA", "#2E7300",
    "#9A8116", "#FED426"),
  c("Noto", "Caud.NP", "pre-FP", "NT", "FP", "PSM",
    "SM", "SM", "SM", "SCL", "SCL", "END",
    "SYN", "SYN", "MYO", "MYO", "SE", "Int.Meso",
    "EC", "Neuron")
)


DimPlot(c,reduction = 'umap.subset',group.by = 'group')
DimPlot(c,reduction = 'umap.subset',group.by = 'celltype')+scale_color_manual(values = D6_Cell)
ggsave(filename = 'output/05_celltype_subset/01_Dimplot_group.pdf',height =5.2,width = 6.8,dpi = 300 )
ggsave(filename = 'output/05_celltype_subset/01_Dimplot_celltype.pdf',height =5.2,width = 6.8,dpi = 300 )

##feature plot
Markers <- c("TBX6", "MESP2", "DLL1", "DLL3", "LFNG", "HES7", "NOTCH1",
             "ITGA5", "FN1", "PTK2", "PTK7", "DCHS1", "DCHS2",
             "CELSR1", "CELSR2", "CELSR3", "PRICKLE1", "PRICKLE2",
             "DVL2", "DVL3", "FZD3", "FZD6", "TBX18", "UNCX",
             "WNT5A", "WNT5B", "ROR2", "VANGL1",
             "WNT3A", "LEF1", "FGF4", "FGF8", "FGF17",
             "BMP4", "BMP5", "BMP7")
for (marker in Markers) {
  filename <- paste0("output/05_celltype_subset/02_featureplot/", marker, ".tiff")
  tiff(filename = filename, width = 8, height = 3.3, units = "in", res = 300)
  print(FeaturePlot(c,features = marker,reduction = 'umap.subset',pt.size = 0.2,order = T,split.by = 'group'))
  dev.off()
  message("Saved: ", filename)  # 打印保存成功信息
}

tiff(filename = "output/05_celltype_subset/02_featureplot/whole.tiff", 
     width = 3000, height = 3300, units = "px")
FeaturePlot(c, features = unique(Markers), pt.size = 0.2,reduction = "umap.subset",
            ncol = 5, order = TRUE,split.by = 'group')
dev.off()

********************************************
********************Dimplot*****************
********************************************
library(Seurat)
library(ggplot2)
library(bbknnR)
library(tidyverse)
library(cowplot)
library(RColorBrewer)
D6 <- qs::qread(file = 'output/03_Figure6_Task/D6_annotation.qs')
D6_Noto <- subset(D6, subset = group %in% c('NTC', 'NOTO-KO', 'SANT1'))
N <- subset (D6_Noto,subset = celltype%in% c('Noto'))
N=RunBBKNN(N,reduction = 'pca',run_UMAP = TRUE, UMAP_name="umap.bbknn",run_TSNE= FALSE,
           batch_key = 'group')

brewer.pal(n=3,name="Set1")
group_colors <- c("NTC" = "#E41A1C", "NOTO-KO" = "#377EB8", "SANT1" = "#4DAF4A")
Idents(N) <- N$group
DimPlot(N, reduction = "umap.bbknn")+scale_colour_manual(values = group_colors)
ggsave(filename = 'output/03_Figure6_Task/04_Noto_dim.tiff',height =3.9,width = 4.8,dpi = 300 )
DimPlot(N, reduction = "umap.bbknn",split.by = 'group')+scale_colour_manual(values = group_colors)
ggsave(filename = 'output/03_Figure6_Task/04_Noto_split.tiff',height =3.9,width = 10,dpi = 300 )


GOI <- c("SOX2", "TBXT", "FOXA1", "FOXA2", "NKX6-1", 
         "NFIA", "SOX9", "FOXJ1", "RFX2", "MNX1", 
         "TCTEX1D1", "SHH", "NOG", "CHRD", "SEMA3C")

pdf("output/03_Figure6_Task/04_Noto_feature_split.pdf", width = 7, height = 2.5) 
for (gene in GOI) {
  p <- FeaturePlot(
    N, 
    features = c(gene), 
    pt.size = 2, 
    reduction = "umap.bbknn", 
    split.by = 'group', 
    order = TRUE
  )
  print(p)
}

dev.off()


N[['RNA']] <- N[['decontX']]
N2 <- N[rownames(N[['RNA']]@scale.data),]
DefaultAssay(N2) ='RNA'
N2[['decontX']] <- NULL
sceasy::convertFormat(N2, from = "seurat", to = "anndata", main_layer = "scale.data", 
                      outFile = file.path('output/03_Figure6_Task/04_Noto_scaled.h5ad'))



#### MISC ####

root.cells <- CellSelector(FeaturePlot(N, reduction = "umap.bbknn", features = "SEMA3C"))
root.cells <- "v3-d6-SANT1_CACGATTGTCTAGCCA-1"
terminal.p <- c("v3-d6-NOTO_CTAGCGCGTGAGACAC-1")
DimPlot(NTSM, reduction = "umap.bbknn", cells.highlight = root.cells)
DimPlot(NTSM, reduction = "tsne.bbknn", cells.highlight = root.cells)

********************************************
****************Violin_plot*****************
********************************************
library(Seurat)
library(ggplot2)
library(bbknnR)
library(tidyverse)
library(cowplot)
library(RColorBrewer)

D6 <- qs::qread(file = 'output/03_Figure6_Task/D6_annotation.qs')
D6_Noto <- subset(D6, subset = group %in% c('NTC', 'NOTO-KO', 'SANT1'))
N <- subset (D6_Noto,subset = celltype%in% c('Noto'))

genes <- c("SHH", "FOXJ1", "RFX2", "TFF3", "TPPP3", 
           "TCTEX1D1", "C1orf189", "SDR16C5", "VEGFA")
#foxj1 example
foxji <- subset(N, subset = FOXJ1 > 0)
VlnPlot(foxji, features =c("FOXJ1"), split.by = 'group',split.plot = TRUE, cols = c("#E41A1C", "#377EB8", "#4DAF4A"),pt.size = 1)+
  stat_summary(fun = "mean", geom = "crossbar", size = 0.5,color = "black",position = position_dodge(width = 0.9))+
  guides(fill = guide_legend(override.aes = list(color = NA)))
ggsave(filename = 'output/03_Figure6_Task/Vinloplot/foxj1.pdf',width = 4.6,height = 3.2,dpi = 300)


********************************************
********************Heatmap*****************
********************************************

library(Seurat)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(Cairo)
scale_to_range <- function(x, new_min, new_max) {
  old_min <- min(x, na.rm = TRUE) 
  old_max <- max(x, na.rm = TRUE) 
  scaled_x <- ((x - old_min) / (old_max - old_min)) * (new_max - new_min) + new_min
  return(scaled_x)
}

plot_heatmap <- function(gene_list, data, title) {
  data_subset <- data %>% select(celltype, all_of(gene_list))
  rownames(data_subset) <- data_subset$celltype
  data_subset <- data_subset[, -1]
  data_matrix <- as.matrix(data_subset)
  data_scaled <- apply(data_matrix, 2, scale_to_range, new_min = -3, new_max = 3)
  pheatmap(
    data_scaled,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = title,
    fontsize_row = 10,
    fontsize_col = 10,
    show_rownames=TRUE,
    row_dend_side = "left"
  )
}

##markers to plot
GOI_1<- c('SIX1', 'LEF1', 'TCF15', 'IGFBP5', 'MEOX1', 'HES4', 'PAX3', 'TBX18', 'UNCX', 'BMP4')
GOI_2<- c('PAX1', 'PAX9', 'NKX3-2', 'GLI1', 'COL3A1', 'HAS2', 'TWIST1', 'SNAI1', 'FOXC2', 'FRZB')
GOI_3 <- c('SCX', 'LOX', 'COL1A1', 'FMOD', 'TNMD','EBF1', 'EBF2', 'EBF3', 'BMP5', 'TGFBI', 'MEF2C', 'NEFM', 'PDGFRA', 'CXCL12', 'ALCAM')
GOI_4 <-c('PAX7', 'MYOG', 'MEOX2', 'TTN','MYF6', 'SIX4', 'MYF5')
GOI_5 <-c('CD34', 'KDR', 'ETV2', 'SOX17', 'VEGFA', 'CLDN5')
GOI_6 <-c('HOXC5', 'HOXC9', 'HOXA3', 'HOXA5','HOXA7', 'HOXA9', 'HOXB1', 'HOXB2', 'HOXB3','HOXB4', 'HOXB5', 'HOXB8', 'HOXB9','HOXD1', 'HOXD3', 'HOXD4')
rm(list = ls(pattern = "^GOI"))
GOI_EMT <- c('CXCL12', 'CXCR4', 'TGFBR1', 'TGFBI', 'SPRY2', 'BMP4', 'BMP7', 'BMPR1A', 'FST', 'FGF8', 'FGFR1','CTNNB1', 'WNT5A', 'SMAD3', 'SMAD4', 'MSX1', 'NFATC1', 'WWTR1', 'ZEB1', 'FZD10')
GOI_NOTCH <- c('NOTCH1', 'NOTCH2', 'NOTCH3', 'DLL1', 'DLL3','HES1', 'HES4', 'HES7', 'LFNG', 'DLK1', 'WDR12', 'CCND1', 'HEY2', 'HEYL','KRT19', 'ZBTB7A')
GOI_somitogenesis <- c('ALDH2', 'ALDH1A2', 'CDC42BPA','CDC42BPB', 'ITGA5', 'FN1', 'EFNB2','EPHA1', 'EPHA2', 'EPHA4','MESP2', 'LEF1', 'FZD7')
GOI_PCP <- c('FZD3', 'FZD6', 'CELSR1', 'CELSR2', 'CELSR3','PRICKLE1', 'PRICKLE2', 'ROR2', 'PTK7', 'DVL1','DVL2', 'DVL3', 'DCHS1', 'GPC3', 'CTHRC1','SFRP1', 'SFRP2', 'PTK2')

D6 <- qs::qread(file = 'output/03_Figure6_Task/D6_annotation.qs')

####NTC####
NTC <- subset(D6, subset = group == 'NTC')
NTC_type <- subset(NTC, subset = celltype %in% c('SCL','SYN','MYO','END','EC'))
NTC_DF <- FetchData(NTC_type,vars = c('celltype',GOI_1,GOI_2,GOI_3,GOI_4,GOI_5,GOI_6))
NTC_DF$celltype <- as.factor(NTC_DF$celltype)
NTC_DF$celltype <- factor(NTC_DF$celltype,levels = c('SCL','SYN','MYO','END','EC'))

NTC_DF <- NTC_DF %>%
  group_by(celltype) %>%
  summarise(across(where(is.numeric), mean)) %>% 
  as.data.frame()
p1 = plot_heatmap(GOI_1, NTC_DF,'GOI_1 Heatmap')
p2 = plot_heatmap(GOI_2, NTC_DF, "GOI_2 Heatmap")
p3 = plot_heatmap(GOI_3, NTC_DF, "GOI_3 Heatmap")
p4 = plot_heatmap(GOI_4, NTC_DF, "GOI_4 Heatmap")
p5 = plot_heatmap(GOI_5, NTC_DF, "GOI_5 Heatmap")
p6 = plot_heatmap(GOI_6, NTC_DF, "GOI_6 Heatmap")


blank <- grid.rect(gp = gpar(col = NA))
CairoPDF('output/03_Figure6_Task/05_heatmap/NTC_heatmap.pdf',width = 18.8,height = 2.6)
grid.arrange(
  p1$gtable, blank, p2$gtable, blank, p3$gtable,blank,p4$gtable,blank,p5$gtable,blank,p6$gtable,
  ncol = 11,  # 每行 5 个元素
  widths = c(1, 0.1, 1, 0.1, 1, 0.1, 1, 0.1, 1, 0.1, 1)
)
dev.off()


####NOTO_K0#####
NOTO_KO <- subset(D6, subset = group == 'NOTO-KO')
NOTO_type <- subset(NOTO_KO, subset = celltype %in% c('SCL','SYN','MYO','END','EC'))
NOTO_DF <- FetchData(NOTO_type,vars = c('celltype',GOI_1,GOI_2,GOI_3,GOI_4,GOI_5,GOI_6))
NOTO_DF$celltype <- as.factor(NOTO_DF$celltype)
NOTO_DF$celltype <- factor(NOTO_DF$celltype,levels = c('SCL','SYN','MYO','END','EC'))

NOTO_DF <- NOTO_DF %>%
  group_by(celltype) %>%
  summarise(across(where(is.numeric), mean)) %>% 
  as.data.frame()
p1 = plot_heatmap(GOI_1, NOTO_DF,'GOI_1 Heatmap')
p2 = plot_heatmap(GOI_2, NOTO_DF, "GOI_2 Heatmap")
p3 = plot_heatmap(GOI_3, NOTO_DF, "GOI_3 Heatmap")
p4 = plot_heatmap(GOI_4, NOTO_DF, "GOI_4 Heatmap")
p5 = plot_heatmap(GOI_5, NOTO_DF, "GOI_5 Heatmap")
p6 = plot_heatmap(GOI_6, NOTO_DF, "GOI_6 Heatmap")


blank <- grid.rect(gp = gpar(col = NA))
CairoPDF('output/03_Figure6_Task/05_heatmap/NOTO_KO_heatmap.pdf',width = 18.8,height = 2.6)
grid.arrange(
  p1$gtable, blank, p2$gtable, blank, p3$gtable,blank,p4$gtable,blank,p5$gtable,blank,p6$gtable,
  ncol = 11, 
  widths = c(1, 0.1, 1, 0.1, 1, 0.1, 1, 0.1, 1, 0.1, 1) 
)
dev.off()

######SANT1######
SANT1 <- subset(D6, subset = group == 'SANT1')
SANT1_type <- subset(SANT1, subset = celltype %in% c('SCL','SYN','MYO','END','EC'))
SANT1_DF <- FetchData(SANT1_type,vars = c('celltype',GOI_1,GOI_2,GOI_3,GOI_4,GOI_5,GOI_6))
SANT1_DF$celltype <- as.factor(SANT1_DF$celltype)
SANT1_DF$celltype <- factor(SANT1_DF$celltype,levels = c('SCL','SYN','MYO','END','EC'))

SANT1_DF <- SANT1_DF %>%
  group_by(celltype) %>%
  summarise(across(where(is.numeric), mean)) %>% 
  as.data.frame()
p1 = plot_heatmap(GOI_1, SANT1_DF,'GOI_1 Heatmap')
p2 = plot_heatmap(GOI_2, SANT1_DF, "GOI_2 Heatmap")
p3 = plot_heatmap(GOI_3, SANT1_DF, "GOI_3 Heatmap")
p4 = plot_heatmap(GOI_4, SANT1_DF, "GOI_4 Heatmap")
p5 = plot_heatmap(GOI_5, SANT1_DF, "GOI_5 Heatmap")
p6 = plot_heatmap(GOI_6, SANT1_DF, "GOI_6 Heatmap")


blank <- grid.rect(gp = gpar(col = NA))
CairoPDF('output/03_Figure6_Task/05_heatmap/SANT1_heatmap.pdf',width = 18.8,height = 2.6)
grid.arrange(
  p1$gtable, blank, p2$gtable, blank, p3$gtable,blank,p4$gtable,blank,p5$gtable,blank,p6$gtable,
  ncol = 11, 
  widths = c(1, 0.1, 1, 0.1, 1, 0.1, 1, 0.1, 1, 0.1, 1) 
)
dev.off()

#####D6_pertub####

c <- subset(D6, subset = celltype %in% c('PSM','SM','Caud.NP','NT','pre-FP','FP'))

##NTC####
NTC <- subset(D6 ,subset = group == 'NTC')
NTC_type <- subset(NTC, subset = celltype %in% c('PSM','SM','Caud.NP','NT','pre-FP','FP'))
NTC_DF <- FetchData(NTC_type, vars = c('celltype',GOI_EMT,GOI_NOTCH,GOI_somitogenesis,GOI_PCP))
NTC_DF$celltype <- as.factor(NTC_DF$celltype)
NTC_DF$celltype <- factor(NTC_DF$celltype,levels =c('PSM','SM','Caud.NP','NT','pre-FP','FP') )

NTC_DF <- NTC_DF %>%
  group_by(celltype) %>%
  summarise(across(where(is.numeric), mean)) %>% 
  as.data.frame()
p1 = plot_heatmap(GOI_EMT, NTC_DF,'EMT')
p2 = plot_heatmap(GOI_NOTCH, NTC_DF, "NOTCH")
p3 = plot_heatmap(GOI_somitogenesis, NTC_DF, "Somitogenesis")
p4 = plot_heatmap(GOI_PCP, NTC_DF, "PCP")

blank <- grid.rect(gp = gpar(col = NA))
CairoPDF('output/03_Figure6_Task/05_heatmap/D6_NTC.pdf',width = 23,height = 3)
grid.arrange(
  p1$gtable, blank, p2$gtable, blank, p3$gtable,blank,p4$gtable,
  ncol = 7,
  widths = c(1, 0.1, 1, 0.1, 1, 0.1, 1)
)
dev.off()

####VANGLE-KO
V <- subset(D6,subset = group == 'VANGL2-KO')
V_type <- subset(V, subset = celltype %in% c('PSM','SM','Caud.NP','NT','pre-FP','FP'))
V_DF <- FetchData(V_type, vars = c('celltype',GOI_EMT,GOI_NOTCH,GOI_somitogenesis,GOI_PCP))
V_DF$celltype <- as.factor(V_DF$celltype)
V_DF$celltype <- factor(V_DF$celltype,levels =c('PSM','SM','Caud.NP','NT','pre-FP','FP') )

V_DF <- V_DF %>%
  group_by(celltype) %>%
  summarise(across(where(is.numeric), mean)) %>% 
  as.data.frame()
p1 = plot_heatmap(GOI_EMT, V_DF,'EMT')
p2 = plot_heatmap(GOI_NOTCH, V_DF, "NOTCH")
p3 = plot_heatmap(GOI_somitogenesis, V_DF, "Somitogenesis")
p4 = plot_heatmap(GOI_PCP, V_DF, "PCP")

blank <- grid.rect(gp = gpar(col = NA))
CairoPDF('output/03_Figure6_Task/05_heatmap/D6_VANGL2-KO.pdf',width = 23,height = 3)
grid.arrange(
  p1$gtable, blank, p2$gtable, blank, p3$gtable,blank,p4$gtable,
  ncol = 7,  
  widths = c(1, 0.1, 1, 0.1, 1, 0.1, 1) 
)
dev.off()


********************************************
*****************other_plot*****************
********************************************

# r-seurat4.3
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(reshape2)
library(future)
library(ComplexHeatmap)
library(qs)
library(presto)
library(circlize)
library(ggrepel)

###############
#### DAY 6 ####
###############
srt <- qread("data/D6_integration.qs")
srt <- FindClusters(srt, resolution = 1.2, graph.name = 'bbknn')

cluster_to_celltype <- c(
  "14" = "Noto",
  "1" = "NMP-neural",
  "18" = "Caud. NP",
  "6" = "NT",
  "15" = "FP",
  "11" = "PSM",
  "2" = "SM",
  "5" = "SM",
  "9" = "SM",
  "10" = "SM",
  "4" = "SCL",
  "13" = "SCL",
  "8" = "Endotome",
  "3" = "Syndetome",
  "17" = "Dermomyotome",
  "0" = "Myogenic",
  "16" = "Caud. Meso",
  "12" = "Int. Meso",
  "7" = "EC",
  "19" = "Neuron"
)

cluster_to_detailed_celltype <- c(
  "14" = "Noto",
  "1" = "NMP-neural",
  "18" = "Caud. NP",
  "6" = "NT",
  "15" = "FP",
  "11" = "PSM",
  "2" = "SM.1",
  "5" = "SM.2",
  "9" = "SM.3",
  "10" = "SM.4",
  "4" = "SCL.1",
  "13" = "SCL.2",
  "8" = "Endotome",
  "3" = "Syndetome",
  "17" = "Dermomyotome",
  "0" = "Myogenic",
  "16" = "Caud. Meso",
  "12" = "Int. Meso",
  "7" = "EC",
  "19" = "Neurons"
)

detailed_celltype_to_color <- c(
  "Noto" = "#FF1677",
  "NMP-neural" = "#FAA31C",
  "Caud. NP" = "#CE471C",
  "NT" = "#D8EA60",
  "FP" = "#D58A81",
  "PSM" = "#00FB0D",
  "SM.1" = "#0D8AFF",
  "SM.2" = "#0D8AFF",
  "SM.3" = "#0D8AFF",
  "SM.4" = "#0D8AFF",
  "SCL.1" = "#16FEC8",
  "SCL.2" = "#16FEC8",
  "Endotome" = "#5335C1",
  "Syndetome" = "#D194EE",
  "Dermomyotome" = "#F96ADB",
  "Myogenic" = "#2A40FE",
  "Caud. Meso" = "#B6CAFA",
  "Int. Meso" = "#2E7300",
  "EC" = "#9A8116",
  "Neurons" = "#FED426"
)

celltype_to_color <- c(
  "Noto" = "#FF1677",
  "NMP-neural" = "#FAA31C",
  "Caud. NP" = "#CE471C",
  "NT" = "#D8EA60",
  "FP" = "#D58A81",
  "PSM" = "#00FB0D",
  "SM" = "#0D8AFF",
  "SCL" = "#16FEC8",
  "Endotome" = "#5335C1",
  "Syndetome" = "#D194EE",
  "Dermomyotome" = "#F96ADB",
  "Myogenic" = "#2A40FE",
  "Caud. Meso" = "#B6CAFA",
  "Int. Meso" = "#2E7300",
  "EC" = "#9A8116",
  "Neuron" = "#FED426"
)

cluster_to_color <- c(
  "14" = "#FF1677",
  "1" = "#FAA31C",
  "18" = "#CE471C",
  "6" = "#D8EA60",
  "15" = "#D58A81",
  "11" = "#00FB0D",
  "2" = "#0D8AFF",
  "5" = "#0D8AFF",
  "9" = "#0D8AFF",
  "10" = "#0D8AFF",
  "4" = "#16FEC8",
  "13" = "#16FEC8",
  "8" = "#5335C1",
  "3" = "#D194EE",
  "17" = "#F96ADB",
  "0" = "#2A40FE",
  "16" = "#B6CAFA",
  "12" = "#2E7300",
  "7" = "#9A8116",
  "19" = "#FED426"
)

cell_order <- c(
  "Noto",
  "NMP-neural",
  "Caud. NP",
  "NT",
  "FP",
  "PSM",
  "SM",
  "SCL",
  "Endotome",
  "Syndetome",
  "Dermomyotome",
  "Myogenic",
  "Caud. Meso",
  "Int. Meso",
  "EC",
  "Neuron"
)

detailed_celltype_order <- c(
  "Noto",
  "NMP-neural",
  "Caud. NP",
  "NT",
  "FP",
  "PSM",
  "SM.1",
  "SM.2",
  "SM.3",
  "SM.4",
  "SCL.1",
  "SCL.2",
  "Endotome",
  "Syndetome",
  "Dermomyotome",
  "Myogenic",
  "Caud. Meso",
  "Int. Meso",
  "EC",
  "Neurons"
)

fig6_to_celltype <- c(
  "14" = "Noto",
  "1" = "Caud. NP",
  "18" = "Pre-FP",
  "6" = "NT",
  "15" = "FP",
  "11" = "PSM",
  "9" = "SM",
  "10" = "SM",
  "17" = "SM",
  "4" = "SCL",
  "13" = "SCL",
  "8" = "END",
  "2" = "SYN",
  "3" = "SYN",
  "5" = "MYO",
  "0" = "MYO",
  "16" = "SE",
  "12" = "Int. Meso",
  "7" = "EC",
  "19" = "Neuron"
)

fig6_colors <- c(
  "Noto" = "#FF1677",
  "Caud. NP" = "#FAA31C",
  "Pre-FP" = "#CE471C",
  "NT" = "#D8EA60",
  "FP" = "#D58A81",
  "PSM" = "#00FB0D",
  "SM" = "#0D8AFF",
  "SM_duplicate1" = "#0D8AFF",
  "SM_duplicate2" = "#F96ADB",
  "SCL" = "#0D8AFF",
  "SCL_alt" = "#16FEC8",
  "END" = "#5335C1",
  "SYN" = "#D194EE",
  "SYN_duplicate" = "#D194EE",
  "MYO" = "#2A40FE",
  "MYO_duplicate" = "#2A40FE",
  "SE" = "#B6CAFA",
  "Int. Meso" = "#2E7300",
  "EC" = "#9A8116",
  "Neuron" = "#FED426"
)

srt@meta.data$celltype <- cluster_to_celltype[as.character(srt@meta.data$seurat_clusters)]
srt@meta.data$detailed_celltype <- cluster_to_detailed_celltype[as.character(srt@meta.data$seurat_clusters)]
srt@meta.data$fig6_celltype <- fig6_to_celltype[as.character(srt@meta.data$seurat_clusters)]

srt@meta.data$celltype <- factor(srt@meta.data$celltype, levels = cell_order)
srt@meta.data$detailed_celltype <- factor(srt@meta.data$detailed_celltype, levels = detailed_celltype_order)

options(repr.plot.width = 8, repr.plot.height = 8)
tiff('fig_6_1_dimplot.tiff', width = 8, height = 8, res = 300, units = 'in')
DimPlot(srt, reduction = 'umap.bbknn', label = T, repel = T, group.by = 'seurat_clusters', cols = cluster_to_color)
dev.off()

tiff('fig_6_1_dimplot_celltype.tiff', width = 9, height = 8, res = 300, units = 'in')
DimPlot(srt, reduction = 'umap.bbknn', group.by = 'celltype', cols = celltype_to_color)
dev.off()

tiff('fig_6_1_dimplot_detailed_celltype.tiff', width = 9, height = 8, res = 300, units = 'in')
DimPlot(srt, reduction = 'umap.bbknn', group.by = 'detailed_celltype', cols = detailed_celltype_to_color)
dev.off()

#### Cell type bar ####
table(srt@meta.data$orig.ident, srt@meta.data$fig6_celltype) %>% data.frame() -> cellcount_df
colnames(cellcount_df) <- c('group', 'celltype', 'count')
# write.csv(cellcount_df, 'fig_6_table_6_4_1_cell_count.csv')

table(srt@meta.data$orig.ident, srt@meta.data$fig6_celltype) -> celltype_prop
celltype_prop <- celltype_prop/rowSums(celltype_prop)
df <- data.frame(celltype_prop)
colnames(df) <- c('group', 'celltype', 'percent')
# write.csv(df, 'fig_6_table_6_4_2_cell_percent.csv')

options(repr.plot.width = 4, repr.plot.height = 8)
tiff('fig_6_3_cellcount_bar.tiff', width = 4, height = 8, res = 300, units = 'in')
print(df %>% ggplot(aes(x = group, y = percent)) + 
    geom_col(aes(fill = celltype), colour="black") + 
    theme_minimal() + 
    ylim(0, 1) + 
    scale_y_continuous(
        limits = c(0, 1.01),  # Set the y-axis range
        expand = c(0, 0)     # No padding
    ) +
    theme(panel.border = element_blank(),
        panel.background = element_blank(),  # Remove background
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),  # Remove minor gridlines
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_manual(values = celltype_to_color) + 
    labs(y = 'Percent', x = 'Group', fill = 'Celltype'))
dev.off()

#### Cell type ratio
df_ntc <- df %>% filter(group == 'v3_d6') %>% arrange(celltype)
df_noto <- df %>% filter(group == 'v3-d6-NOTO') %>% arrange(celltype)
df_log_ratio <- data.frame(celltype = levels(df$celltype), ratio = log2(df_noto$percent/df_ntc$percent))
df_log_ratio$celltype <- factor(df_log_ratio$celltype, levels = rev(levels(df$celltype)))

options(repr.plot.width = 6, repr.plot.height = 2.75)

cell_use <- c('Noto', 'Caud. NP','NT','FP','PSM','SM','SCL','SYN','MYO','END','EC')

df_log_ratio %>% filter(celltype %in% cell_use) -> plot_df
plot_df$celltype <- factor(plot_df$celltype, levels = rev(cell_use))
plot_df %>%    ggplot(aes(x = ratio, y = celltype)) + 
    geom_rect(aes(xmin=-1.5, xmax=1.5, ymin=0, ymax=Inf), color='#E5F4FB', fill='#E5F4FB') +
    geom_point(size = 4, aes(fill = celltype), stroke=0.5, colour="black", shape = 21) + 
    geom_vline(xintercept = 0) +
    theme(panel.border = element_blank(),
        panel.background = element_blank(),  # Remove background
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),  # Remove minor gridlines
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = 'na'
        ) + scale_fill_manual(values = fig6_colors) + 
        labs(x = 'log2(NOTO-KO-percent/NTC-percent)', y = 'Celltype') + 
        xlim(-3.9, 3) -> p1

tiff('fig_6B_4_ratio_main_NOTO_vs_NTC.tiff', width = 6, height = 2.75, res = 300, units = 'in')
print(p1)
dev.off()

df_ntc <- df %>% filter(group == 'v3_d6') %>% arrange(celltype)
df_sant <- df %>% filter(group == 'v3-d6-SANT1') %>% arrange(celltype)
df_log_ratio <- data.frame(celltype = levels(df$celltype), ratio = log2(df_sant$percent/df_ntc$percent))
df_log_ratio$celltype <- factor(df_log_ratio$celltype, levels = rev(levels(df$celltype)))

cell_use <- c('Noto', 'Caud. NP','NT','FP','PSM','SM','SCL','SYN','MYO','END','EC')

df_log_ratio %>% filter(celltype %in% cell_use) -> plot_df
plot_df$celltype <- factor(plot_df$celltype, levels = rev(cell_use))
plot_df %>%
    ggplot(aes(x = ratio, y = celltype)) + 
    geom_rect(aes(xmin=-1.5, xmax=1.5, ymin=0, ymax=Inf), color='#E5F4FB', fill='#E5F4FB') +
    geom_point(size = 4, aes(fill = celltype), stroke=0.5, colour="black", shape = 21) + 
    geom_vline(xintercept = 0) +
    theme(panel.border = element_blank(),
        panel.background = element_blank(),  # Remove background
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),  # Remove minor gridlines
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = 'na'
        ) + scale_fill_manual(values = fig6_colors) + 
        labs(x = 'log2(SANT1-percent/NTC-percent)', y = 'Celltype') + 
        xlim(-3.9, 3) -> p1
p1
tiff('fig_6B_4_ratio_main_SANT1_vs_NTC.tiff', width = 6, height = 2.75, res = 300, units = 'in')
print(p1)
dev.off()

df_ntc <- df %>% filter(group == 'v3_d6') %>% arrange(celltype)
df_sant <- df %>% filter(group == 'v3-d6-VANGLE2') %>% arrange(celltype)
df_log_ratio <- data.frame(celltype = levels(df$celltype), ratio = log2(df_sant$percent/df_ntc$percent))
df_log_ratio$celltype <- factor(df_log_ratio$celltype, levels = rev(levels(df$celltype)))
options(repr.plot.width = 6, repr.plot.height = 2.5)

cell_use <- c('Noto', 'Caud. NP','NT','FP','PSM','SM','SCL','SYN','MYO','END','EC')

df_log_ratio %>% filter(celltype %in% cell_use) -> plot_df
plot_df$celltype <- factor(plot_df$celltype, levels = rev(cell_use))
plot_df %>%
    ggplot(aes(x = ratio, y = celltype)) + 
    geom_rect(aes(xmin=-1.5, xmax=1.5, ymin=0, ymax=Inf), color='#E5F4FB', fill='#E5F4FB') +
    geom_point(size = 4, aes(fill = celltype), stroke=0.5, colour="black", shape = 21) + 
    geom_vline(xintercept = 0) +
    theme(panel.border = element_blank(),
        panel.background = element_blank(),  # Remove background
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),  # Remove minor gridlines
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = 'na'
        ) + scale_fill_manual(values = fig6_colors) + 
        labs(x = 'log2(VANGL2-percent/NTC-percent)', y = 'Celltype') + 
        xlim(-3.9, 3) -> p1
tiff('fig_6B_4_ratio_main_VANGLE2_vs_NTC.tiff', width = 6, height = 2.75, res = 300, units = 'in')
print(p1)
dev.off()

#### Heatmap ####
GOI <- c(
  "SOX2", "TBXT", "CHRD", "SHH", "GLI1", "NOG", "SEMA3C",
  "PAX6", "MSX1", "DBX2", "IRX3", "EGF", "NKX6-2", "NKX6-1", "OLIG2", "OLIG1", "FOXA1", "FOXA2",
  "FGF3", "FGF4", "FGF8", "FGF17", "FGF18", "FGFR1", "FGFR2", "FGFR3",
  "HES1", "HES4", "HES5", "HES7", "DLL1", "DLL3", "LFNG", "NOTCH1", "NOTCH2", "NOTCH3"
)
length(GOI)
GOI_group <- c(rep("1", 7), 
    rep('2', 11), 
    rep('3', 8), rep('4', 10))

cell_order_heatmap <- c('Noto', 'Caud. NP', 'NT', 'Pre-FP', 'FP', 'PSM')

srt@meta.data$condition <- paste0(srt@meta.data$fig6_celltype, "_", srt@meta.data$orig.ident)

Idents(srt) <- 'condition'
mat_df <- data.frame(index = 1:length(unique(GOI)))

s1 <- 'v3_d6'
s2 <- 'v3-d6-NOTO'

for (c in cell_order_heatmap){
    id1 <- paste0(c, '_', s1)
    id2 <- paste0(c, '_', s2)
    log2fc_results <- FindMarkers(
        object = srt, # Your Seurat object
        ident.1 = id2,   # First cluster or condition
        ident.2 = id1,   # Second cluster or condition
        features = unique(GOI), # Specific genes to analyze
        logfc.threshold = 0,
        min.pct = 0,
        only.pos = F,
        group_by = 'condition'
    )
    tdf <- log2fc_results %>% dplyr::select(avg_log2FC)
    tdf['gene'] <- rownames(tdf)
    tdf <- tdf %>% arrange(gene) %>% dplyr::select(avg_log2FC)
    colnames(tdf) <- c(paste0('NOTO-KO_vs_NTC_', c))
    mat_df <- cbind(mat_df, tdf)
}

s1 <- 'v3_d6'
s2 <- 'v3-d6-SANT1'

for (c in cell_order_heatmap){
    id1 <- paste0(c, '_', s1)
    id2 <- paste0(c, '_', s2)
    log2fc_results <- FindMarkers(
        object = srt, # Your Seurat object
        ident.1 = id2,   # First cluster or condition
        ident.2 = id1,   # Second cluster or condition
        features = unique(GOI), # Specific genes to analyze
        logfc.threshold = 0,
        min.pct = 0,
        only.pos = F,
        group_by = 'condition'
    )
    tdf <- log2fc_results %>% dplyr::select(avg_log2FC)
    tdf['gene'] <- rownames(tdf)
    tdf <- tdf %>% arrange(gene) %>% dplyr::select(avg_log2FC)
    colnames(tdf) <- c(paste0('SANT1_vs_NTC_', c))
    mat_df <- cbind(mat_df, tdf)
}

s1 <- 'v3-d6-SANT1'
s2 <- 'v3-d6-NOTO'

for (c in cell_order_heatmap){
    id1 <- paste0(c, '_', s1)
    id2 <- paste0(c, '_', s2)
    log2fc_results <- FindMarkers(
        object = srt, # Your Seurat object
        ident.1 = id2,   # First cluster or condition
        ident.2 = id1,   # Second cluster or condition
        features = unique(GOI), # Specific genes to analyze
        logfc.threshold = 0,
        min.pct = 0,
        only.pos = F,
        group_by = 'condition'
    )
    tdf <- log2fc_results %>% dplyr::select(avg_log2FC)
    tdf['gene'] <- rownames(tdf)
    tdf <- tdf %>% arrange(gene) %>% dplyr::select(avg_log2FC)
    colnames(tdf) <- c(paste0('NOTO_vs_SANT1_', c))
    mat_df <- cbind(mat_df, tdf)
}

mat_df <- mat_df %>% dplyr::select(-index) 

mat_df %>% as.matrix() -> mat
mat <- mat[GOI,]

plot_mat <- mat
# colnames(plot_mat) <- rep(cell_order_heatmap, 2)
row_split_order <- c(rep('NOTO-KO vs NTC', length(cell_order_heatmap)), rep('SANT1 vs NTC', length(cell_order_heatmap)), rep('NOTO-KO vs SANT1', length(cell_order_heatmap)))
row_split_order <- factor(row_split_order, levels = c('NOTO-KO vs NTC', 'SANT1 vs NTC', 'NOTO-KO vs SANT1'))
col_fun = colorRamp2(c(-0.58, 0, 0.58), c("blue", "white", "red"))
tiff('fig_6E_5_NOTO_SANT1_heatmap.tiff', width = 12, height = 7, res=300, units = 'in')
print(Heatmap(t(plot_mat), 
    name = 'Log2FC', 
    col = col_fun,
    row_split = row_split_order,
    row_order = colnames(plot_mat),
    column_order = GOI,
    column_split = GOI_group, 
    border = TRUE, 
    row_dend_reorder = F,
    column_dend_reorder = F,
    heatmap_legend_param = list(
        at = c(-0.58, 0, 0.58),          # Specify tick marks on the legend
        labels = c("-0.58", "0", "0.58") # Corresponding labels
    ),
    row_labels = rep(cell_order_heatmap, 3)
    ))
dev.off()

###############
#### DAY 7 ####
###############
srt <- qread("data/D7_integration.qs")
srt <- FindClusters(srt, resolution = 1.2, graph.name = 'bbknn')

cluster_to_celltype <- c(
  "17" = "Noto",
  "1" = "NMP-neural",
  "5" = "Caud. NP",
  "3" = "NT",
  "13" = "NT",
  "11" = "FP",
  "6" = "PSM",
  "0" = "SM",
  "8" = "SM",
  "10" = "SM",
  "14" = "SM",
  "4" = "SCL",
  "7" = "Endotome",
  "2" = "Syndetome",
  "12" = "Dermomyotome/Myogenic",
  "9" = "EC",
  "15" = "Neuron",
  "16" = "SE"
)

cluster_to_detailed_celltype <- c(
  "17" = "Notochord",
  "1" = "NMP-neural",
  "5" = "Caud. NP",
  "3" = "Neural Tube.1",
  "13" = "Neural Tube.2",
  "11" = "Floor Plate",
  "6" = "PSM",
  "0" = "Somitic Meso.1",
  "8" = "Somitic Meso.2",
  "10" = "Somitic Meso.3",
  "14" = "Somitic Meso.4",
  "4" = "Sclerotome",
  "7" = "Endotome",
  "2" = "Syndetome",
  "12" = "Dermomyotome/Myogenic",
  "9" = "Endothelial Cells",
  "15" = "Neurons",
  "16" = "Surface Ectoderm"
)

celltype_to_color <- c(
  "Noto" = "#FF1677",
  "NMP-neural" = "#FAA31C",
  "Caud. NP" = "#CE471C",
  "NT" = "#D8EA60",
  "FP" = "#D58A81",
  "PSM" = "#00FB0D",
  "SM" = "#0D8AFF",
  "SCL" = "#16FEC8",
  "Endotome" = "#5335C1",
  "Syndetome" = "#D194EE",
  "Dermomyotome/Myogenic" = "#F96ADB",
  "EC" = "#9A8116",
  "Neuron" = "#FED426",
  "SE" = "#B6CAFA"
)

detailed_celltype_to_color <- c(
  "Notochord" = "#FF1677",
  "NMP-neural" = "#FAA31C",
  "Caud. NP" = "#CE471C",
  "Neural Tube.1" = "#D8EA60",
  "Neural Tube.2" = "#D8EA60",
  "Floor Plate" = "#D58A81",
  "PSM" = "#00FB0D",
  "Somitic Meso.1" = "#0D8AFF",
  "Somitic Meso.2" = "#0D8AFF",
  "Somitic Meso.3" = "#0D8AFF",
  "Somitic Meso.4" = "#0D8AFF",
  "Sclerotome" = "#16FEC8",
  "Endotome" = "#5335C1",
  "Syndetome" = "#D194EE",
  "Dermomyotome/Myogenic" = "#F96ADB",
  "Endothelial Cells" = "#9A8116",
  "Neurons" = "#FED426",
  "Surface Ectoderm" = "#B6CAFA"
)

cluster_to_color <- c(
  "17" = "#FF1677",
  "1" = "#FAA31C",
  "5" = "#CE471C",
  "3" = "#D8EA60",
  "13" = "#D8EA60",
  "11" = "#D58A81",
  "6" = "#00FB0D",
  "0" = "#0D8AFF",
  "8" = "#0D8AFF",
  "10" = "#0D8AFF",
  "14" = "#0D8AFF",
  "4" = "#16FEC8",
  "7" = "#5335C1",
  "2" = "#D194EE",
  "12" = "#F96ADB",
  "9" = "#9A8116",
  "15" = "#FED426",
  "16" = "#B6CAFA"
)

cell_order <- c(
  "Noto",
  "NMP-neural",
  "Caud. NP",
  "NT",
  "FP",
  "PSM",
  "SM",
  "SCL",
  "Endotome",
  "Syndetome",
  "Dermomyotome/Myogenic",
  "EC",
  "Neuron",
  "SE"
)

detailed_celltype_order <- c(
  "Notochord",
  "NMP-neural",
  "Caud. NP",
  "Neural Tube.1",
  "Neural Tube.2",
  "Floor Plate",
  "PSM",
  "Somitic Meso.1",
  "Somitic Meso.2",
  "Somitic Meso.3",
  "Somitic Meso.4",
  "Sclerotome",
  "Endotome",
  "Syndetome",
  "Dermomyotome/Myogenic",
  "Endothelial Cells",
  "Neurons",
  "Surface Ectoderm"
)

srt@meta.data$celltype <- cluster_to_celltype[as.character(srt@meta.data$seurat_clusters)]
srt@meta.data$detailed_celltype <- cluster_to_detailed_celltype[as.character(srt@meta.data$seurat_clusters)]

srt@meta.data$celltype <- factor(srt@meta.data$celltype, levels = cell_order)
srt@meta.data$detailed_celltype <- factor(srt@meta.data$detailed_celltype, levels = detailed_celltype_order)

options(repr.plot.width = 8, repr.plot.height = 8)
tiff('fig_6_1_dimplot.tiff', width = 8, height = 8, res = 300, units = 'in')
DimPlot(srt, reduction = 'umap.bbknn', label = T, repel = T, group.by = 'seurat_clusters', cols = cluster_to_color)
dev.off()

tiff('fig_6_1_dimplot_celltype.tiff', width = 9, height = 8, res = 300, units = 'in')
DimPlot(srt, reduction = 'umap.bbknn', group.by = 'celltype', cols = celltype_to_color)
dev.off()

tiff('fig_6_1_dimplot_detailed_celltype.tiff', width = 9, height = 8, res = 300, units = 'in')
DimPlot(srt, reduction = 'umap.bbknn', group.by = 'detailed_celltype', cols = detailed_celltype_to_color)
dev.off()

# Stacked bar celltype
table(srt@meta.data$orig.ident, srt@meta.data$celltype) %>% data.frame() -> cellcount_df
colnames(cellcount_df) <- c('group', 'celltype', 'count')
write.csv(cellcount_df, 'table_6_3_1_cell_count.csv')

table(srt@meta.data$orig.ident, srt@meta.data$celltype) -> celltype_prop
celltype_prop <- celltype_prop/rowSums(celltype_prop)
df <- data.frame(celltype_prop)
colnames(df) <- c('group', 'celltype', 'percent')
write.csv(df, 'table_6_3_2_cell_percent.csv')

options(repr.plot.width = 4, repr.plot.height = 8)
tiff('supp_fig_6_3_cellcount_bar.tiff', width = 4, height = 8, res = 300, units = 'in')
print(df %>% ggplot(aes(x = group, y = percent)) + 
    geom_col(aes(fill = celltype), colour="black") + 
    theme_minimal() + 
    ylim(0, 1) + 
    scale_y_continuous(
        limits = c(0, 1.01),  # Set the y-axis range
        expand = c(0, 0)     # No padding
    ) +
    theme(panel.border = element_blank(),
        panel.background = element_blank(),  # Remove background
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),  # Remove minor gridlines
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_manual(values = celltype_to_color) + 
    labs(y = 'Percent', x = 'Group', fill = 'Celltype'))
dev.off()

###################################
#### Cell proportion pie chart ####
###################################
celltype_pie <- c('Noto','FP','SM','SCL','Dermomyotome/Myogenic','Endotome')
celltype_pie <- c('Noto', '','FP','SM','SCL','Dermomyotome/Myogenic','Endotome')
celltype_index <- c(1, 2, 3, 4, 5, 6, 7, 8)

color_list <- c( "v1-D7" = "#5335C1",
  "v2-D7" = "#FF1677",
  "v3-D7" = "#FAA31C")

meta_df <- srt@meta.data
table(srt@meta.data$orig.ident, srt@meta.data$celltype) %>% data.frame() -> cellcount_df
colnames(cellcount_df) <- c('group', 'celltype', 'count')

table(srt@meta.data$orig.ident, srt@meta.data$celltype) -> celltype_prop
celltype_prop <- celltype_prop/rowSums(celltype_prop)
df <- data.frame(celltype_prop)
colnames(df) <- c('group', 'celltype', 'percent')

plot_list <- c()

c = 'Noto(SHH+/CHRD+)'

expr_data <- FetchData(object = srt, vars = c('SHH', 'CHRD'), layer = 'data', assay = 'decontX')
expr_data <- cbind(meta_df, expr_data[rownames(meta_df),])

expr_data['expr_type_1'] <- 'NA'
expr_data['expr_type_2'] <- 'NA'
expr_data['expr_type_1'] <- ifelse(expr_data$SHH > 0, 'SHH+', 'SHH-')
expr_data['expr_type_2'] <- ifelse(expr_data['CHRD'] > 0, 'CHRD+', 'CHRD-')
expr_data['expr_type'] <- paste0(expr_data$expr_type_1, "/", expr_data$expr_type_2)
expr_data['expr_type'] <- paste0(expr_data$celltype,"(",expr_data$expr_type,")")

table(expr_data$orig.ident, expr_data$expr_type) %>% data.frame() -> cellcount_df
colnames(cellcount_df) <- c('group', 'celltype', 'count')
table(expr_data$orig.ident) %>% data.frame() -> sample_cellcount_df
rownames(sample_cellcount_df) <- sample_cellcount_df$Var1
cellcount_df['sample_number'] <- sample_cellcount_df[cellcount_df$group, 'Freq']
cellcount_df['percent_per_sample'] = cellcount_df$count*100/cellcount_df$sample_number
cellcount_df1 <- cellcount_df %>% filter(celltype == 'NT(PAX6+/NKX6-1-)')
cellcount_df1['percent'] <- cellcount_df1$percent_per_sample/sum(cellcount_df1$percent_per_sample)
cellcount_df2 <- cellcount_df %>% filter(celltype == 'NT(PAX6-/NKX6-1+)')
cellcount_df2['percent'] <- cellcount_df2$percent_per_sample/sum(cellcount_df2$percent_per_sample)
cellcount_df <- rbind(cellcount_df1, cellcount_df2)
write.csv(cellcount_df, 'revised_table_pie_1_shh_chrd_count.csv')
table(expr_data$orig.ident, expr_data$expr_type) -> celltype_prop
celltype_prop <- celltype_prop/rowSums(celltype_prop)
df <- data.frame(celltype_prop)
colnames(df) <- c('group', 'celltype', 'percent')
# write.csv(df, 'revised_table_pie_1_shh_percent.csv')

df %>% 
  filter(celltype == 'Noto(SHH+/CHRD+)') %>%
  mutate(percent_ratio = percent/sum(percent)) %>%
  mutate(label_position = cumsum(percent_ratio) - percent_ratio/2) %>%
  mutate(group = factor(group, levels = c('v3-D7', 'v2-D7', 'v1-D7'))) %>%
  ggplot(aes(x="", y=percent_ratio, fill=group)) +
    geom_bar(stat="identity", width=1, color='white') +
    coord_polar("y", start=0) +
    scale_fill_manual(values = color_list) + 
    theme_void() + 
    labs(fill = 'group', title = 'Noto(SHH+/CHRD+)') + 
    theme(legend.position = 'na', plot.title = element_text(hjust = 0.5)) +
      geom_text_repel(aes(y = label_position, label = paste0(round(percent*100, 2), "%")), 
                vjust = 0.5, hjust = 0.5, size = 4, color = "white", box.padding = 0.01) -> plot_list[[1]]

expr_data <- FetchData(object = srt, vars = c('PAX6', 'NKX6-1'), layer = 'data', assay = 'decontX')
expr_data <- cbind(meta_df, expr_data[rownames(meta_df),])

expr_data['expr_type_1'] <- 'NA'
expr_data['expr_type_2'] <- 'NA'
expr_data['expr_type_1'] <- ifelse(expr_data$PAX6 > 0, 'PAX6+', 'PAX6-')
expr_data['expr_type_2'] <- ifelse(expr_data['NKX6-1'] > 0, 'NKX6-1+', 'NKX6-1-')
expr_data['expr_type'] <- paste0(expr_data$expr_type_1, "/", expr_data$expr_type_2)
expr_data['expr_type'] <- paste0(expr_data$celltype,"(",expr_data$expr_type,")")

table(expr_data$orig.ident, expr_data$expr_type) %>% data.frame() -> cellcount_df
colnames(cellcount_df) <- c('group', 'celltype', 'count')
table(expr_data$orig.ident) %>% data.frame() -> sample_cellcount_df
rownames(sample_cellcount_df) <- sample_cellcount_df$Var1
cellcount_df['sample_number'] <- sample_cellcount_df[cellcount_df$group, 'Freq']
cellcount_df['percent_per_sample'] = cellcount_df$count*100/cellcount_df$sample_number
cellcount_df1 <- cellcount_df %>% filter(celltype == 'NT(PAX6+/NKX6-1-)')
cellcount_df1['percent'] <- cellcount_df1$percent_per_sample/sum(cellcount_df1$percent_per_sample)
cellcount_df2 <- cellcount_df %>% filter(celltype == 'NT(PAX6-/NKX6-1+)')
cellcount_df2['percent'] <- cellcount_df2$percent_per_sample/sum(cellcount_df2$percent_per_sample)
cellcount_df <- rbind(cellcount_df1, cellcount_df2)
write.csv(cellcount_df, 'revised_table_pie_2_pax6_nkx6-1_count.csv')
table(expr_data$orig.ident, expr_data$expr_type) -> celltype_prop
celltype_prop <- celltype_prop/rowSums(celltype_prop)
df <- data.frame(celltype_prop)
colnames(df) <- c('group', 'celltype', 'percent')
# write.csv(df, 'revised_table_pie_2_pax6_nkx6-1_percent.csv')

df %>% 
  filter(celltype == 'NT(PAX6+/NKX6-1-)') %>%
  mutate(percent_ratio = percent/sum(percent)) %>%
  mutate(label_position = cumsum(percent_ratio) - percent_ratio/2) %>%
  mutate(group = factor(group, levels = c('v3-D7', 'v2-D7', 'v1-D7'))) %>%
  ggplot(aes(x="", y=percent_ratio, fill=group)) +
    geom_bar(stat="identity", width=1, color='white') +
    coord_polar("y", start=0) +
    scale_fill_manual(values = color_list) + 
    theme_void() + 
    labs(fill = 'group', title = 'NT(PAX6+/NKX6-1-)') + 
    theme(legend.position = 'na', plot.title = element_text(hjust = 0.5)) +
      geom_text_repel(aes(y = label_position, label = paste0(round(percent*100, 1), "%")), 
                vjust = 0.5, hjust = 0.5, size = 4, color = "white", box.padding = 0.01) -> plot_list[[2]]

df %>% 
  filter(celltype == 'NT(PAX6-/NKX6-1+)') %>%
  mutate(percent_ratio = percent/sum(percent)) %>%
  mutate(label_position = cumsum(percent_ratio) - percent_ratio/2) %>%
  mutate(group = factor(group, levels = c('v3-D7', 'v2-D7', 'v1-D7'))) %>%
  ggplot(aes(x="", y=percent_ratio, fill=group)) +
    geom_bar(stat="identity", width=1, color='white') +
    coord_polar("y", start=0) +
    scale_fill_manual(values = color_list) + 
    theme_void() + 
    labs(fill = 'Sample', title = 'NT(PAX6-/NKX6-1+)') + 
    theme(plot.title = element_text(hjust = 0.5)) +
      geom_text_repel(aes(y = label_position, label = paste0(round(percent*100, 1), "%")), 
                vjust = 0.5, hjust = 0.5, size = 4, color = "white", box.padding = 0.01) -> plot_list[[3]]

c = 'FP(FOXA1+)'
gene <- 'FOXA1'
expr_data <- FetchData(object = srt, vars = c(gene), layer = 'data', assay = 'decontX')
expr_data <- cbind(meta_df, expr_data[rownames(meta_df), gene])
colnames(expr_data) <- c(colnames(meta_df), gene)

expr_data['expr_type'] <- 'NA'
expr_data['expr_type'] <- ifelse(expr_data[gene] > 0, paste0(gene,'+'), paste0(gene,'-'))
expr_data['expr_type'] <- paste0(expr_data$celltype,"(",expr_data$expr_type,")")

table(expr_data$orig.ident, expr_data$expr_type) %>% data.frame() -> cellcount_df
colnames(cellcount_df) <- c('group', 'celltype', 'count')
table(expr_data$orig.ident) %>% data.frame() -> sample_cellcount_df
rownames(sample_cellcount_df) <- sample_cellcount_df$Var1
cellcount_df['sample_number'] <- sample_cellcount_df[cellcount_df$group, 'Freq']
cellcount_df['percent_per_sample'] = cellcount_df$count*100/cellcount_df$sample_number
cellcount_df <- cellcount_df %>% filter(celltype == c)
cellcount_df['percent'] <- cellcount_df$percent_per_sample/sum(cellcount_df$percent_per_sample)
write.csv(cellcount_df, 'revised_table_pie_4_foxa1_count.csv')
table(expr_data$orig.ident, expr_data$expr_type) -> celltype_prop
celltype_prop <- celltype_prop/rowSums(celltype_prop)
df <- data.frame(celltype_prop)
colnames(df) <- c('group', 'celltype', 'percent')
# write.csv(df, 'table_pie_4_foxa1_percent.csv')

df %>% 
  filter(celltype == c) %>%
  mutate(percent_ratio = percent/sum(percent)) %>%
  mutate(label_position = cumsum(percent_ratio) - percent_ratio/2) %>%
  mutate(group = factor(group, levels = c('v3-D7', 'v2-D7', 'v1-D7'))) %>%
    ggplot(aes(x="", y=percent_ratio, fill=group)) +
      geom_bar(stat="identity", width=1, color='white') +
      coord_polar("y", start=0) +
      scale_fill_manual(values = color_list) + 
      theme_void() +
      labs(fill = 'group', title = c) +
      theme(legend.position = 'na', plot.title = element_text(hjust = 0.5)) +
      geom_text_repel(aes(y = label_position, label = paste0(round(percent*100, 1), "%")), 
                vjust = 0.5, hjust = 0.5, size = 4, color = "white", box.padding = 0.01) -> plot_list[[4]]

c = 'SM(PAX3+)'
gene <- 'PAX3'
expr_data <- FetchData(object = srt, vars = c(gene), layer = 'data', assay = 'decontX')
expr_data <- cbind(meta_df, expr_data[rownames(meta_df), gene])
colnames(expr_data) <- c(colnames(meta_df), gene)

expr_data['expr_type'] <- 'NA'
expr_data['expr_type'] <- ifelse(expr_data[gene] > 0, paste0(gene,'+'), paste0(gene,'-'))
expr_data['expr_type'] <- paste0(expr_data$celltype,"(",expr_data$expr_type,")")

table(expr_data$orig.ident, expr_data$expr_type) %>% data.frame() -> cellcount_df
colnames(cellcount_df) <- c('group', 'celltype', 'count')
table(expr_data$orig.ident) %>% data.frame() -> sample_cellcount_df
rownames(sample_cellcount_df) <- sample_cellcount_df$Var1
cellcount_df['sample_number'] <- sample_cellcount_df[cellcount_df$group, 'Freq']
cellcount_df['percent_per_sample'] = cellcount_df$count*100/cellcount_df$sample_number
cellcount_df <- cellcount_df %>% filter(celltype == c)
cellcount_df['percent'] <- cellcount_df$percent_per_sample/sum(cellcount_df$percent_per_sample)
write.csv(cellcount_df, 'revised_table_pie_5_pax3_count.csv')
table(expr_data$orig.ident, expr_data$expr_type) -> celltype_prop
celltype_prop <- celltype_prop/rowSums(celltype_prop)
df <- data.frame(celltype_prop)
colnames(df) <- c('group', 'celltype', 'percent')
# write.csv(df, 'table_pie_5_pax3_percent.csv')

df %>% 
  filter(celltype == c) %>%
    mutate(percent_ratio = percent/sum(percent)) %>%
    mutate(label_position = cumsum(percent_ratio) - percent_ratio/2) %>%
    mutate(group = factor(group, levels = c('v3-D7', 'v2-D7', 'v1-D7'))) %>%
    ggplot(aes(x="", y=percent_ratio, fill=group)) +
      geom_bar(stat="identity", width=1, color='white') +
      coord_polar("y", start=0) +
      scale_fill_manual(values = color_list) + 
      theme_void() +
      labs(fill = 'group', title = c) +
      theme(legend.position = 'na', plot.title = element_text(hjust = 0.5)) +
      geom_text_repel(aes(y = label_position, label = paste0(round(percent*100, 1), "%")), 
                vjust = 0.5, hjust = 0.5, size = 4, color = "white", box.padding = 0.01) -> plot_list[[5]]

c = 'SCL(PAX1+ or PAX9+)'
expr_data <- FetchData(object = srt, vars = c('PAX1', 'PAX9'), layer = 'data', assay = 'decontX')
expr_data <- cbind(meta_df, expr_data[rownames(meta_df),])

expr_data['expr_type_1'] <- 'NA'
expr_data['expr_type_2'] <- 'NA'
expr_data['expr_type_1'] <- ifelse(expr_data$PAX1 > 0, 'PAX1+', 'PAX1-')
expr_data['expr_type_2'] <- ifelse(expr_data$PAX9 > 0, 'PAX9+', 'PAX9-')
expr_data['expr_type'] <- ifelse(expr_data$expr_type_1 == 'PAX1+' | expr_data$expr_type_2 == 'PAX9+', 'PAX1+ or PAX9+', '0')
expr_data['expr_type'] <- paste0(expr_data$celltype,"(",expr_data$expr_type,")")

table(expr_data$orig.ident, expr_data$expr_type) %>% data.frame() -> cellcount_df
colnames(cellcount_df) <- c('group', 'celltype', 'count')
table(expr_data$orig.ident) %>% data.frame() -> sample_cellcount_df
rownames(sample_cellcount_df) <- sample_cellcount_df$Var1
cellcount_df['sample_number'] <- sample_cellcount_df[cellcount_df$group, 'Freq']
cellcount_df['percent_per_sample'] = cellcount_df$count*100/cellcount_df$sample_number
cellcount_df <- cellcount_df %>% filter(celltype == c)
cellcount_df['percent'] <- cellcount_df$percent_per_sample/sum(cellcount_df$percent_per_sample)
write.csv(cellcount_df, 'revised_table_pie_6_pax1_pax9_count.csv')
table(expr_data$orig.ident, expr_data$expr_type) -> celltype_prop
celltype_prop <- celltype_prop/rowSums(celltype_prop)
df <- data.frame(celltype_prop)
colnames(df) <- c('group', 'celltype', 'percent')
# write.csv(df, 'table_pie_6_pax1_pax9_percent.csv')

df %>% 
  filter(celltype == c) %>%
    mutate(percent_ratio = percent/sum(percent)) %>%
    mutate(label_position = cumsum(percent_ratio) - percent_ratio/2) %>%
    mutate(group = factor(group, levels = c('v3-D7', 'v2-D7', 'v1-D7'))) %>%
    ggplot(aes(x="", y=percent_ratio, fill=group)) +
      geom_bar(stat="identity", width=1, color='white') +
      coord_polar("y", start=0) +
      scale_fill_manual(values = color_list) + 
      theme_void() +
      labs(fill = 'group', title = 'SCL\n(PAX1+ or PAX9+)') + 
      theme(legend.position = 'na', plot.title = element_text(hjust = 0.5)) +
      geom_text_repel(aes(y = label_position, label = paste0(round(percent*100, 1), "%")), 
                vjust = 0.5, hjust = 0.5, size = 4, color = "white", box.padding = 0.01) -> plot_list[[6]]

c = 'Dermomyotome/Myogenic(MYF5+ or PAX7+)'
expr_data <- FetchData(object = srt, vars = c('MYF5', 'PAX7'), layer = 'data', assay = 'decontX')
expr_data <- cbind(meta_df, expr_data[rownames(meta_df),])

expr_data['expr_type_1'] <- 'NA'
expr_data['expr_type_2'] <- 'NA'
expr_data['expr_type_1'] <- ifelse(expr_data$MYF5 > 0, 'MYF5+', 'MYF5-')
expr_data['expr_type_2'] <- ifelse(expr_data$PAX7 > 0, 'PAX7+', 'PAX7-')
expr_data['expr_type'] <- ifelse(expr_data$expr_type_1 == 'MYF5+' | expr_data$expr_type_2 == 'PAX7+', 'MYF5+ or PAX7+', '0')
expr_data['expr_type'] <- paste0(expr_data$celltype,"(",expr_data$expr_type,")")

table(expr_data$orig.ident, expr_data$expr_type) %>% data.frame() -> cellcount_df
colnames(cellcount_df) <- c('group', 'celltype', 'count')
table(expr_data$orig.ident) %>% data.frame() -> sample_cellcount_df
rownames(sample_cellcount_df) <- sample_cellcount_df$Var1
cellcount_df['sample_number'] <- sample_cellcount_df[cellcount_df$group, 'Freq']
cellcount_df['percent_per_sample'] = cellcount_df$count*100/cellcount_df$sample_number
cellcount_df <- cellcount_df %>% filter(celltype == c)
cellcount_df['percent'] <- cellcount_df$percent_per_sample/sum(cellcount_df$percent_per_sample)
write.csv(cellcount_df, 'revised_table_pie_7_myf5_pax7_count.csv')
table(expr_data$orig.ident, expr_data$expr_type) -> celltype_prop
celltype_prop <- celltype_prop/rowSums(celltype_prop)
df <- data.frame(celltype_prop)
colnames(df) <- c('group', 'celltype', 'percent')
# write.csv(df, 'table_pie_7_myf5_pax7_percent.csv')

df %>% 
  filter(celltype == c) %>%
    mutate(percent_ratio = percent/sum(percent)) %>%
    mutate(label_position = cumsum(percent_ratio) - percent_ratio/2) %>%
    mutate(group = factor(group, levels = c('v3-D7', 'v2-D7', 'v1-D7'))) %>%
    ggplot(aes(x="", y=percent_ratio, fill=group)) +
      geom_bar(stat="identity", width=1, color='white') +
      coord_polar("y", start=0) +
      scale_fill_manual(values = color_list) + 
      theme_void() +
      labs(fill = 'group', title = 'Dermomyotome/Myogenic\n(MYF5+ or PAX7+)') + 
      theme(legend.position = 'na', plot.title = element_text(hjust = 0.5)) +
      geom_text_repel(aes(y = label_position, label = paste0(round(percent*100, 1), "%")), 
                vjust = 0.5, hjust = 0.5, size = 4, color = "white", box.padding = 0.01) -> plot_list[[7]]

c = 'Endotome(EBF2+)'
gene <- 'EBF2'
expr_data <- FetchData(object = srt, vars = c(gene), layer = 'data', assay = 'decontX')
expr_data <- cbind(meta_df, expr_data[rownames(meta_df), gene])
colnames(expr_data) <- c(colnames(meta_df), gene)

expr_data['expr_type'] <- 'NA'
expr_data['expr_type'] <- ifelse(expr_data[gene] > 0, paste0(gene,'+'), paste0(gene,'-'))
expr_data['expr_type'] <- paste0(expr_data$celltype,"(",expr_data$expr_type,")")

table(expr_data$orig.ident, expr_data$expr_type) %>% data.frame() -> cellcount_df
colnames(cellcount_df) <- c('group', 'celltype', 'count')
table(expr_data$orig.ident) %>% data.frame() -> sample_cellcount_df
rownames(sample_cellcount_df) <- sample_cellcount_df$Var1
cellcount_df['sample_number'] <- sample_cellcount_df[cellcount_df$group, 'Freq']
cellcount_df['percent_per_sample'] = cellcount_df$count*100/cellcount_df$sample_number
cellcount_df <- cellcount_df %>% filter(celltype == c)
cellcount_df['percent'] <- cellcount_df$percent_per_sample/sum(cellcount_df$percent_per_sample)
write.csv(cellcount_df, 'revised_table_pie_8_ebf2_count.csv')
table(expr_data$orig.ident, expr_data$expr_type) -> celltype_prop
celltype_prop <- celltype_prop/rowSums(celltype_prop)
df <- data.frame(celltype_prop)
colnames(df) <- c('group', 'celltype', 'percent')
# write.csv(df, 'table_pie_8_ebf2_percent.csv')

df %>% 
  filter(celltype == c) %>%
    mutate(percent_ratio = percent/sum(percent)) %>%
    mutate(label_position = cumsum(percent_ratio) - percent_ratio/2) %>%
    mutate(group = factor(group, levels = c('v3-D7', 'v2-D7', 'v1-D7'))) %>%
    ggplot(aes(x="", y=percent_ratio, fill=group)) +
      geom_bar(stat="identity", width=1, color='white') +
      coord_polar("y", start=0) +
      scale_fill_manual(values = color_list) + 
      theme_void() +
      labs(fill = 'group', title = c) + 
      theme(legend.position = 'na', plot.title = element_text(hjust = 0.5)) +
      geom_text_repel(aes(y = label_position, label = paste0(round(percent*100, 1), "%")), 
                vjust = 0.5, hjust = 0.5, size = 4, color = "white", box.padding = 0.01) -> plot_list[[8]]

options(repr.plot.width = 10, repr.plot.height = 5)
tiff('revised_fig_6_3_piechart.tiff', width = 10, height = 5, res = 300, units = 'in')
print(wrap_plots(plot_list, ncol = 4) +
  plot_layout(guides = "collect"))
dev.off()

###############################
#### Barplot gene positive ####
###############################
genes <- c('SHH', 'FOXA2', 'CHRD', 'NOG', 'PAX6', 'IRX3', 'DBX2', 'NKX6-1', 'OLIG2', 
    'FOXA1', 'FRZB', 'NKX6-1', 'PAX3', 'MEOX1', 'MYF5', 'PAX7', 
    'SCX', 'PAX1', 'PAX9', 'EBF2', 'ETV2')
celltypes <- c("Noto", 'Noto', "Noto", 'Noto', 'NT', 'NT', 'NT', 'NT', 'NT', 'FP', 'FP', 'FP', 
    'SM', 'SM', 'Dermomyotome/Myogenic', 'Dermomyotome/Myogenic', 'Syndetome',
    'SCL', 'SCL', 'Endotome', 'Endotome')
plots <- c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 6, 6, 7, 8, 8, 9, 9)

plots_df <- data.frame(gene = genes, celltype = celltypes, plot = plots)

p <- 1

current_plot_df <- plots_df %>% filter(plot == p)
c <- unique(current_plot_df$celltype)

# for (p in unique(plots_df$plot)){
for (p in c(1)){
    current_plot_df <- plots_df %>% filter(plot == p)
    c <- unique(current_plot_df$celltype)

    expr_data <- FetchData(object = srt, vars = current_plot_df$gene, layer = 'data', assay = 'decontX')
    all_df <- data.frame()
    for (g in unique(current_plot_df$gene)){
        meta_df <- srt@meta.data
        meta_df[g] <- expr_data[rownames(meta_df), g]
        meta_df[paste0('expr_type_', g)] <- ifelse(meta_df[g] > 0, paste0(g, '+'), paste0(g, '-'))
        meta_df[paste0('expr_type_', g)] <- paste0(meta_df$celltype," (",meta_df[,paste0('expr_type_', g)],")")
        table(meta_df$orig.ident, meta_df[,paste0('expr_type_', g)]) -> celltype_count
        # celltype_prop <- celltype_prop/rowSums(celltype_prop)
        df <- data.frame(celltype_count)
        colnames(df) <- c('group', 'celltype', 'count')
        # df <- df %>% filter(celltype == paste0(c," (",g,"+)"))
        df['gene'] <- g
        # write.csv(df, paste0('revised_supp_fig_6_5_cellbar_',g,'.csv'))
        all_df <- rbind(all_df, df)
    }

    # df['celltype'] <- NA
    # df['count'] <- NA
    # df['gene'] <- NA
    # all_df <- rbind(all_df, df)

    all_df$gene <- factor(all_df$gene, levels = c(unique(current_plot_df$gene), NA))
    all_df <- all_df %>% arrange(group)
    all_df['condition'] <- paste0(all_df$group, '_', all_df$celltype)
    all_df$condition <- factor(all_df$condition, levels = unique(all_df$condition))

    all_df['label'] <- as.character(all_df$gene)
    all_df['label'] <- paste0(all_df$label,"+")
    all_df[all_df$label == 'NA+', 'label'] <- ''
    all_df <- all_df %>% filter(grepl(c, celltype))
    all_df <- all_df %>% filter(!grepl('-)', condition))
    all_df %>% select(group) %>% unique() -> space_df
    space_df['celltype'] <- 'space'
    space_df['count'] <- 0
    space_df['gene'] <- 'space'
    space_df['condition'] <- 'space'
    space_df$condition <- paste0(space_df$group, "_", space_df$condition)
    space_df['label'] <- ''

    all_df <- rbind(all_df, space_df)
    all_df <- all_df %>% arrange(group, gene)
    rownames(all_df) <- 1:nrow(all_df)
    all_df$condition <- factor(all_df$condition, levels = unique(all_df$condition))

    max_y <- max(na.omit(all_df$count))
    segment_y <- max_y*1.1
    segment_x_1 <- length(current_plot_df$gene)
    segment_x_2 <- segment_x_1 + 1
    named_list <- setNames(all_df$label, as.character(all_df$condition))
    custom_colors <- setNames(color_list[1:segment_x_1], unique(current_plot_df$gene))

    options(repr.plot.width = 5, repr.plot.height = 3)
    all_df %>% ggplot(aes(x = condition, y = count)) + 
        geom_col(aes(fill = gene), color = 'black') +
        theme(panel.border = element_blank(),
            panel.background = element_blank(),  # Remove background
            panel.grid.major = element_blank(),  # Remove major gridlines
            panel.grid.minor = element_blank(),  # Remove minor gridlines
            axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
            axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = 'na',
            plot.title = element_text(hjust = 0.5)) + 
        scale_x_discrete(labels = named_list) +
        scale_fill_manual(values = custom_colors) +
        labs(y = 'Count', x = '', fill = 'Gene', title = c) + 
        scale_y_continuous(expand = c(0,0), limits = c(0, max_y*1.14)) +
        geom_segment(aes(x = 0.6, xend = 0.6+segment_x_1-0.2, y = max_y*1.03, yend = max_y*1.03), size = 0.5) +
        annotate("text", x = (0.6+(0.6+segment_x_1-0.2))/2, y = max_y*1.09, label = "v1-D7", size = 4) +
        geom_segment(aes(x = 0.6+segment_x_2, xend = 0.6+segment_x_2+segment_x_1-0.2, y = max_y*1.03, yend = max_y*1.03), size = 0.5) +
        annotate("text", x = ((0.6+segment_x_2)+(0.6+segment_x_2+segment_x_1-0.2))/2, y = max_y*1.09, label = "v2-D7", size = 4) +
        geom_segment(aes(x = 0.6+segment_x_2*2, xend = 0.6+segment_x_2*2+segment_x_1-0.2, y = max_y*1.03, yend = max_y*1.03), size = 0.5) +
        annotate("text", x = ((0.6+segment_x_2*2)+(0.6+segment_x_2*2+segment_x_1-0.2))/2, y = max_y*1.09, label = "v3-D7", size = 4) -> p1
    
    tiff(paste0('revised_supp_fig_6_5_cellbar_',p,'.tiff'), width = 4, height = 3, res = 300, units = 'in')
    print(p1)
    dev.off()
}