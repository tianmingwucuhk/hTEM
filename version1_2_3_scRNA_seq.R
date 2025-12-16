*********************************************************************
********************V1_QC&integration&annotation*********************
*********************************************************************


library(Seurat)
library(tidyverse)
library(scutilsR)

samples <- list.dirs("data/10x_matrix/", full.names = F, recursive = F)
seu.list <- pbapply::pblapply(samples, function(sn) {
  counts <- Read10X(file.path("data/10x_matrix/", sn))
  sn <- gsub("_", "-", sn) # 注意"_"在`CreateSeuratObject()`里有特殊的意义
  colnames(counts) <- paste(sn, colnames(counts), sep = "_")
  seu <- CreateSeuratObject(counts = counts,min.cells = 3, min.features = 200)
  return(seu)
})

NTSM <- base::Reduce(f = merge, x = seu.list)


rm(seu.list)
gc()

#### QC ####
table(NTSM$orig.ident)

NTSM <- PercentageFeatureSet(NTSM, pattern = "^MT-", col.name = "percent.mt")
NTSM <- PercentageFeatureSet(NTSM, pattern = "^RPS", col.name = "percent.rps")
NTSM <- PercentageFeatureSet(NTSM, pattern = "^RPL", col.name = "percent.rpl")

VlnPlot(NTSM, features = c("nFeature_RNA"), group.by = "orig.ident", pt.size = 0) + NoLegend() +
  geom_hline(yintercept = c(1000, 5000), linetype = "dashed", color = "blue")

VlnPlot(NTSM, features = c("nCount_RNA"), group.by = "orig.ident", pt.size = 0, log = T) + NoLegend() + 
  geom_hline(yintercept = c(3000, 50000), linetype = "dashed", color = "blue") + NoLegend()

VlnPlot(NTSM, features = c("percent.mt"), group.by = "orig.ident", pt.size = 0, log = T) +
  geom_hline(yintercept = c(10), linetype = "dashed", color = "blue")

NTSM <- subset(NTSM, subset = nFeature_RNA > 200 & percent.mt < 5)

NTSM <- NormalizeData(NTSM)
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
NTSM <- CellCycleScoring(NTSM, s.features = s.genes, g2m.features = g2m.genes)

NTSM <- scutilsR::MarkDoublets(NTSM, split.by = "orig.ident")

#NTSM <- subset(NTSM,subset = DF.classifications=='Singlet')

QuickCluster <- function(object) {
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object, nfeatures = 2000)
  object <- ScaleData(object)
  object <- RunPCA(object)
  object <- FindNeighbors(object, reduction = "pca", dims = 1:30)
  object <- FindClusters(object)
  return(object)
}
seu.list <- SplitObject(NTSM, split.by = "orig.ident")
seu.list <- lapply(seu.list, QuickCluster)
clusters <- lapply(seu.list, function(xx) xx$seurat_clusters) %>% base::Reduce(c, .)
NTSM$quick_clusters <- clusters[rownames(NTSM@meta.data)]
rm(seu.list)
gc()

NTSM <- scutilsR::RemoveAmbientRNAs(NTSM, split.by = "orig.ident", cluster.name = "quick_clusters")
NTSM 
head(NTSM@meta.data)

qs::qsave(NTSM, file = 'output/01_output/NTSM_downstream_tmp_adjustment.qs')

##################batch correction and integration
NTSM <- qs::qread(file = 'output/01_output/NTSM_downstream_tmp_adjustment.qs')

###remove doublets
NTSM <- subset(NTSM,subset = DF.classifications=='Singlet')
NTSM <- NormalizeData(NTSM,normalization.method = 'LogNormalize',scale.factor = 10000,verbose = TRUE)
NTSM <- FindVariableFeatures(NTSM, selection.method = "vst", nfeatures = 2000, verbose=FALSE)
NTSM <- CellCycleScoring(NTSM, s.features = s.genes, g2m.features = g2m.genes)
NTSM <- ScaleData(NTSM, vars.to.regress = c("S.Score", "G2M.Score"))
NTSM <- RunPCA(NTSM)

NTSM <- RunUMAP(NTSM,dims = 1:30,reduction.name='umap.naive')

#integrate with harmony
library(harmony)
NTSM <- RunHarmony(NTSM, group.by.vars=c('orig.ident'),reduction.use='pca',dim.use=1:50)
NTSM <- RunUMAP(NTSM,reduction = 'harmony',dims = 1:30,n.neighbors = 200,reduction.name = 'umap.harmony')

#integrate with bbknn
library(bbknnR)
NTSM <- RunBBKNN(NTSM, reduction ='pca',run_UMAP = TRUE, run_TSNE = FALSE, UMAP_name ='umap.bbknn',
                 batch_key='orig.ident')
p1 =DimPlot(NTSM,reduction = 'umap.harmony',group.by = 'orig.ident')
p2 =DimPlot(NTSM,reduction = 'umap.bbknn', group.by = 'orig.ident')
p1+p2


GOI <- c('NOTO','FOXA1','FOXA2','TBXT','SOX9', 'FOXJ1','CDX2','HES7',
         'CHRD','NOG','SHH','CER1','NODAL','EGF','LEFTY1','LEFTY2',
         'BMP7','WNT3A','FGF8','ARL13B', 'KRT18', 'KRT19', 'CDH1', 'TIMP3', 'SEMA3C','SEMA3E',
         'SOX2','TBXT','PAX3','CYP26A1','EPHA1','FGF3','FGF4','FGF8','FGF17','FGF19',
         'WNT3A','WNT5A','WNT8A','NEFM',
         'TBXT', 'TBX6','HES7','MSGN1','DLL1','DLL3','LEF1','DKK1','FGF18',
         'MESP1','MESP2','RIPPLY2','CER1','DKK1','DLL1','DLL3', 'LEF1','TCF15','MEOX1','BMP4',
         'UNCX','LFNG','TBX18','FST','SIX1','EYA1','MEOX1','DCHS1','ALDH1A2',
         'BMP5','PAX3','MYF5','CRABP2','RXRG','RARA',
         'PAX3','PAX7','EBF2','TGFB1','TGFBI','FOXC2','SCX','ALDH1A2','TIMP3','RARA',
         'FRZB','PAX1','PAX9','NKX3-2','SOX9','SEMA3C','TGFBI',
         "TWIST1",'SNAI1','FOXC2','WNT5B','GLI1',
         'TBXT','SOX2','HES7','NKX1-2','NEFM','FGF8','FGF17','FGF19',
         'BMP7','WNT3A','WNT5A','WNT8A','SPRY2','SOX2','WNT5B','RARB','CRABP2','HES4','HES5','CCND1','CELSR1',
         'PAX6','IRX3','IRX5','DBX1','DBX2','WNT1','PAX7','NKX6-1','NKX6-2','OLIG2','NKX2-2',
         'FOXA1','FOXA2','SHH','TIMP3','FOXJ1','SEMA3E','DHRS3','RARB',
         'CRABP1','HES4','HES5','DLL1','EBF2','DLL3','ISL1',
         'FGF8','NEFM','CYP26A1','NKX1-2','HES7','WNT5A','WNT8A',
         'CXCR4','CDX2','MSX1','RARG','DHRS3',
         'PAX2','PAX8','OSR1','LHX1','NEFM','DKK1','KRT18','FGF8','RDH10',
         'CD34','KDR','ETV2','SOX17','CLDN5','TIMP3','FOXC2','TGFB1')

tiff(filename = "output/01_output/tmp/v1_integrate_adjustment_RNA_bbknn_nodoublets.tiff", 
     width = 3000, height = 3300, units = "px")
FeaturePlot(NTSM, features = unique(GOI), pt.size = 0.2,reduction = 'umap.bbknn',
            ncol = 10, order = TRUE)
dev.off()

qs::qsave(NTSM, file = 'output/01_output/02_NTSM_downstream_RNA_nodoublets.qs')
NTSM <- qs::qread(file = 'output/01_output/02_NTSM_downstream_RNA_nodoublets.qs')


#####annotation
######

NTSM <- qs::qread(file = 'output/01_output/02_NTSM_downstream_RNA_nodoublets.qs')

resolutions <- c(0.6,1,1.2,1.5,2.0,2.2,2.4)

NTSM <- FindNeighbors(NTSM, reduction = "harmony", dims = 1:20, k.param = 20)
NTSM <- FindClusters(NTSM, resolution = resolutions)

DimPlot(NTSM, reduction = "umap.bbknn", group.by = paste0("RNA_snn_res.", resolutions), ncol = 3, label = T) & NoLegend()

ggsave(filename = 'v1_integration_clusters2_bbknn.tiff',width = 14,height =11,dpi = 300,path = 'output/02_output/Figure/')

FeaturePlot(NTSM, reduction = 'umap.harmony',split.by = 'RNA_snn_res.3',features = c('NOTO'))
Idents(NTSM) <- NTSM$RNA_snn_res.0.6

NTSM <- RenameIdents(NTSM,
                     '0' = 'Neural Tube',
                     '10'  = 'Neural Tube',
                     '2' = 'Neural Tube',
                     '4'  = 'Somite',
                     '6' = 'Somite',
                     '5' = 'Somite',
                     '7' = 'Somite',
                     '8' = 'Somite',
                     '1' = 'Somite',
                     '3' = 'Somite',
                     '9' = 'Endothelium')

DimPlot(NTSM, reduction = "umap.harmony",label = T)
FeaturePlot(NTSM, reduction = 'umap.bbknn',features = c('decontX_contamination'))
NTSM$celltype=Idents(NTSM)
NTSM$celltype=as.character(NTSM$celltype)

qs::qsave(NTSM, file = 'output/02_output/01_NTSM_annotation.qs')
saveRDS(NTSM,file = 'output/02_output/01_NTSM_v1_annotation.rds')

NTSM <- qs::qread(file = 'output/02_output/01_NTSM_annotation.qs')



###features##
Markers<- c("TBXT", "HES7", "TBX6", "DLL1", "MSGN1", "PAX3", "MEOX1", "SOX2", "NKX1-2", "HES4", "PAX6", "ONECUT1", "NOTO", "FOXA2", "SHH",
            "ETV2", "MEF2C", "NOTO", "CHRD", "FOXA2", "FOXA1", "WNT8A", "FGF3", "FGF4", "FGF8", "FGF17", "WNT3A", "LEF1", "WNT5A", "WNT5B", 
            "CYP26A1", "MESP1", "MESP2", "CER1", "RIPPLY2", "TBX18", "UNCX", "HES1", "FOXC1", "FOXC2", "ALDH1A2", "IRX3", "IRX5", "GRHL2", 
            "SOX1", "SOX21", "RARB", "RXRG", "SOX17", "CD34", "CLDN5", "KDR")

for (marker in Markers) {
  filename <- paste0("output/03_output/featureplot_tiff/", marker, ".tiff")
  tiff(filename = filename, width = 1320, height = 1200,units = 'px',res = 300)
  print(FeaturePlot(NTSM,features = marker,reduction = 'umap.harmony',pt.size = 0.2,order = T))
  dev.off()
  message("Saved: ", filename)  # 打印保存成功信息
}

tiff(filename = "output/03_output/featureplot_tiff/all_markers.tiff", 
     width = 9000, height = 9900, units = "px",res = 300)
FeaturePlot(NTSM, features = unique(Markers), pt.size = 0.2,reduction = "umap.harmony",
            ncol = 7, order = TRUE)
dev.off()


********************V2_QC&integration&annotation*********************

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

NTSM <- base::Reduce(f = merge, x = seu.list)

rm(seu.list)
gc()

table(NTSM$orig.ident)
#########qc################
NTSM <- PercentageFeatureSet(NTSM, pattern = "^MT-", col.name = "percent.mt")
NTSM <- PercentageFeatureSet(NTSM, pattern = "^RPS", col.name = "percent.rps")
NTSM <- PercentageFeatureSet(NTSM, pattern = "^RPL", col.name = "percent.rpl")

VlnPlot(NTSM, features = c("nFeature_RNA"), group.by = "orig.ident", pt.size = 0) + NoLegend() +
  geom_hline(yintercept = c(1000, 5000), linetype = "dashed", color = "blue")

VlnPlot(NTSM, features = c("nCount_RNA"), group.by = "orig.ident", pt.size = 0, log = T) + NoLegend() + 
  geom_hline(yintercept = c(3000, 50000), linetype = "dashed", color = "blue") + NoLegend()


VlnPlot(NTSM, features = c("percent.mt"), group.by = "orig.ident", pt.size = 0, log = T) +
  geom_hline(yintercept = c(10), linetype = "dashed", color = "blue")

NTSM <- subset(NTSM, subset = nFeature_RNA > 200 & percent.mt < 5)

NTSM <- NormalizeData(NTSM)
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
NTSM <- CellCycleScoring(NTSM, s.features = s.genes, g2m.features = g2m.genes)

NTSM <- scutilsR::MarkDoublets(NTSM, split.by = "orig.ident")

QuickCluster <- function(object) {
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object, nfeatures = 2000)
  object <- ScaleData(object)
  object <- RunPCA(object)
  object <- FindNeighbors(object, reduction = "pca", dims = 1:30)
  object <- FindClusters(object)
  return(object)
}
seu.list <- SplitObject(NTSM, split.by = "orig.ident")
seu.list <- lapply(seu.list, QuickCluster)
clusters <- lapply(seu.list, function(xx) xx$seurat_clusters) %>% base::Reduce(c, .)
NTSM$quick_clusters <- clusters[rownames(NTSM@meta.data)]
rm(seu.list)
gc()

NTSM <- scutilsR::RemoveAmbientRNAs(NTSM, split.by = "orig.ident", cluster.name = "quick_clusters")
NTSM 
head(NTSM@meta.data)

qs::qsave(NTSM, file = 'output/NTSM_downstream_tmp_adjustment.qs')

######batch correction and integration
library(Seurat)
library(tidyverse)
library(scutilsR)
library(harmony)

NTSM <- qs::qread(file = 'output/01_output/01_NTSM_downstream_tmp_adjustment.qs')

# gene complexity score
NTSM <- NTSM %>%
  AddMetaData(
    metadata = data.frame(
      GeneComplexityScore = NTSM$nFeature_RNA / NTSM$nCount_RNA))

NTSM <- NormalizeData(NTSM,normalization.method = 'LogNormalize',scale.factor = 10000,verbose = TRUE)
NTSM <- FindVariableFeatures(NTSM, selection.method = "vst", nfeatures = 1500, verbose=FALSE)

NTSM <- ScaleData(NTSM, vars.to.regress = c("S.Score", "G2M.Score","GeneComplexityScore","percent.rps","percent.rpl"))
NTSM <- RunPCA(NTSM)
NTSM <- RunHarmony(NTSM, group.by.vars=c('batch'),reduction.use='pca',dim.use=1:50)
NTSM <- RunUMAP(NTSM,reduction = 'harmony',dims = 1:30,n.neighbors = 200,reduction.name = 'umap.harmony')
library(bbknnR)
NTSM <- RunBBKNN(NTSM, reduction ='pca',run_UMAP = TRUE, run_TSNE = FALSE, UMAP_name ='umap.bbknn',
                 batch_key='batch')

p1 =DimPlot(NTSM,reduction = 'umap.harmony',group.by = 'orig.ident')
p2 =DimPlot(NTSM,reduction = 'umap.bbknn', group.by = 'orig.ident')
p1+p2
DimPlot(NTSM, reduction = 'umap.harmony',split.by = 'orig.ident')

tiff(filename = "output/GOI/v2_integrate_adjustment_harmony_nodoublets_GCP_rp_regress2.tiff", 
     width = 3000, height = 3300, units = "px")
FeaturePlot(NTSM, features = unique(GOI), pt.size = 0.2,reduction = 'umap.harmony',
            ncol = 10, order = TRUE)
dev.off()

qs::qsave(NTSM,file = 'output/01_output/v2_integration.qs')

###check markers
NTSM2 <- qs::qread(file = 'output/01_output/02_v2_integration.qs')
GOI <- c('NOTO','FOXA1','FOXA2','TBXT','SOX9', 'FOXJ1','CDX2','HES7',
         'CHRD','NOG','SHH','CER1','NODAL','EGF','LEFTY1','LEFTY2',
         'BMP7','WNT3A','FGF8','ARL13B', 'KRT18', 'KRT19', 'CDH1', 'TIMP3', 'SEMA3C','SEMA3E',
         'SOX2','TBXT','PAX3','CYP26A1','EPHA1','FGF3','FGF4','FGF8','FGF17','FGF19',
         'WNT3A','WNT5A','WNT8A','NEFM',
         'TBXT', 'TBX6','HES7','MSGN1','DLL1','DLL3','LEF1','DKK1','FGF18',
         'MESP1','MESP2','RIPPLY2','CER1','DKK1','DLL1','DLL3', 'LEF1','TCF15','MEOX1','BMP4',
         'UNCX','LFNG','TBX18','FST','SIX1','EYA1','MEOX1','DCHS1','ALDH1A2',
         'BMP5','PAX3','MYF5','CRABP2','RXRG','RARA',
         'PAX3','PAX7','EBF2','TGFB1','TGFBI','FOXC2','SCX','ALDH1A2','TIMP3','RARA',
         'FRZB','PAX1','PAX9','NKX3-2','SOX9','SEMA3C','TGFBI',
         "TWIST1",'SNAI1','FOXC2','WNT5B','GLI1',
         'TBXT','SOX2','HES7','NKX1-2','NEFM','FGF8','FGF17','FGF19',
         'BMP7','WNT3A','WNT5A','WNT8A','SPRY2','SOX2','WNT5B','RARB','CRABP2','HES4','HES5','CCND1','CELSR1',
         'PAX6','IRX3','IRX5','DBX1','DBX2','WNT1','PAX7','NKX6-1','NKX6-2','OLIG2','NKX2-2',
         'FOXA1','FOXA2','SHH','TIMP3','FOXJ1','SEMA3E','DHRS3','RARB',
         'CRABP1','HES4','HES5','DLL1','EBF2','DLL3','ISL1',
         'FGF8','NEFM','CYP26A1','NKX1-2','HES7','WNT5A','WNT8A',
         'CXCR4','CDX2','MSX1','RARG','DHRS3',
         'PAX2','PAX8','OSR1','LHX1','NEFM','DKK1','KRT18','FGF8','RDH10',
         'CD34','KDR','ETV2','SOX17','CLDN5','TIMP3','FOXC2','TGFB1')

for (marker in GOI) {
  filename <- paste0("output/01_output/Featureplot/", marker, ".tiff")
  tiff(filename = filename, width = 5.4, height = 4.5, units = "in", res = 300)
  print(FeaturePlot(hTLOV2,features = marker,reduction = 'umap.harmony',pt.size = 0.2,order = T))
  dev.off()
  message("Saved: ", filename)  # 打印保存成功信息
}


#cluster
resolutions <- c(2.0,2.2,2.4,2.6,2.8,3.0,3.4)

NTSM <- FindNeighbors(NTSM, reduction = "harmony", dims = 1:20, k.param = 20)
NTSM <- FindClusters(NTSM, resolution = resolutions)


DimPlot(NTSM, reduction = "umap.harmony", group.by = paste0("RNA_snn_res.", resolutions), ncol = 3, label = T) & NoLegend()
ggsave(filename = 'v2_integration_clusters_harmony.tiff',width = 14,height =11,dpi = 300,path = 'output/01_output/Figures/')

qs::qsave(NTSM,file = 'output/01_output/v2_integration.qs')
saveRDS(NTSM,file='output/01_output/v2_integration.rds')

#
DimPlot(NTSM,reduction = 'umap.harmony',split.by = 'orig.ident',group.by = 'orig.ident')

###############annotation##########
hTLOV2=qs::qread(file = 'output/01_output/02_v2_integration.qs')
metadata <- read.csv(file = 'data/20250730_htlov2_metadata.csv',row.names = 1)

m <- metadata[,c(25,26)]

v2 <- AddMetaData(hTLOV2,metadata = m)

DimPlot(v2, group.by = 'bbknn_res.1.6',label = T)
Idents(v2) = v2$bbknn_res.1.6
v2 <- RenameIdents(v2,
                   '24' = 'Node-like',
                   '17' = 'NMP-Meso/pPSM',
                   '8' = 'NMP-Neural',
                   '2' = 'Caud.Meso',
                   '18' = 'Caud.NP',
                   '3' = 'Dorsal NT',
                   '11' = 'Ventral NT',
                   '10' = 'aPSM',
                   '0' = 'E-SM',
                   '1' = 'E-SM',
                   '5' = 'E-SM',
                   '6' = 'E-SM',
                   '13' = 'E-SM',
                   '7' = 'E-SM',
                   '20'='E-SM',
                   '14'='M-SM',
                   '16'='M-SM',
                   '23' = 'M-SM',
                   '4'='L-SM',
                   '9' = 'L-SM',
                   '19' = 'L-SM',
                   '12' = 'SAG induced',
                   '15' = 'SAG induced',
                   '22' = 'Neuron',
                   '21' = 'EC')

v2_Cell <- setNames(
  c("#FF1677", "#00FB0D", "#FAA31C", "#B6CAFA", "#CE471C", "#D8EA60", 
    "#D58A81", "#F96ADB", "#A0EBF1", "#A0EBF1", "#A0EBF1", "#A0EBF1", 
    "#A0EBF1", "#A0EBF1", "#A0EBF1", "#0D8AFF", "#0D8AFF", "#0D8AFF", 
    "#16478B", "#16478B", "#16478B", "#2E7300", "#2E7300", "#FED426", 
    "#9A8116"),
  c("Node-like", "NMP-Meso/pPSM", "NMP-Neural", "Caud.Meso", "Caud.NP", "Dorsal NT", 
    "Ventral NT", "aPSM", "E-SM", "E-SM", "E-SM", "E-SM", 
    "E-SM", "E-SM", "E-SM", "M-SM", "M-SM", "M-SM", 
    "L-SM", "L-SM", "L-SM", "SAG induced", "SAG induced", "Neuron", "EC")
)

DimPlot(v2, reduction = "umap.harmony") +scale_color_manual(values = v2_Cell)
ggsave(filename = 'output/01_output/03_v2_annotation/v2_dimplot.tiff',height =5.3,width = 8,dpi = 300 )

DimPlot(v2, reduction = "umap.harmony",split.by = 'orig.ident') +scale_color_manual(values = v2_Cell)
ggsave(filename = 'output/01_output/03_v2_annotation/v2_dimplot_split.tiff',height =5.3,width = 17,dpi = 300 )

qs::qsave(v2,file = 'output/01_output/03_v2_annotation/v2_annotation.qs')

v2 <- qs::qread(file = 'output/01_output/03_v2_annotation/v2_annotation.qs')
saveRDS(v2, file = 'output/01_output/03_v2_annotation/v2_annotation.rds')
metadata <- v2@meta.data
write.csv(metadata,file = 'output/01_output/03_v2_annotation/v2_metadata.csv')

#########dotplot###########

s <- subset(v2, orig.ident %in% c('v2-D5','v2-D7'))
ss <- subset(s,celltype %in% c('Node-like','E-SM','M-SM','L-SM','Dorsal NT','Ventral NT/FP'))
ss$celltype <- factor(ss$celltype,levels = c('Node-like','E-SM','M-SM','L-SM','Dorsal NT','Ventral NT/FP'))
GOI <- c(
  "NOTO", "CHRD", "TBXT", "SOX2",
  "TCF15", "MEOX1", "COL5A2", "NR2F2", "UNCX", "PAX3", "SIX1", "SNAI2", 
  "EBF2", "TFBI", "TWIST1", "COL1A1", "SCX", 
  "PAX1", "PAX9", "NKX3-2",
  "PAX7", "MYF5", "HES1", "MSX1", "DBX2", "PAX6", "IRX3", "GLI1", "HES4", 
  "NKX6-1", "NKX2-8", "OLIG1", "OLIG2", "FOXA2", "SHH"
)


dotplot <- DotPlot(ss, features =GOI, group.by = 'celltype',cols = c('lightblue','brown'),dot.scale = 10) + RotatedAxis() + ggtitle('v2')+
  theme(plot.title = element_text(hjust = 0.5))

dotplot_data <- dotplot$data
ggplot(dotplot_data, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_size_continuous(name = "Percent Expressed", range = c(1, 10),
                        limits = c(5,100),
                        breaks = c(5,25,50,75)) +
  scale_color_gradient(low = "lightblue", high = "brown") +
  labs(x = "Features", y = "Identity",  title = "v2") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(color = "black",size = 10),                    # 坐标轴刻度文字颜色
    axis.title = element_text(color = "black"),                   # 坐标轴标题文字颜色
    legend.text = element_text(color = "black"),                  # 图例文字颜色
    legend.title = element_text(color = "black"), 
    panel.border = element_rect(color = "grey", fill = NA, linewidth = 1),
    panel.grid.major = element_line(color = "grey", linewidth = 0.1, linetype = "solid"),
    panel.grid.minor = element_line(color = "grey", linewidth = 0.1, linetype = "solid"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave(filename="output/01_output/Figures/v2_subseted_dotplot.pdf",width=15,height=3.9,dpi=300)

*********************************************************************
********************V3_QC&integration&annotation&analysis************
*********************************************************************


******************************************************************
*************************01.V3_QC_integration_annotation**********
******************************************************************
library(Seurat)
library(tidyverse)
library(scutilsR)

#####QC
samples <- list.dirs("data/10X_matrix/", full.names = F, recursive = F)
seu.list <- pbapply::pblapply(samples, function(sn) {
  counts <- Read10X(file.path("data/10X_matrix/", sn))
  sn <- gsub("_", "-", sn)
  colnames(counts) <- paste(sn, colnames(counts), sep = "_")
  seu <- CreateSeuratObject(counts = counts,min.cells = 3, min.features = 200)
  return(seu)
})
v3 <- base::Reduce(f = merge, x = seu.list)
rm(seu.list)
gc()
v3 <- PercentageFeatureSet(v3, pattern = "^MT-", col.name = "percent.mt")
v3 <- PercentageFeatureSet(v3, pattern = "^RPS", col.name = "percent.rps")
v3 <- PercentageFeatureSet(v3, pattern = "^RPL", col.name = "percent.rpl")

v3 <- subset(v3, subset = nFeature_RNA > 200 & percent.mt < 5)
##cellcycle
v3 <- NormalizeData(v3)
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
v3 <- CellCycleScoring(v3, s.features = s.genes, g2m.features = g2m.genes)
##doublets
v3 <- scutilsR::MarkDoublets(v3, split.by = "orig.ident")
##RNA contamination
QuickCluster <- function(object) {
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object, nfeatures = 2000)
  object <- ScaleData(object)
  object <- RunPCA(object)
  object <- FindNeighbors(object, reduction = "pca", dims = 1:30)
  object <- FindClusters(object)
  return(object)
}
seu.list <- SplitObject(v3, split.by = "orig.ident")
seu.list <- lapply(seu.list, QuickCluster)
clusters <- lapply(seu.list, function(xx) xx$seurat_clusters) %>% base::Reduce(c, .)
NTSM$quick_clusters <- clusters[rownames(NTSM@meta.data)]
rm(seu.list)
gc()
v3 <- scutilsR::RemoveAmbientRNAs(v3, split.by = "orig.ident", cluster.name = "quick_clusters")
qs::qsave(v3, file = 'output/NTSM_downstream_tmp_adjustment.qs')


##batch correction and integration
v3 <- subset(v3,subset = DF.classifications=='Singlet')
DefaultAssay(v3)<- 'decontX'

v3 <- PercentageFeatureSet(v3, pattern = "^MT-", col.name = "percent.mt")
v3 <- PercentageFeatureSet(v3, pattern = "^RPS", col.name = "percent.rps")
v3 <- PercentageFeatureSet(v3, pattern = "^RPL", col.name = "percent.rpl")
v3 <- subset(v3, subset = nFeature_decontX > 200 & percent.mt < 5)

v3 <- NormalizeData(v3,normalization.method = 'LogNormalize',scale.factor = 10000,verbose = TRUE)
v3 <- FindVariableFeatures(v3, selection.method = "vst", nfeatures = 2000, verbose=FALSE)
v3 <- CellCycleScoring(v3, s.features = s.genes, g2m.features = g2m.genes)

v3 <- ScaleData(v3, vars.to.regress = c("S.Score", "G2M.Score"))
v3 <- RunPCA(v3)

v3 <- RunUMAP(v3,dims = 1:30,reduction.name='umap.naive')
DimPlot(v3,reduction = "umap.naive")
#integrate with harmony
library(harmony)
v3 <- RunHarmony(v3, group.by.vars=c('orig.ident'),reduction.use='pca',dim.use=1:50)
v3 <- RunUMAP(v3,reduction = 'harmony',dims = 1:30,n.neighbors = 200,reduction.name = 'umap_no_doublets_decontX.harmony')

#integrate with bbknn
library(bbknnR)
v3 <- RunBBKNN(v3, reduction ='pca',run_UMAP = TRUE, run_TSNE = TRUE, UMAP_name ='umap.bbknn',TSNE_name = 'tsne.bbknn',
                 batch_key='orig.ident')

p1 =DimPlot(v3,reduction = 'tsne.bbknn',group.by = 'orig.ident')
p2 =DimPlot(v3,reduction = 'umap.bbknn', group.by = 'orig.ident')
p1+p2


#check markers
GOI <- c('NOTO','FOXA1','FOXA2','TBXT','SOX9', 'FOXJ1','CDX2','HES7',
         'CHRD','NOG','SHH','CER1','NODAL','EGF','LEFTY1','LEFTY2',
         'BMP7','WNT3A','FGF8','ARL13B', 'KRT18', 'KRT19', 'CDH1', 'TIMP3', 'SEMA3C','SEMA3E',
         'SOX2','TBXT','PAX3','CYP26A1','EPHA1','FGF3','FGF4','FGF8','FGF17','FGF19',
         'WNT3A','WNT5A','WNT8A','NEFM',
         'TBXT', 'TBX6','HES7','MSGN1','DLL1','DLL3','LEF1','DKK1','FGF18',
         'MESP1','MESP2','RIPPLY2','CER1','DKK1','DLL1','DLL3', 'LEF1','TCF15','MEOX1','BMP4',
         'UNCX','LFNG','TBX18','FST','SIX1','EYA1','MEOX1','DCHS1','ALDH1A2',
         'BMP5','PAX3','MYF5','CRABP2','RXRG','RARA',
         'PAX3','PAX7','EBF2','TGFB1','TGFBI','FOXC2','SCX','ALDH1A2','TIMP3','RARA',
         'FRZB','PAX1','PAX9','NKX3-2','SOX9','SEMA3C','TGFBI',
         "TWIST1",'SNAI1','FOXC2','WNT5B','GLI1',
         'TBXT','SOX2','HES7','NKX1-2','NEFM','FGF8','FGF17','FGF19',
         'BMP7','WNT3A','WNT5A','WNT8A','SPRY2','SOX2','WNT5B','RARB','CRABP2','HES4','HES5','CCND1','CELSR1',
         'PAX6','IRX3','IRX5','DBX1','DBX2','WNT1','PAX7','NKX6-1','NKX6-2','OLIG2','NKX2-2',
         'FOXA1','FOXA2','SHH','TIMP3','FOXJ1','SEMA3E','DHRS3','RARB',
         'CRABP1','HES4','HES5','DLL1','EBF2','DLL3','ISL1',
         'FGF8','NEFM','CYP26A1','NKX1-2','HES7','WNT5A','WNT8A',
         'CXCR4','CDX2','MSX1','RARG','DHRS3',
         'PAX2','PAX8','OSR1','LHX1','NEFM','DKK1','KRT18','FGF8','RDH10',
         'CD34','KDR','ETV2','SOX17','CLDN5','TIMP3','FOXC2','TGFB1')

tiff(filename = "output/GOI/integrate_adjustment_nodoublets_noAmbien_bbknn.tiff", 
     width = 3000, height = 3300, units = "px")
FeaturePlot(v3, features = unique(GOI), pt.size = 0.2,reduction = "umap.bbknn",
            ncol = 10, order = TRUE)
dev.off()

v3[["RNA"]]=NULL
qs::qsave(v3, file = 'output/01_output/NTSM_downstream_decontX_nodoublets.qs')


##annotation
v3 <- qs::qread(file = 'output/01_output/NTSM_downstream_decontX_nodoublets.qs')

resolutions <- c(1.2,1.4,1.6,1.8,2.0,2.4,2.8,3.0)

v3 <- FindClusters(v3, resolution = resolutions,graph.name = 'bbknn')

p1=DimPlot(v3, reduction = "tsne.bbknn", group.by = paste0("bbknn_res.", resolutions), ncol = 3, label = T) & NoLegend()
p2=DimPlot(v3,reduction = 'umap.bbknn',group.by =  paste0("bbknn_res.", resolutions), ncol = 3, label = T) & NoLegend()

DimPlot(v3,reduction = 'umap.bbknn',group.by =  'bbknn_res.3', label = T) 
Idents(v3) <- v3$bbknn_res.3

v3 <- RenameIdents(v3,
                     '28' = 'Notochord',
                     '31' = 'NMP-Neural',
                     '5' = 'NMP-Meso',
                     '16' = 'Caudal Neural Plate',
                     '12' = 'Early NT',
                     '32' = 'Early FP',
                     '1' = 'Neural Tube',
                     '20' = 'Floor Plate',
                     '0' = 'pPSM',
                     '17' = 'pPSM',
                     '2' = 'pPSM',
                     '25' = 'pPSM',
                     '33' = 'pPSM',
                     '19' = 'pPSM',
                     '29' = 'aPSM',
                     '23' = 'aPSM',
                     '3' = 'E-SM',
                     '10' = 'E-SM',
                     '15' = 'E-SM',
                     '27' = 'E-SM',
                     '8' = 'E-SM',
                     '18' = 'E-SM',
                     '21' = 'M-SM',
                     '7' = 'M-SM',
                     '22' = 'M-SM',
                     '24' = 'Sclerotome',
                     '4' = 'Syndetome',
                     '6' = 'Endotome',
                     '14' = 'Myogenic Progenitors',
                     '11' = 'Myogenic Progenitors',
                     '34' = 'Dermomyotome',
                     '26' = 'Endothelial Cells',
                     '36' = 'Neuroectoderm',
                     '30' = 'Suface Ectoderm',
                     '13' = 'Caudal Meso Progenitors',
                     '9' = 'Intermediate Meso',
                     '35' = 'Neurons')
qs::qsave(v3, file = 'output/02_output/01_NTSM_integrated_decontX.qs')

******************************************************************
*************************02.V3_sublcuster*************************
******************************************************************

v3 <- qs::qread(file = 'output/02_output/01_NTSM_integrated_decontX.qs')
Idents(v3)=v3$seurat_clusters
#### neural tube ####
##anterior_to_posterior
NT_AP_cluster = c('31','16','32','12','1','20')
NT_AP_counts <- v3[['decontX']]@counts[,which(Idents(v3) %in% NT_AP_cluster)]
NT_AP <- CreateSeuratObject(counts = NT_AP_counts, min.cells = 3, min.features = 200, project = 'Neural Tube_AP')

NT_AP <- NormalizeData(NT_AP,normalization.method = 'LogNormalize',scale.factor = 10000,verbose = TRUE)
NT_AP<- FindVariableFeatures(NT_AP, selection.method = "vst", nfeatures = 2000, verbose=FALSE)

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

NT_AP <- CellCycleScoring(NT_AP, s.features = s.genes, g2m.features = g2m.genes)

NT_AP <- ScaleData(NT_AP, vars.to.regress = c("S.Score", "G2M.Score"))
NT_AP <- RunPCA(NT_AP)

#integrate with bbknn
library(bbknnR)
NT_AP <- RunBBKNN(NT_AP, reduction ='pca',run_UMAP = TRUE, run_TSNE = TRUE, UMAP_name ='umap.bbknn',TSNE_name = 'tsne.bbknn',
                  batch_key='orig.ident')

p1 =DimPlot(NT_AP,reduction = 'tsne.bbknn',group.by = 'orig.ident')
p2 =DimPlot(NT_AP,reduction = 'umap.bbknn', group.by = 'orig.ident')
p1+p2


##add metadata
v3_cluster=FetchData(v3,vars = c('bbknn_res.3','celltype','detailed_celltype','detailed_celltype_2nd'))


for (i in 1:nrow(v3_cluster)){
  print(i)
  type=as.character(v3_cluster[i,1])
  if (type %in% NT_AP_cluster){v3_cluster[i,5]='NT_AP'}
  else{v3_cluster[i,5]='Others'}
}
colnames(v3_cluster)=c('cluster','celltype','detailed_celltype','detailed_celltype_2nd','type')
v3_NT_AP=v3_cluster %>% filter(type == 'NT_AP') %>% select(-type)
NT_AP=AddMetaData(NT_AP,metadata = v3_NT_AP)
p1=DimPlot(NT_AP,reduction = 'umap.bbknn',group.by =  'cluster', label = T) & NoLegend()
p2=DimPlot(NT_AP,reduction = 'umap.bbknn',group.by =  'celltype', label = T) & NoLegend()
p3=DimPlot(NT_AP,reduction = 'umap.bbknn',group.by =  'detailed_celltype', label = T) & NoLegend()
p1+p2+p3

qs::qsave(NT_AP, file = 'output/02_output/01_NT_AP.qs')


###NT dorsal to ventral
NT_DV_cluster = c('1','20','32','12','16')
#### neural tube ####
NT_DV_counts <- v3[['decontX']]@counts[,which(Idents(v3) %in% NT_DV_cluster)]
NT_DV <- CreateSeuratObject(counts = NT_DV_counts, min.cells = 3, min.features = 200, project = 'Neural Tube_DV')

NT_DV <- NormalizeData(NT_DV,normalization.method = 'LogNormalize',scale.factor = 10000,verbose = TRUE)
NT_DV <- FindVariableFeatures(NT_DV, selection.method = "vst", nfeatures = 2000, verbose=FALSE)

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

NT_DV <- CellCycleScoring(NT_DV, s.features = s.genes, g2m.features = g2m.genes)

NT_DV <- ScaleData(NT_DV, vars.to.regress = c("S.Score", "G2M.Score"))
NT_DV <- RunPCA(NT_DV)

#integrate with bbknn
library(bbknnR)
NT_DV <- RunBBKNN(NT_DV, reduction ='pca',run_UMAP = TRUE, run_TSNE = TRUE, UMAP_name ='umap.bbknn',TSNE_name = 'tsne.bbknn',
                  batch_key='orig.ident')

p1 =DimPlot(NT_DV,reduction = 'tsne.bbknn',group.by = 'orig.ident')
p2 =DimPlot(NT_DV,reduction = 'umap.bbknn', group.by = 'orig.ident')
p1+p2


##add metadata
v3_cluster=FetchData(v3,vars = c('bbknn_res.3','celltype','detailed_celltype','detailed_celltype_2nd'))

for (i in 1:nrow(v3_cluster)){
  print(i)
  type=as.character(v3_cluster[i,1])
  if (type %in% NT_DV_cluster){v3_cluster[i,5]='NT_DV'}
  else{v3_cluster[i,5]='Others'}
}
colnames(v3_cluster)=c('cluster','celltype','detailed_celltype','detailed_celltype_2nd','type')
v3_NT=v3_cluster %>% filter(type == 'NT_DV') %>% select(-type)
NT_DV=AddMetaData(NT_DV,metadata = v3_NT)
p1=DimPlot(NT_DV,reduction = 'umap.bbknn',group.by =  'cluster', label = T) & NoLegend()
p2=DimPlot(NT_DV,reduction = 'umap.bbknn',group.by =  'celltype', label = T) & NoLegend()
p3=DimPlot(NT_DV,reduction = 'umap.bbknn',group.by =  'detailed_celltype', label = T) & NoLegend()
p1+p2+p3

qs::qsave(NT_DV, file = 'output/02_output/01_NT_DV.qs')
####### Somite ####
###anterior_to_posterior
SM_AP_cluster = c('5','2','25','33','19','0','17','29','23','3','10','15','27','8','18','7','21','22','4','6','14','11','34','24')


SM_AP_counts <- v3[['decontX']]@counts[,which(Idents(v3) %in% SM_AP_cluster)]
SM_AP <- CreateSeuratObject(counts = SM_AP_counts, min.cells = 3, min.features = 200, project = 'Somite_AP')

SM_AP <- NormalizeData(SM_AP,normalization.method = 'LogNormalize',scale.factor = 10000,verbose = TRUE)
SM_AP <- FindVariableFeatures(SM_AP, selection.method = "vst", nfeatures = 2000, verbose=FALSE)

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

SM_AP <- CellCycleScoring(SM_AP, s.features = s.genes, g2m.features = g2m.genes)

SM_AP <- ScaleData(SM_AP, vars.to.regress = c("S.Score", "G2M.Score"))
SM_AP <- RunPCA(SM_AP)

#integrate with bbknn
library(bbknnR)
SM_AP <- RunBBKNN(SM_AP, reduction ='pca',run_UMAP = TRUE, run_TSNE = TRUE, UMAP_name ='umap.bbknn',TSNE_name = 'tsne.bbknn',
                  batch_key='orig.ident')

p1 =DimPlot(SM_AP,reduction = 'tsne.bbknn',group.by = 'orig.ident')
p2 =DimPlot(SM_AP,reduction = 'umap.bbknn', group.by = 'orig.ident')
p1+p2

##add metadata
v3_cluster=FetchData(v3,vars = c('bbknn_res.3','celltype','detailed_celltype','detailed_celltype_2nd'))


for (i in 1:nrow(v3_cluster)){
  print(i)
  type=as.character(v3_cluster[i,1])
  if (type %in% SM_AP_cluster){v3_cluster[i,5]='SM_AP'}
  else{v3_cluster[i,5]='Others'}
}
colnames(v3_cluster)=c('cluster','celltype','detailed_celltype','detailed_celltype_2nd','type')
v3_SM_AP=v3_cluster %>% filter(type == 'SM_AP') %>% select(-type)
SM_AP=AddMetaData(SM_AP,metadata = v3_SM_AP)
p1=DimPlot(SM_AP,reduction = 'umap.bbknn',group.by =  'cluster', label = T) & NoLegend()
p2=DimPlot(SM_AP,reduction = 'umap.bbknn',group.by =  'celltype', label = T) & NoLegend()
p3=DimPlot(SM_AP,reduction = 'umap.bbknn',group.by =  'detailed_celltype', label = T,label.size = 3) & NoLegend()
p1+p2+p3

qs::qsave(SM_AP, file = 'output/02_output/01_SM_AP.qs')

##dorsal to ventral
SM_DV_cluster = c('14','11','34','24','7','21','22','4','6','3','10','15','27','8','18')
#### Somite ####
SM_DV_counts <- v3[['decontX']]@counts[,which(Idents(v3) %in% SM_DV_cluster)]
SM_DV <- CreateSeuratObject(counts = SM_DV_counts, min.cells = 3, min.features = 200, project = 'Somite_DV')

SM_DV <- NormalizeData(SM_DV,normalization.method = 'LogNormalize',scale.factor = 10000,verbose = TRUE)
SM_DV <- FindVariableFeatures(SM_DV, selection.method = "vst", nfeatures = 2000, verbose=FALSE)

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

SM_DV <- CellCycleScoring(SM_DV, s.features = s.genes, g2m.features = g2m.genes)

AP <- c('TBX18','UNCX')
SM_DV <- AddModuleScore(SM_DV, features = list(AP),name = 'AP_enriched')

SM_DV <- ScaleData(SM_DV, vars.to.regress = c("S.Score", "G2M.Score",'AP_enriched1'))
SM_DV <- RunPCA(SM_DV)

#integrate with bbknn
library(bbknnR)
SM_DV <- RunBBKNN(SM_DV, reduction ='pca',run_UMAP = TRUE, run_TSNE = TRUE, UMAP_name ='umap.bbknn',TSNE_name = 'tsne.bbknn',
                  batch_key='orig.ident')

p1 =DimPlot(SM_DV,reduction = 'tsne.bbknn',group.by = 'orig.ident')
p2 =DimPlot(SM_DV,reduction = 'umap.bbknn', group.by = 'orig.ident')
p1+p2


##add metadata
v3_cluster=FetchData(v3,vars = c('bbknn_res.3','celltype','detailed_celltype','detailed_celltype_2nd'))


for (i in 1:nrow(v3_cluster)){
  print(i)
  type=as.character(v3_cluster[i,1])
  if (type %in% SM_DV_cluster){v3_cluster[i,5]='SM_DV'}
  else{v3_cluster[i,5]='Others'}
}
colnames(v3_cluster)=c('cluster','celltype','detailed_celltype','detailed_celltype_2nd','type')
v3_SM=v3_cluster %>% filter(type == 'SM_DV') %>% select(-type)
SM_DV=AddMetaData(SM_DV,metadata = v3_SM)
p1=DimPlot(SM_DV,reduction = 'umap.bbknn',group.by =  'cluster', label = T) & NoLegend()
p2=DimPlot(SM_DV,reduction = 'umap.bbknn',group.by =  'celltype', label = T) & NoLegend()
p3=DimPlot(SM,reduction = 'umap.bbknn',group.by =  'detailed_celltype', label = T) & NoLegend()
p1+p2+p3

qs::qsave(SM_DV, file = 'output/02_output/01_SM_DV.qs')


##noto
Noto_cluster = c('28')
#### noto ####
Noto_counts <- v3[['decontX']]@counts[,which(Idents(v3) %in% Noto_cluster)]
Noto <- CreateSeuratObject(counts = Noto_counts, min.cells = 3, min.features = 200, project = 'Notochord')

Noto <- NormalizeData(Noto,normalization.method = 'LogNormalize',scale.factor = 10000,verbose = TRUE)
Noto<- FindVariableFeatures(Noto, selection.method = "vst", nfeatures = 2000, verbose=FALSE)

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

Noto <- CellCycleScoring(Noto, s.features = s.genes, g2m.features = g2m.genes)

Noto<- ScaleData(Noto, vars.to.regress = c("S.Score", "G2M.Score"))
Noto <- RunPCA(Noto)
#integrate with bbknn
library(bbknnR)
Noto <- RunBBKNN(Noto, reduction ='pca',run_UMAP = TRUE, run_TSNE = TRUE, UMAP_name ='umap.bbknn',TSNE_name = 'tsne.bbknn',
                 batch_key='orig.ident')

Noto <- FindClusters(Noto, resolution = 0.6,graph.name = 'bbknn')
DimPlot(Noto,reduction = 'umap.bbknn',group.by =  "bbknn_res.0.6",label = T) & NoLegend()
ggsave(filename = 'output/02_output/Figure/01_NOTO_resolution.tiff',height = 3.3, width = 4.7,dpi = 300)

Idents(Noto) <- Noto$bbknn_res.0.6
Noto <- RenameIdents(Noto,
                     '0'= 'Noto.4',
                     '1'= 'Noto.2',
                     '2'= 'Noto.3',
                     '3'= 'Noto.1')
Noto$celltype2<- Idents(Noto)
##add metadata
v3_cluster=FetchData(v3,vars = c('bbknn_res.3','celltype','detailed_celltype'))


for (i in 1:nrow(v3_cluster)){
  print(i)
  type=as.character(v3_cluster[i,1])
  if (type %in% Noto_cluster){v3_cluster[i,4]='Noto'}
  else{v3_cluster[i,4]='Others'}
}
colnames(v3_cluster)=c('cluster','celltype','detailed_celltype','type')
v3_Noto=v3_cluster %>% filter(type == 'Noto') %>% select(-type)
Noto=AddMetaData(Noto,metadata = v3_Noto)
p1=DimPlot(Noto,reduction = 'umap.bbknn',group.by =  'seurat_clusters', label = T) & NoLegend()
p2=DimPlot(Noto,reduction = 'umap.bbknn',group.by =  'celltype2', label = T) & NoLegend()
p3=DimPlot(Noto,reduction = 'umap.bbknn',group.by =  'detailed_celltype', label = T) & NoLegend()
p1+p2+p3


##remove high pseudotime cells whihc might be NMP-NEURAL
index <- c('NTSM-V4-D3_AATCTAGCATCACTGC-1', 'NTSM-V4-D3_AATGGGACAAGTGAAG-1',
           'NTSM-V4-D3_ACAACAGGTAACGCTT-1', 'NTSM-V4-D3_ACAATAGAGGCTTAAT-1',
           'NTSM-V4-D3_ACAGCGGGTTCGTATA-1', 'NTSM-V4-D3_ACATTCCTCTATCAGG-1',
           'NTSM-V4-D3_ACCTAACGTTGGATCT-1', 'NTSM-V4-D3_ACCTTCACACTAGGAT-1',
           'NTSM-V4-D3_AGCCCCTTCAGGTTCA-1', 'NTSM-V4-D3_AGCCCTAGTAAGGTAT-1',
           'NTSM-V4-D3_AGGCACGGTGAGACTA-1', 'NTSM-V4-D3_AGGTAATGTTAACCCG-1',
           'NTSM-V4-D3_AGTCGATCACCTTGGT-1', 'NTSM-V4-D3_ATCCTGAAGTACCATG-1',
           'NTSM-V4-D3_ATGCGGTTCCTGAATG-1', 'NTSM-V4-D3_CAAACCATCCCATAGA-1',
           'NTSM-V4-D3_CAAAGCTCATCATCAT-1', 'NTSM-V4-D3_CACAAGTCAGTTGGTG-1',
           'NTSM-V4-D3_CACACAGCAATCTCTG-1', 'NTSM-V4-D3_CAGGGTTTCCATAGGC-1',
           'NTSM-V4-D3_CAGTATTAGCGGACTT-1', 'NTSM-V4-D3_CCAAATTAGAAGCCCT-1',
           'NTSM-V4-D3_CCGTCCACATAAGGTC-1', 'NTSM-V4-D3_CCGTGCAAGGTTCCAA-1',
           'NTSM-V4-D3_CCTACAGGTCAGCCGT-1', 'NTSM-V4-D3_CCTAGCTCATAAGGCT-1',
           'NTSM-V4-D3_CCTATTCAGGTACCGT-1', 'NTSM-V4-D3_CGATGTCAGTATCTGG-1',
           'NTSM-V4-D3_CGCAACTAGTACGCCG-1', 'NTSM-V4-D3_CGCGTTAGTTAGGGTA-1',
           'NTSM-V4-D3_CGGCTTGAGGCGAAGC-1', 'NTSM-V4-D3_CTACGTTTCGTACTAG-1',
           'NTSM-V4-D3_CTCGATCCATAAGTTG-1', 'NTSM-V4-D3_CTCTGATTCAGTACCA-1',
           'NTSM-V4-D3_CTGAACATCCTGAGTA-1', 'NTSM-V4-D3_CTGTACCCAAGTGCCA-1',
           'NTSM-V4-D3_CTTCAGCCATTAGCGG-1', 'NTSM-V4-D3_GCAAGTGAGGATTACT-1',
           'NTSM-V4-D3_GCATGTTGTCTAATCA-1', 'NTSM-V4-D3_GCATTCAAGCCACTCC-1',
           'NTSM-V4-D3_GCCCCAAAGTAGCCGT-1', 'NTSM-V4-D3_GCTCATTTCAGTAGCC-1',
           'NTSM-V4-D3_GCTCTCATCTGTCATC-1', 'NTSM-V4-D3_GGATTAAAGGATCTCA-1',
           'NTSM-V4-D3_GGCGTTACACCTCGTG-1', 'NTSM-V4-D3_GGGCGTAGTAGTGACA-1',
           'NTSM-V4-D3_GGTCAGCCAATGTCGG-1', 'NTSM-V4-D3_GTATCGATCCAGATGC-1',
           'NTSM-V4-D3_GTATTCAGTGGCTCTA-1', 'NTSM-V4-D3_GTCCTTTAGCCGACCA-1',
           'NTSM-V4-D3_GTCTCAAAGTGAGGCA-1', 'NTSM-V4-D3_GTGTAAGGTAGGTTAG-1',
           'NTSM-V4-D3_GTTAAGCCATCCATTG-1', 'NTSM-V4-D3_TAGGCTTGTTTAAGCC-1',
           'NTSM-V4-D3_TATCTTGTCACGCTAA-1', 'NTSM-V4-D3_TCATTCGCAGTTGCAG-1',
           'NTSM-V4-D3_TCCTTAATCGCGCGTA-1', 'NTSM-V4-D3_TCGATTTCACAAGGAA-1',
           'NTSM-V4-D3_TCGCTCGCAGGACATT-1', 'NTSM-V4-D3_TCTCAATCAAGCCGTA-1',
           'NTSM-V4-D3_TGACTAACATTACTAC-1', 'NTSM-V4-D3_TGAGGTTCACATAGCG-1',
           'NTSM-V4-D3_TGAGTCGCACTTACAT-1', 'NTSM-V4-D3_TGGTCCTCACAATCCC-1',
           'NTSM-V4-D3_TGTAAGCCATAGATGC-1', 'NTSM-V4-D4_ATTGACTAGAACGGGC-1',
           'NTSM-V4-D4_GCATAGTCAATTCGCC-1', 'NTSM-V4-D4_TCGGTCAGTGATGGGT-1',
           'NTSM-V4-D5_CCCAGTGAGCTGGACT-1', 'NTSM-V4-D5_TCACTTGCATGCGTCT-1',
           'NTSM-V4-D5_TGTGCAATCCGTGCAG-1', 'NTSM-V4-D6_ACGGGATGTTGAGTAG-1',
           'NTSM-V4-D6_ACTAAGGGTTGCACGC-1', 'NTSM-V4-D6_GGGAACATCCATCGCC-1',
           'NTSM-V4-D7_CTCTCAACACTCATCA-1', 'NTSM-V4-D7_CTGCATTTCACGTGCG-1',
           'NTSM-V4-D7_CTTCAATCATGGCTTC-1',
           'NTSM-V4-D3_AAAGTAGTCCACCAGG-1', 'NTSM-V4-D3_AATGGCGAGTATGGTT-1',
           'NTSM-V4-D3_ACACCCAGTGATTAGG-1', 'NTSM-V4-D3_ACCACACAGCGGTAGA-1',
           'NTSM-V4-D3_ACCCATCCAGGGACTC-1', 'NTSM-V4-D3_ACGCTCAGTAGTCAGA-1',
           'NTSM-V4-D3_ACGGTCACAAGTTCAA-1', 'NTSM-V4-D3_AGGAAATCAGTTGCAG-1',
           'NTSM-V4-D3_AGGAGATTCCTCTCAG-1', 'NTSM-V4-D3_ATCCGACAGGCGAAGC-1',
           'NTSM-V4-D3_ATGCCTGGTAGGTCGG-1', 'NTSM-V4-D3_CAAAGGTAGTGCAATT-1',
           'NTSM-V4-D3_CATATTCCATAACCGG-1', 'NTSM-V4-D3_CCACTAAAGATAGGAC-1',
           'NTSM-V4-D3_CCAGTCCTCACGTGTA-1', 'NTSM-V4-D3_CCGCTAAGTTTGAGGA-1',
           'NTSM-V4-D3_CCGTTGCTCGCTCCAC-1', 'NTSM-V4-D3_CCTACGGAGGGGATAT-1',
           'NTSM-V4-D3_CCTCTCAAGCAGGCTT-1', 'NTSM-V4-D3_GCAACGTTCCTTGAGA-1',
           'NTSM-V4-D3_GCCATAGCATCAAGTA-1', 'NTSM-V4-D3_GCCATCCTCAATTCCC-1',
           'NTSM-V4-D3_GCCCCATTCGTGAATT-1', 'NTSM-V4-D3_GGCATTCAGTGATACT-1',
           'NTSM-V4-D3_GGGCGATAGTTGCGAC-1', 'NTSM-V4-D3_GTAGCAAGTAATGAGT-1',
           'NTSM-V4-D3_GTGTTACGTATTGAGG-1', 'NTSM-V4-D3_GTTGGTTGTCAGGCGA-1',
           'NTSM-V4-D3_TACGTAGTCAATCGGC-1', 'NTSM-V4-D3_TGCAAAGGTTACTGCG-1',
           'NTSM-V4-D3_TGCCGTAGTGCACGCT-1', 'NTSM-V4-D4_ATACACTAGTAGCATT-1',
           'NTSM-V4-D4_CGTACCTCATAAGTGC-1', 'NTSM-V4-D4_CTTACGAAGAATCGAA-1',
           'NTSM-V4-D4_TGTCTCACAAGCCTGA-1', 'NTSM-V4-D6_TGCCGTTCAACATGCG-1',
           'NTSM-V4-D6_TGGTTCCCACATTGAA-1')

Noto2 <- subset(Noto, cells = setdiff(Cells(Noto), index))
qs::qsave(Noto2, file = 'output/02_output/01_Noto.qs')

******************************************************************
****************03.V3_trajectory_pesedotime_analysis**************
******************************************************************

library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(magrittr)

setwd(here::here())
source("R/TA/IO.R")
source("R/TA/preprocess.R")

#### Load data
v3 <- qs::qread("output/02_output/01_NTSM_integrated_decontX.qs")
v3$celltype <- as.character(v3$celltype)
v3$detailed_celltype <- as.character(v3$detailed_celltype)
v3$detailed_celltype_2nd <-as.character(v3$detailed_celltype_2nd)

##trajectory
samples <- list.files("data/loom/")
samples <- samples[1:5]

ldat <- lapply(samples, function(fn) {
  message(glue::glue("Loading {fn} ..."))
  read.loom.matrices(file.path("data/loom/", fn))
})

names(ldat) <- sub(".loom", "", gsub("_", "-", samples))
sapply(ldat[[1]], dim)

matrix.name <- names(ldat[[1]])
matrix.name

ldat.merged <- lapply(matrix.name, function(mn){
  mat.list <- lapply(ldat, function(xx) xx[[mn]])
  do.call(cbind, mat.list) ## merge by columns (cell.ID)
})
names(ldat.merged) <- matrix.name
sapply(ldat.merged, dim)

head(colnames(ldat.merged$spliced))
head(colnames(v3))

fix.cellID <- . %>% gsub("_", "-", .) %>% sub(":", "_", .) %>% sub("x$", "-1", .) %>% sub("cellranger-mergedout-v",'NTSM-V',.) %>% 
  sub("-withBAM",'',.)

fix.cellID(colnames(ldat.merged$spliced)) %>% head()

for (i in seq_along(ldat.merged)) {
  colnames(ldat.merged[[i]]) %<>% fix.cellID()
}

colnames(v3) %in% colnames(ldat.merged$spliced) %>% all() 

for (i in seq_along(ldat.merged)) {
  ldat.merged[[i]] <- ldat.merged[[i]][, colnames(v3)]
}
sapply(ldat.merged, dim)

emat <- ldat.merged$spliced 
nmat <- ldat.merged$unspliced
percent.intron <- colSums(nmat) / (colSums(nmat) + colSums(emat))

clusters <- v3$bbknn_res.3
emat <- filter.genes.by.cluster.expression(emat, clusters, min.max.cluster.average = 0.1)
nmat <- filter.genes.by.cluster.expression(nmat, clusters, min.max.cluster.average = 0.01)

genes.use <- intersect(rownames(emat), rownames(nmat))
length(genes.use)

emat <- emat[genes.use, ]
nmat <- nmat[genes.use, ]

ldat2 <- list(
  spliced = emat,
  unspliced = nmat
)

sapply(ldat2, dim)

v32 <- as.Seurat(x = ldat2)
v32
v32[["decontX"]] <- v32[["spliced"]]
DefaultAssay(v32) <- "decontX"

v32@meta.data <- v3@meta.data

for (rd.name in Reductions(v3)) {
  v32[[rd.name]] <- v3[[rd.name]]
}

v32[["pca"]]@feature.loadings <- matrix()

SeuratDisk::SaveH5Seurat(v32, filename = "output/02_output/03_NTSM_decontX.h5Seurat", overwrite = TRUE)
SeuratDisk::Convert("output/02_output/03_NTSM_decontX.h5Seurat", dest = "h5ad", overwrite = TRUE)


#########pseudotime
#V3
v3[['RNA']] <- v3[['decontX']]
v32 <- v3[rownames(v3[['RNA']]@scale.data),]
DefaultAssay(v32) ='RNA'
v32[['decontX']] <- NULL
sceasy::convertFormat(v32, from = "seurat", to = "anndata", main_layer = "scale.data", 
                      outFile = file.path('output/02_output/01_NTSM_scaled.h5ad'))

##NT
NT_AP <- qs::qread(file = 'output/02_output/01_NT_AP.qs')
NT_AP2 <- NT_AP[rownames(NT_AP[['RNA']]@scale.data),]
sceasy::convertFormat(NT_AP2, from = "seurat", to = "anndata", main_layer = "scale.data", 
                      outFile = file.path('output/02_output/01_NT_AP_scaled.h5ad'))

##SM
SM_AP <- qs::qread(file = 'output/02_output/01_SM_AP.qs')
SM2 <- SM_AP[rownames(SM_AP[['RNA']]@scale.data),]
sceasy::convertFormat(SM2, from = "seurat", to = "anndata", main_layer = "scale.data", 
                      outFile = file.path('output/02_output/01_SM_AP_scaled.h5ad'))

##Noto
Noto <- qs::qread(file = 'output/02_output/01_Noto.qs')
Noto2 <- Noto[rownames(Noto[['RNA']]@scale.data),]
sceasy::convertFormat(Noto2, from = "seurat", to = "anndata", main_layer = "scale.data", 
                      outFile = file.path('output/02_output/01_Noto_scaled.h5ad'))

*********trajectory.ipynb********
{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "66c5ec9b-278b-4682-8716-7d5f454d4c91",
   "metadata": {},
   "outputs": [],
   "source": [
    "import scvelo as scv\n",
    "import scanpy as sc\n",
    "import pandas as pd\n",
    "scv.set_figure_params('scvelo')\n",
    "scv.settings.presenter_view = True"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "e582ae06-2ac2-46e4-ae61-b48b956b3702",
   "metadata": {},
   "source": [
    "# trajectory"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "5c9728db-e6f9-444a-9cc6-ff06c9d0fb7a",
   "metadata": {},
   "outputs": [],
   "source": [
    "adata = scv.read(\"../../output/02_output/03_NTSM_decontX.h5ad\")\n",
    "adata"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "41df2846-49eb-4bfc-9609-16016169ef2a",
   "metadata": {},
   "outputs": [],
   "source": [
    "adata.obs['celltype'] = adata.obs['celltype'].astype('category')\n",
    "scv.pl.proportions(adata, groupby = \"celltype\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "e1b3fb6a-bfc1-4322-a707-41d595af819a",
   "metadata": {},
   "outputs": [],
   "source": [
    "scv.pp.filter_and_normalize(adata, min_shared_counts=50, n_top_genes=2000)\n",
    "scv.pp.moments(adata, n_pcs=30, n_neighbors=30)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "209c5608-7873-466e-8668-773d4e203445",
   "metadata": {},
   "outputs": [],
   "source": [
    "scv.tl.velocity(adata,mode='stochastic')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "553ffa3f-d1a5-4cb5-b43c-593fc9d7f52e",
   "metadata": {},
   "outputs": [],
   "source": [
    "scv.tl.velocity_graph(adata)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "110a0117-1402-4264-9e1d-b4bddd7c3b2b",
   "metadata": {},
   "outputs": [],
   "source": [
    "new_order = [\n",
    "    'Notochord', 'NMP-Neural', 'NMP-Meso', 'Caudal Neural Plate',\n",
    "    'Early NT', 'Early FP', 'Neural Tube', 'Floor Plate',\n",
    "    'Early pPSM.1', 'Early pPSM.2', 'Late pPSM.1', 'Late pPSM.2',\n",
    "    'Late pPSM.3', 'Late pPSM.4', 'Early aPSM', 'Late aPSM',\n",
    "    'E-SM.1', 'E-SM.2', 'E-SM.3', 'E-SM.4', 'E-SM.5', 'E-SM.6',\n",
    "    'M-SM.1', 'M-SM.2', 'M-SM.3', 'Sclerotome', 'Syndetome',\n",
    "    'Endotome', 'Myogenic Progenitors', 'Dermomyotome',\n",
    "    'Endothelial Cells', 'Neuroectoderm', 'Suface Ectoderm',\n",
    "    'Caudal Meso Progenitors', 'Intermediate Meso', 'Neurons'\n",
    "]\n",
    "adata.obs['detailed_celltype'] =pd.Categorical(adata.obs['detailed_celltype'],categories=new_order,ordered=True)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "7b26a2fa-96c1-454e-a141-cc4d1ac3ccb8",
   "metadata": {},
   "outputs": [],
   "source": [
    "scv.pl.velocity_embedding_stream(adata, basis=\"X_umap.bbknn\", color=\"detailed_celltype\",title='hTLO_v3',size=20,fontsize=15,legend_loc= 'right margin',legend_fontsize=10,arrow_size=1,save='../../output/02_output/Figure/04_NTSM_velocity_detailed_celltype.png')"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "a483361e-2955-475c-aa0a-4bfea8eaacb5",
   "metadata": {},
   "source": [
    "### Neural tube DV"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "bfdbd7e6-9be0-468a-a647-0abebeac0cab",
   "metadata": {},
   "outputs": [],
   "source": [
    "adata = scv.read(\"../../output/02_output/03_NT_DV.h5ad\")\n",
    "adata.obs['celltype'] = adata.obs['celltype'].astype('category')\n",
    "scv.pl.proportions(adata, groupby = \"celltype\")\n",
    "scv.pp.filter_and_normalize(adata, min_shared_counts=50, n_top_genes=2000)\n",
    "scv.pp.moments(adata, n_pcs=30, n_neighbors=30)\n",
    "scv.tl.velocity(adata,mode='stochastic')\n",
    "scv.tl.velocity_graph(adata)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "c8eae16e-6267-4567-835d-c7444e173710",
   "metadata": {},
   "outputs": [],
   "source": [
    "scv.pl.velocity_embedding_stream(adata, basis=\"X_umap.bbknn\", color=\"celltype\",title='NT',size=20,fontsize=15,legend_loc= 'right margin',legend_fontsize=10,arrow_size=1,save='../../output/02_output/Figure/04_NT_DV_velocity_celltype.png')"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "64035666-4dc2-4dea-ab17-bebebbb91d5b",
   "metadata": {},
   "source": [
    "### Somite DV "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "32d83bc7-58fa-4a7e-ba38-bbdd4b85a4e1",
   "metadata": {},
   "outputs": [],
   "source": [
    "adata = scv.read(\"../../output/02_output/03_SM_DV.h5ad\")\n",
    "adata.obs['celltype'] = adata.obs['celltype'].astype('category')\n",
    "scv.pl.proportions(adata, groupby = \"celltype\")\n",
    "scv.pp.filter_and_normalize(adata, min_shared_counts=50, n_top_genes=2000)\n",
    "scv.pp.moments(adata, n_pcs=30, n_neighbors=30)\n",
    "scv.tl.velocity(adata,mode='stochastic')\n",
    "scv.tl.velocity_graph(adata)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "a89ad917-1dbe-45da-bd8f-e5d3ead52300",
   "metadata": {},
   "outputs": [],
   "source": [
    "new_order = ['E-SM','M-SM','Sclerotome', 'Syndetome','Endotome','Myogenic Progenitors','Dermomyotome']\n",
    "adata.obs['detailed_celltype_2nd'] =pd.Categorical(adata.obs['detailed_celltype_2nd'],categories=new_order,ordered=True)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "7418dae6-d556-4059-92f9-dc68965b6796",
   "metadata": {},
   "outputs": [],
   "source": [
    "scv.pl.velocity_embedding_stream(adata, basis=\"X_umap.bbknn\", color=\"detailed_celltype_2nd\",title='SM',size=20,fontsize=15,legend_loc= 'right margin',legend_fontsize=10,arrow_size=1,save='../../output/02_output/Figure/04_SM_DV_velocity_detailed_celltype.png')"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "fb884965-fcdf-40fd-8ede-88dc0bf924af",
   "metadata": {},
   "source": [
    "### Noto"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "7a2378a8-7d5a-4704-b407-e6ceee08c932",
   "metadata": {},
   "outputs": [],
   "source": [
    "adata = scv.read(\"../../output/02_output/03_Noto.h5ad\")\n",
    "adata.obs['celltype'] = adata.obs['celltype'].astype('category')\n",
    "scv.pl.proportions(adata, groupby = \"celltype\")\n",
    "scv.pp.filter_and_normalize(adata, min_shared_counts=50, n_top_genes=2000)\n",
    "scv.pp.moments(adata, n_pcs=30, n_neighbors=30)\n",
    "scv.tl.velocity(adata,mode='stochastic')\n",
    "scv.tl.velocity_graph(adata)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "002b6dbd-3727-450d-b8c3-948d4da5b73f",
   "metadata": {},
   "outputs": [],
   "source": [
    "scv.pl.velocity_embedding_stream(adata, basis=\"X_umap.bbknn\", color=\"bbknn_res.0.6\",title='Notochord',size=50,fontsize=15,legend_loc= 'right margin',legend_fontsize=10,arrow_size=1,figsize=(5.7,5.7),save='../../output/02_output/Figure/04_NOTO_velocity.png')"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "velocyto-env",
   "language": "python",
   "name": "velocyto-env"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.7.12"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}

*******pseudotime.ipynb*******
{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "a475cec1-fab6-4bb9-b008-0f708ce82476",
   "metadata": {},
   "outputs": [],
   "source": [
    "import palantir\n",
    "import scanpy as sc\n",
    "import numpy as np\n",
    "import pandas as pd\n",
    "import os\n",
    "\n",
    "# Plotting \n",
    "import matplotlib\n",
    "import matplotlib.pyplot as plt\n",
    "import seaborn as sns\n",
    "import warnings\n",
    "\n",
    "# Inline plotting\n",
    "%matplotlib inline\n",
    "\n",
    "sns.set_style('ticks')\n",
    "matplotlib.rcParams['figure.figsize'] = [4, 4]\n",
    "matplotlib.rcParams['figure.dpi'] = 100\n",
    "matplotlib.rcParams['image.cmap'] = 'Spectral_r'\n",
    "warnings.filterwarnings(action=\"ignore\", module=\"matplotlib\", message=\"findfont\")\n",
    "\n",
    "# Reset random seed\n",
    "np.random.seed(5)"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "cf01f3d3-2af9-49ae-baa3-35f28fa6d9c1",
   "metadata": {},
   "source": [
    "# pseudotime"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "0df7f0d4-3449-4392-b2c4-c4fbd8b9e469",
   "metadata": {},
   "source": [
    "### Neural tube AP"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "fe3b1765-ccd5-4663-84b7-39e154dbe382",
   "metadata": {},
   "outputs": [],
   "source": [
    "adata = sc.read_h5ad(\"../../output/02_output/01_NT_AP_scaled.h5ad\")\n",
    "adata.obs.celltype.value_counts()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "65dda009-50b0-49bf-9786-48eea25f9552",
   "metadata": {},
   "outputs": [],
   "source": [
    "sc.tl.pca(adata, n_comps=50)\n",
    "sc.external.pp.harmony_integrate(adata, key=\"orig.ident\", basis='X_pca')\n",
    "pca_projections = pd.DataFrame(adata.obsm['X_pca_harmony'], index=adata.obs_names)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "9390a47a-f633-45b3-b726-371ccebd1b7c",
   "metadata": {},
   "outputs": [],
   "source": [
    "dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=10)\n",
    "dm_res['EigenVectors'].shape\n",
    "ms_data = palantir.utils.determine_multiscale_space(dm_res)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "df88d22f-6103-4d6f-b6a0-31df9b6d0a79",
   "metadata": {},
   "outputs": [],
   "source": [
    "early_cell = 'NTSM-V4-D4_CCGTTCATCCATGCCC-1'\n",
    "pr_res = palantir.core.run_palantir(ms_data, early_cell=early_cell, use_early_cell_as_start=True, num_waypoints=500, knn=30)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "248586a3-d95b-4a65-a531-0677e29a86df",
   "metadata": {},
   "outputs": [],
   "source": [
    "adata.obs['pseudotime'] = pr_res.pseudotime"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "c2fe3781-d3c4-4d18-adb2-6add9117b62f",
   "metadata": {},
   "outputs": [],
   "source": [
    "sc.pl.embedding(adata, basis=\"X_umap.bbknn\", color=['pseudotime'])"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "934a57f1-e88f-4b7f-82a1-a89fb7ea7a3b",
   "metadata": {},
   "outputs": [],
   "source": [
    "os.chdir('/home/yh001/NTSM_project/NTSM_V4_TA/output/02_output/')\n",
    "os.getcwd()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "0efcf23f-b1e9-478b-b15b-f8bb20c5a306",
   "metadata": {},
   "outputs": [],
   "source": [
    "sc.pl.embedding(adata, basis=\"X_umap.bbknn\", color=['pseudotime'],show=False)\n",
    "plt.savefig(\"06_NT_AP_pseudotime.tiff\", dpi=300, format=\"tiff\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "16340fb9-4579-4d6a-abbe-b168b59470f6",
   "metadata": {},
   "outputs": [],
   "source": [
    "for i in pr_res.branch_probs.columns.tolist():\n",
    "    adata.obs[i] = pr_res.branch_probs[i]\n",
    "sc.pl.embedding(adata, basis=\"X_umap.bbknn\", color=pr_res.branch_probs.columns)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "a2f1db17-d77f-4984-9ca5-a06e408567c0",
   "metadata": {},
   "outputs": [],
   "source": [
    "adata.obs.to_csv(\"../02_output/06_NT_AP_palantir_result.csv\")"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "2d392274-fd30-4baa-887a-ac6624540844",
   "metadata": {},
   "source": [
    "### Somite AP"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "58f4a05d-9fe4-4eb4-9aae-49a9bb3558a4",
   "metadata": {},
   "outputs": [],
   "source": [
    "adata = sc.read_h5ad(\"../../output/02_output/01_SM_AP_scaled.h5ad\")\n",
    "adata.obs.celltype.value_counts()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "9731a653-c24d-40e9-ab0f-89590ba52938",
   "metadata": {},
   "outputs": [],
   "source": [
    "sc.tl.pca(adata, n_comps=50)\n",
    "sc.external.pp.harmony_integrate(adata, key=\"orig.ident\", basis='X_pca')\n",
    "pca_projections = pd.DataFrame(adata.obsm['X_pca_harmony'], index=adata.obs_names)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "d1f196f4-381c-47e1-be3a-e9842b0eb8f1",
   "metadata": {},
   "outputs": [],
   "source": [
    "dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=10)\n",
    "dm_res['EigenVectors'].shape\n",
    "ms_data = palantir.utils.determine_multiscale_space(dm_res)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "310c15c2-37fc-46e4-ae11-22a917a3ba74",
   "metadata": {},
   "outputs": [],
   "source": [
    "early_cell = 'NTSM-V4-D3_AACGGCAAGTTACCGA-1'\n",
    "pr_res = palantir.core.run_palantir(ms_data, early_cell=early_cell, use_early_cell_as_start=True, num_waypoints=500, knn=30)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "36d9b008-8caa-417a-b29c-2fa325a12f8b",
   "metadata": {},
   "outputs": [],
   "source": [
    "adata.obs['pseudotime'] = pr_res.pseudotime"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "37b633c1-73b3-4111-8a28-7fff82c7e7bd",
   "metadata": {},
   "outputs": [],
   "source": [
    "os.chdir('/home/yh001/NTSM_project/NTSM_V4_TA/output/02_output/')\n",
    "os.getcwd()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "726492a2-e46c-4281-87d6-084ab8583fb1",
   "metadata": {},
   "outputs": [],
   "source": [
    "sc.pl.embedding(adata, basis=\"X_umap.bbknn\", color=['pseudotime'],show=False)\n",
    "plt.savefig(\"06_SM_AP_pseudotime.tiff\", dpi=300, format=\"tiff\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "9436af33-f8c8-416f-a058-dfb61fb6b40b",
   "metadata": {},
   "outputs": [],
   "source": [
    "for i in pr_res.branch_probs.columns.tolist():\n",
    "    adata.obs[i] = pr_res.branch_probs[i]\n",
    "sc.pl.embedding(adata, basis=\"X_umap.bbknn\", color=pr_res.branch_probs.columns)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "eac98406-d77a-4865-998d-0636370ad28c",
   "metadata": {},
   "outputs": [],
   "source": [
    "adata.obs.to_csv(\"../02_output/06_SM_AP_palantir_result.csv)"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "dcc20652-d2c2-4be8-a4da-59ddecf1fbc8",
   "metadata": {},
   "source": [
    "### NOTO"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "1a6f45b5-3cbd-41a4-becf-acc74508f466",
   "metadata": {},
   "outputs": [],
   "source": [
    "adata = sc.read_h5ad(\"../../output/02_output/01_Noto_scaled.h5ad\")\n",
    "adata.obs.seurat_clusters.value_counts()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "6f27e256-8c16-4027-ad95-ddb0488c8493",
   "metadata": {},
   "outputs": [],
   "source": [
    "sc.tl.pca(adata, n_comps=50)\n",
    "sc.external.pp.harmony_integrate(adata, key=\"orig.ident\", basis='X_pca')\n",
    "pca_projections = pd.DataFrame(adata.obsm['X_pca_harmony'], index=adata.obs_names)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "8357b65b-3991-48c3-af16-0b94fb449ccd",
   "metadata": {},
   "outputs": [],
   "source": [
    "dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=10)\n",
    "dm_res['EigenVectors'].shape\n",
    "ms_data = palantir.utils.determine_multiscale_space(dm_res)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "0f5170fe-7fc8-4e2f-a52a-87d59664ca39",
   "metadata": {},
   "outputs": [],
   "source": [
    "early_cell = 'NTSM-V4-D3_GTTGCCTCAAGTAAGG-1'\n",
    "pr_res = palantir.core.run_palantir(ms_data, early_cell=early_cell, use_early_cell_as_start=True, num_waypoints=500, knn=30)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "e0c90dde-7aa6-4e3a-9b2f-d8e2658353ab",
   "metadata": {},
   "outputs": [],
   "source": [
    "adata.obs['pseudotime'] = pr_res.pseudotime\n",
    "adata.obs[\"DP\"] = pr_res.entropy"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "7c1671c5-fdec-4a07-8550-c7a6d130c127",
   "metadata": {},
   "outputs": [],
   "source": [
    "sc.pl.embedding(adata, basis=\"X_umap.bbknn\", color=['pseudotime', 'DP'])"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "b148305c-cf13-448f-9568-b23c20cdddb0",
   "metadata": {},
   "outputs": [],
   "source": [
    "os.chdir('/home/yh001/NTSM_project/NTSM_V4_TA/output/02_output/')\n",
    "os.getcwd()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "e8504af1-81ae-4376-8f7e-2c39a9740a54",
   "metadata": {},
   "outputs": [],
   "source": [
    "sc.pl.embedding(adata, basis=\"X_umap.bbknn\", color=['pseudotime'],save='06_NOTO_pesudotime.png')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "7b5c98a6-6339-427d-96f2-402b90f43ec7",
   "metadata": {},
   "outputs": [],
   "source": [
    "for i in pr_res.branch_probs.columns.tolist():\n",
    "    adata.obs[i] = pr_res.branch_probs[i]"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "28c7c482-4411-4ba0-ad42-67ad49533fc8",
   "metadata": {},
   "outputs": [],
   "source": [
    "adata.obs.to_csv('06_Noto_palantir_result.csv')"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "palantir-env",
   "language": "python",
   "name": "palantir-env"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.9.19"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}

******************************************************************
***********************04.V3_scenic_analysis**********************
******************************************************************

library(Seurat)
library(tidydr)
library(reticulate)
library(Matrix)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

library(SCENIC)
library(AUCell)
library(RcisTarget)
library(SCopeLoomR)
library(KernSmooth)
library(BiocParallel)
library(ggplot2)
library(data.table)
library(grid)

srt <- readRDS("finalized_hTLO/01_NTSM_integrated_decontX_annotated_third_rctd.rds")

Idents(combined.sct) <- "detailed_celltype"
exprMat <- as.matrix(GetAssayData(combined.sct, assay = "decontX", slot = "data"))
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
dbDir="ref/SCENIC/cisTarget_databases"
dbName=c("hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather","hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
scenicOptions <- initializeScenic(org = "hgnc",dbDir = dbDir, dbs=dbName, nCores = 10) 
# saveRDS(scenicOptions, file="scenic/savedobj/checkpoint_1_scenicOptions.rds") 

# Rename TBXT to T to match the annotation in SCENIC database
if ("TBXT" %in% rownames(exprMat)) {
    rownames(exprMat)[rownames(exprMat) == "TBXT"] <- "T"
}

# Filter expression matrix
minCountsPerGene <- 40
minSamples <- 30

genesKept <- geneFiltering(exprMat, scenicOptions, minCountsPerGene = minCountsPerGene, minSamples = minSamples)
length(genesKept)
exprMat_filtered <- exprMat[genesKept, ]
# saveRDS(exprMat_filtered, file = "scenic/savedobj/checkpoint_1_exprMat_filtered.rds")

# Run correlation
runCorrelation(exprMat_filtered, scenicOptions)

##### RUN GRNBoost in  4.1_scenic_grnboost.py #####

# Read GRNBoost output and convert to GENIE3 format
GRNBoost_output <- read.delim("scenic/savedobj/intermediate_for_grnboost/output/grn_output.tsv", header=FALSE)
colnames(GRNBoost_output) <- c("TF","Target","weight")
saveRDS(GRNBoost_output, file="int/1.4_GENIE3_linkList.Rds")

# Main scenic computation
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

saveRDS(scenicOptions, file="scenic/savedobj/checkpoint_3_scenicOptions.rds") # To save status
saveRDS(cellInfo, file=getDatasetInfo(scenicOptions, "cellInfo"))


#### Get AUC results ####
cellInfo <- readRDS("scenic/script/int/cellInfo.Rds")
exprMat_filtered <- readRDS("scenic/savedobj/checkpoint_1_exprMat_filtered.rds")
scenicOptions <- readRDS("scenic/script/int/scenicOptions.Rds")

cellInfo['detailed_celltype'] <- combined.sct@meta.data[rownames(cellInfo), 'detailed_celltype']

library(SCopeLoomR)
scenicLoomPath <- getOutName(scenicOptions, "loomFile")
loom <- open_loom(scenicLoomPath)

# Read information from loom file:
regulons_incidMat <- get_regulons(loom)
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom)
regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)

# AUC to umap
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
auc_df = as.data.frame(aucell_regulonAUC@assays@data$AUC)
auc_df <- t(auc_df)

combined.sct@meta.data <- cbind(combined.sct@meta.data, auc_df[rownames(combined.sct@meta.data), ])

# Visualize AUC in umap
pdf('scenic/script/output/regulonAUC_featureplot.pdf', width=10, height=10)
for (regulon in colnames(auc_df)){
   print(FeaturePlot(combined.sct, regulon, order=T, reduction='umap.bbknn', raster = F))
}
dev.off()


*****scenic grnboost.py*******

import pandas as pd
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

in_file  = 'scenic/savedobj/intermediate_for_grnboost/1.1_exprMatrix_filtered_t.txt'
tf_file  = 'scenic/savedobj/intermediate_for_grnboost/1.1_inputTFs.txt'

ex_matrix = pd.read_csv(in_file, sep='\t')
tf_names = load_tf_names(tf_file)

# infer the gene regulatory network
network = grnboost2(expression_data=ex_matrix,
                    tf_names=tf_names)

network.head()

********calculate_rss********

library(Seurat)
library(tidydr)
library(reticulate)
library(Matrix)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

library(SCENIC)
library(AUCell)
library(RcisTarget)
library(SCopeLoomR)
library(KernSmooth)
library(BiocParallel)
library(ggplot2)
library(data.table)
library(grid)

srt <- readRDS("finalized_hTLO/01_NTSM_integrated_decontX_annotated_third_rctd.rds")

celltype_to_detailed_cluster <- c(
  "28" = "Notochord",
  "31" = "NMP-Neural",
  "5" = "NMP-Meso",
  "16" = "Caudal Neural Plate",
  "12" = "Early NT",
  "32" = "Early FP",
  "1" = "Neural Tube",
  "20" = "Floor Plate",
  "0" = "pPSM",
  "17" = "pPSM",
  "2" = "pPSM",
  "25" = "pPSM",
  "33" = "pPSM",
  "19" = "pPSM",
  "29" = "aPSM",
  "23" = "aPSM",
  "3" = "E-SM",
  "10" = "E-SM",
  "15" = "E-SM",
  "27" = "E-SM",
  "8" = "E-SM",
  "18" = "E-SM",
  "21" = "M-SM",
  "7" = "M-SM",
  "22" = "M-SM",
  "24" = "Sclerotome",
  "4" = "Syndetome",
  "6" = "Endotome",
  "14" = "Myogenic Progenitors",
  "11" = "Myogenic Progenitors",
  "34" = "Dermomyotome",
  "26" = "Endothelial Cells",
  "36" = "Neuroectoderm",
  "30" = "Suface Ectoderm",
  "13" = "Caudal Meso Progenitors",
  "9" = "Intermediate Meso",
  "35" = "Neurons"
)
combined.sct@meta.data['detailed_celltype_2'] <- celltype_to_detailed_cluster[combined.sct@meta.data$seurat_clusters]

Idents(combined.sct)="detailed_celltype_2"

exprMat <- as.matrix(GetAssayData(combined.sct, assay = "decontX", slot = "data"))
cellInfo <- data.frame(seuratCluster=Idents(combined.sct))
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
dbDir="ref/SCENIC/cisTarget_databases"
dbName=c("hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather","hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
scenicOptions <- initializeScenic(org = "hgnc",dbDir = dbDir, dbs=dbName, nCores = 10) 

cellInfo <- readRDS("scenic/script/int/cellInfo.Rds")
exprMat_filtered <- readRDS("scenic/savedobj/checkpoint_2_exprMat_filtered.rds")
scenicOptions <- readRDS("scenic/script/int/scenicOptions.Rds")

library(SCopeLoomR)
scenicLoomPath <- getOutName(scenicOptions, "loomFile")
loom <- open_loom(scenicLoomPath)

# Read information from loom file:
regulons_incidMat <- get_regulons(loom)
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom)
regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)

## Internal functions:
.H <- function(pVect){
  pVect <- pVect[pVect>0] # /sum(pVect) ??
  - sum(pVect * log2(pVect))
}

# Jensen-Shannon Divergence (JSD)
calcJSD <- function(pRegulon, pCellType)
{
  (.H((pRegulon+pCellType)/2)) - ((.H(pRegulon)+.H(pCellType))/2)
}

# Regulon specificity score (RSS)
.calcRSS.oneRegulon <- function(pRegulon, pCellType)
{
  jsd <- calcJSD(pRegulon, pCellType)
  1 - sqrt(jsd)
}

AUC <- getAUC(regulonsAUC)

cellAnnotation <- cellInfo %>% select(detailed_celltype)
cellAnnotation <- cellAnnotation[colnames(AUC), 'detailed_celltype']
cellTypes <- unique(cellAnnotation)

normAUC <- AUC/rowSums(AUC)
# 
ctapply <- lapply
rss <- ctapply(cellTypes, function(thisType)
sapply(rownames(normAUC), function(thisRegulon)
    {
        pRegulon <- normAUC[thisRegulon,]
        pCellType <- as.numeric(cellAnnotation==thisType)
        pCellType <- pCellType/sum(pCellType)
        .calcRSS.oneRegulon(pRegulon, pCellType)
    })
)
rss <- do.call(cbind, rss)
colnames(rss) <- cellTypes
# return(rss)

saveRDS(rss, 'scenic/savedobj/scenic_rss.rds')



