********************************************************************
********************************V2_plot*****************************
********************************************************************

********************************************
******************feature_plot**************
********************************************

library(Seurat)
library(SeuratObject)
library(ggplot2)
library(reshape2)
T_39 <- readRDS(file = 'data/spatial/v2_39L_rctd-1_annotated.rds')
T_65 <- readRDS(file = 'data/spatial/v2_65T_rctd-1_annotated.rds')

Idents(T_65) <-T_65$Spatial.008um_snn_res.1.4
T_65 <- RenameIdents(T_65,
                     '2' = 'Tail bud progenitors.1',
                     '10' = 'Tail bud progenitors.2',
                     '4' = 'NMP-Neural',
                     '9' = 'Node-like',
                     '7' = 'Neural Lineage.1',
                     '5' = 'Neural Lineage.2',
                     '8' = 'Neural Lineage.3 (Dorsal)',
                     '3' = 'Neural Lineage.4 (Ventral)',
                     '6' = 'Somitic Lineage.1 (Dorsal)',
                     '0' = 'Somitic Lineage.2 (Ventral)',
                     '1' = 'Somitic Lineage.3',
                     '12' = 'Endothelial Cells',
                     '11' = 'Mixed Meso')
T_65$detailed_celltype <- Idents(T_65)
T65_Cell <- setNames(
  c("#FF1677", "#FF1677", "#FAA31C", "#4F474D", "#16FEC8", "#16FEC8", 
    "#D8EA60", "#D58A81", "#0D8AFF", "#009285", "#A0EBF1", "#9A8116", 
    "#B6CAFA"),
  c("2", "10", "4", "9", "7", "5", "8", "6", "0", "1", "3", "12", "11")
)
DimPlot(T_65, group.by = 'Spatial.008um_snn_res.1.4',label = T) +scale_color_manual(values = T65_Cell)
ggsave(filename = 'output/03_output_spatial/02_T65_dimplot.tiff',height =6.1,width = 7.5,dpi = 300 )


Idents(T_39) <- T_39$Spatial.008um_snn_res.0.8
T_39 <- RenameIdents(T_39,
                     '10' = 'NMP-Neural / Tail bud progenitors',
                     '7' = 'Neural Lineage.1 (Caud. NP)',
                     '2' = 'Neural Lineage.3 (Dorsal)',
                     '14' = 'Neural Lineage.4 (Dorsal)',
                     '9' = 'Neural Lineage.5 (Ventral)',
                     '13' = 'Neural Lineage.2',
                     '11' = 'PSM',
                     '4' = 'Somitic Lineage.1 (Dorsal)',
                     '3' = 'Somitic Lineage.2',
                     '5' = 'Somitic Lineage.3 (Ventral)',
                     '1' = 'Somitic Lineage.4 (Ventral)',
                     '0' = 'Somitic Lineage.5 (Ventral)',
                     '6' = 'Somitic Lineage.5 (Syndetome)',
                     '8' = 'Somitic Lineage.6 (Endotome)',
                     '12' = 'Endothelial Cells')


T39_Cell <- setNames(
  c("#FF1677", "#16FEC8", "#16FEC8", "#16FEC8", "#D58A81", 
    "#16FEC8", "#9F16FB", "#0D8AFF", "#0D8AFF", "#009285", 
    "#009285", "#009285", "#FAA31C", "#F96ADB", "#9A8116"),
  c("10", "7", "2", "14", "9","13", "11", "4", "3", "5","1", "0", "6", "8", "12"))
DimPlot(T_39, group.by = 'Spatial.008um_snn_res.0.8',label = T) +scale_color_manual(values = T39_Cell)
ggsave(filename = 'output/03_output_spatial/02_T39_dimplot.pdf',height =6.1,width = 7.5,dpi = 300 )

********************************************
********************other_plot**************
********************************************
# Seurat4.3
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(reshape2)
library(future)
library(qs)
library(ComplexHeatmap)
library(scCustomize)
library(ggplot2)
library(tidyr)

srt <- readRDS('data/02_v2_integration.rds')

srt <- FindClusters(srt, resolution = 1.6, graph.name = 'bbknn')

# Annotate
cluster_to_celltype <- c(
  "10" = "aPSM",
  "0" = "E-SM",
  "1" = "E-SM",
  "5" = "E-SM",
  "6" = "E-SM",
  "13" = "E-SM",
  "7" = "E-SM",
  "21" = "EC",
  "11" = "Ventral NT/FP",
  "3" = "Dorsal NT",
  "9" = "L-SM",
  "19" = "L-SM",
  "4" = "M-SM",
  "14" = "M-SM",
  "16" = "M-SM",
  "23" = "M-SM",
  "17" = "NMP-Meso/pPSM",
  "8" = "NMP-Neural",
  "2" = "Caud. Meso",
  "24" = "Node-like",
  "22" = "Neuron",
  "12" = "SAG responsive, neural",
  "15" = "SAG responsive, neural",
  "20" = "SAG responsive, somitic",
  "18" = "Caud. NP"
)
srt@meta.data$celltype <- cluster_to_celltype[as.character(srt@meta.data$seurat_clusters)]

# Dotplot
Idents(srt) <- 'celltype'
Idents(srt) <- factor(Idents(srt), levels = c(
    'SAG responsive, somitic',
    'Ventral NT/FP',
    'Dorsal NT',
    'L-SM',
    'M-SM',
    'NMP-Meso/pPSM',
    'E-SM',
    'aPSM',
    'Caud. Meso',
    'Caud. NP',
    'NMP-Neural',
    'SAG responsive, neural',
    'EC',
    'Neuron',
    'Node-like'
))

# Define the gene list as an R character vector
gene_list <- c(
  'LEF1', 'TCF15', 'MEOX1', 'PAX3', 'TBX18', 'UNCX', 'MEIS1', 'EBF2', 'PAX7',  
  'TWIST1', 'SCX', 'PAX1', 'PAX9', 'NKX3-2', 'COL2A1',
  'SOX2', 'PAX6', 'IRX3', 'IRX5', 'NKX6-2', 'OLIG2', 'OLIG3', 'NKX6-1', 'GLI1', 
  'FOXA1', 'FOXA2', 'SHH'
)

options(repr.plot.width = 13, repr.plot.height = 2.8)
tiff('v2_figures/htlov2_dotplot.tiff', width = 13, height = 2.8, res = 300, units = 'in')
print(DotPlot(srt, features = gene_list, idents = c('E-SM', 'M-SM', 'L-SM', 'Dorsal NT', 'Ventral NT/FP')) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  geom_point(aes(fill = ifelse(pct.exp < 5, 'white', NA), alpha = ifelse(pct.exp < 5, 1, 0)), 
             size = 9, shape = 21, color = 'white', na.rm = TRUE) +
  scale_fill_manual(values = c('white' = "white", 'FALSE' = NA), 
                    na.value = "white", guide = "none") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")), alpha = "none") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.box = "horizontal") + 
  labs(size='% Expr.', color = 'Avg. Expr.') + 
  scale_size(
    range = c(1, 8),
    breaks = c(10, 20, 30, 40),
    labels = c("10", "20", "30", '40')
  ))
dev.off()

# Cell counts
df1 <- data.frame()
for (d in c('v2-D5', 'v2-D7')){
    for (c in c("E-SM", "M-SM", "L-SM")){
        genes_of_interest <- c("PAX3", "MYF5", "PAX1")
        data <- FetchData(srt, vars = genes_of_interest)
        data <- cbind(data, srt@meta.data[row.names(data), c('orig.ident', 'celltype')])
        data <- data %>% filter(celltype == c, orig.ident == d) %>% select(genes_of_interest)

        # Step 1: Convert expression values to binary
        binary_data <- data > 0  # Drop the 'cell' column and check for values > 0
        binary_data <- as.data.frame(binary_data)  # Ensure it's a data frame

        # Step 2: Create a binary code for each combination
        binary_data$combination <- apply(binary_data, 1, function(row) paste(as.integer(row), collapse = ""))

        # Step 3: Count the number of cells expressing each combination
        combination_counts <- binary_data %>%
            group_by(combination) %>%
            summarize(count = n())
        combination_counts <- combination_counts %>%
            complete(combination = c("100","010","001","110","101","011", "000","111"), 
                fill = list(value = 0))
        combination_counts['group'] <- d
        combination_counts['celltype'] <- c
        df1 <- rbind(df1, combination_counts)

    }
}
df1[is.na(df1$count), "count"] <- 0
df1$combination <- factor(df1$combination, levels = c("100","010","001","110","101","011", "000","111"))
df1$celltype <- factor(df1$celltype, levels = c("E-SM", "M-SM", "L-SM"))
df1 <- df1 %>% arrange(group, celltype, combination)
df1$condition <- paste0(df1$group, df1$celltype, df1$combination)
df1$condition <- factor(df1$condition, levels = unique(df1$condition))
df1 <- df1 %>% data.frame()

list_values <- c(
  "v2-D5E-SM000" = "PAX3-MYF5-PAX1-",
  "v2-D5E-SM100" = "PAX3+",
  "v2-D5E-SM010" = "MYF5+",
  "v2-D5E-SM001" = "PAX1+",
  "v2-D5E-SM110" = "PAX3+MYF5+",
  "v2-D5E-SM101" = "PAX3+PAX1+",
  "v2-D5E-SM011" = "MYF5+PAX1+",
  "v2-D5E-SM111" = "",
  "v2-D5M-SM000" = "PAX3-MYF5-PAX1-",
  "v2-D5M-SM100" = "PAX3+",
  "v2-D5M-SM010" = "MYF5+",
  "v2-D5M-SM001" = "PAX1+",
  "v2-D5M-SM110" = "PAX3+MYF5+",
  "v2-D5M-SM101" = "PAX3+PAX1+",
  "v2-D5M-SM011" = "MYF5+PAX1+",
  "v2-D5M-SM111" = "",
  "v2-D5L-SM000" = "PAX3-MYF5-PAX1-",
  "v2-D5L-SM100" = "PAX3+",
  "v2-D5L-SM010" = "MYF5+",
  "v2-D5L-SM001" = "PAX1+",
  "v2-D5L-SM110" = "PAX3+MYF5+",
  "v2-D5L-SM101" = "PAX3+PAX1+",
  "v2-D5L-SM011" = "MYF5+PAX1+",
  "v2-D5L-SM111" = "",
  "v2-D7E-SM000" = "PAX3-MYF5-PAX1-",
  "v2-D7E-SM100" = "PAX3+",
  "v2-D7E-SM010" = "MYF5+",
  "v2-D7E-SM001" = "PAX1+",
  "v2-D7E-SM110" = "PAX3+MYF5+",
  "v2-D7E-SM101" = "PAX3+PAX1+",
  "v2-D7E-SM011" = "MYF5+PAX1+",
  "v2-D7E-SM111" = "",
  "v2-D7M-SM000" = "PAX3-MYF5-PAX1-",
  "v2-D7M-SM100" = "PAX3+",
  "v2-D7M-SM010" = "MYF5+",
  "v2-D7M-SM001" = "PAX1+",
  "v2-D7M-SM110" = "PAX3+MYF5+",
  "v2-D7M-SM101" = "PAX3+PAX1+",
  "v2-D7M-SM011" = "MYF5+PAX1+",
  "v2-D7M-SM111" = "",
  "v2-D7L-SM000" = "PAX3-MYF5-PAX1-",
  "v2-D7L-SM100" = "PAX3+",
  "v2-D7L-SM010" = "MYF5+",
  "v2-D7L-SM001" = "PAX1+",
  "v2-D7L-SM110" = "PAX3+MYF5+",
  "v2-D7L-SM101" = "PAX3+PAX1+",
  "v2-D7L-SM011" = "MYF5+PAX1+",
  "v2-D7L-SM111" = ""
)

color_mapping <- c(
  "v2-D5E-SM000" = "#4F474D",
  "v2-D5E-SM001" = "#E2E3E1",
  "v2-D5E-SM010" = "#FD1616",
  "v2-D5E-SM011" = "#FD1CFF",
  "v2-D5E-SM100" = "#00FB0D",
  "v2-D5E-SM101" = "#0D8AFF",
  "v2-D5E-SM110" = "#FAA31C",
  "v2-D5E-SM111" = "white",
  "v2-D5M-SM000" = "#4F474D",
  "v2-D5M-SM001" = "#E2E3E1",
  "v2-D5M-SM010" = "#FD1616",
  "v2-D5M-SM011" = "#FD1CFF",
  "v2-D5M-SM100" = "#00FB0D",
  "v2-D5M-SM101" = "#0D8AFF",
  "v2-D5M-SM110" = "#FAA31C",
  "v2-D5M-SM111" = "white",
  "v2-D5L-SM000" = "#4F474D",
  "v2-D5L-SM001" = "#E2E3E1",
  "v2-D5L-SM010" = "#FD1616",
  "v2-D5L-SM011" = "#FD1CFF",
  "v2-D5L-SM100" = "#00FB0D",
  "v2-D5L-SM101" = "#0D8AFF",
  "v2-D5L-SM110" = "#FAA31C",
  "v2-D5L-SM111" = "white",
  "v2-D7E-SM000" = "#4F474D",
  "v2-D7E-SM001" = "#E2E3E1",
  "v2-D7E-SM010" = "#FD1616",
  "v2-D7E-SM011" = "#FD1CFF",
  "v2-D7E-SM100" = "#00FB0D",
  "v2-D7E-SM101" = "#0D8AFF",
  "v2-D7E-SM110" = "#FAA31C",
  "v2-D7E-SM111" = "white",
  "v2-D7M-SM000" = "#4F474D",
  "v2-D7M-SM001" = "#E2E3E1",
  "v2-D7M-SM010" = "#FD1616",
  "v2-D7M-SM011" = "#FD1CFF",
  "v2-D7M-SM100" = "#00FB0D",
  "v2-D7M-SM101" = "#0D8AFF",
  "v2-D7M-SM110" = "#FAA31C",
  "v2-D7M-SM111" = "white",
  "v2-D7L-SM000" = "#4F474D",
  "v2-D7L-SM001" = "#E2E3E1",
  "v2-D7L-SM010" = "#FD1616",
  "v2-D7L-SM011" = "#FD1CFF",
  "v2-D7L-SM100" = "#00FB0D",
  "v2-D7L-SM101" = "#0D8AFF",
  "v2-D7L-SM110" = "#FAA31C",
  "v2-D7L-SM111" = "white"
)

options(repr.plot.width = 9, repr.plot.height = 4)
df1 %>% ggplot(aes(x = condition, y = count, fill = condition)) + 
    geom_col(color='black') + 
    scale_x_discrete(labels = list_values) +
    coord_cartesian(ylim = c(135, 2850)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none',
        panel.background = element_blank(),  # Remove background
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),  # Remove minor gridlines
        axis.line = element_line(color = "black"),  # Keep axis lines
        axis.line.y.right = element_blank(),  # Remove right axis line
        axis.line.x.top = element_blank(),    # Remove top axis line
        axis.ticks.y.right = element_blank(), # Remove right axis ticks
        axis.ticks.x.top = element_blank(),   # Remove top axis ticks
        axis.text.y.right = element_blank(),  # Remove right axis text
        axis.text.x.top = element_blank()) +     # Remove top axis text) +
    geom_segment(aes(x = 0.6, xend = 7.4, y = 2450, yend = 2450), size = 0.5) +
    annotate("text", x = 4, y = 2550, label = "E-SM", size = 4) +
    geom_segment(aes(x = 8.6, xend = 15.4, y = 2450, yend = 2450), size = 0.5) +
    annotate("text", x = 12, y = 2550, label = "M-SM", size = 4) +
    geom_segment(aes(x = 16.6, xend = 23.4, y = 2450, yend = 2450), size = 0.5) +
    annotate("text", x = 20, y = 2550, label = "L-SM", size = 4) +
    geom_segment(aes(x = 24.6, xend = 31.4, y = 2450, yend = 2450), size = 0.5) +
    annotate("text", x = 28, y = 2550, label = "E-SM", size = 4) +
    geom_segment(aes(x = 32.6, xend = 39.4, y = 2450, yend = 2450), size = 0.5) +
    annotate("text", x = 36, y = 2550, label = "M-SM", size = 4) +
    geom_segment(aes(x = 40.6, xend = 47.4, y = 2450, yend = 2450), size = 0.5) +
    annotate("text", x = 44, y = 2550, label = "L-SM", size = 4) +

    geom_segment(aes(x = 0.6, xend = 23.4, y = 2700, yend = 2700), size = 0.5) +
    annotate("text", x = 12, y = 2800, label = "v2-D5", size = 4) +
    geom_segment(aes(x = 24.6, xend = 47.4, y = 2700, yend = 2700), size = 0.5) +
    annotate("text", x = 36, y = 2800, label = "v2-D7", size = 4) +
    scale_fill_manual(values = color_mapping) + 
    labs(y = 'Count', x = ',') -> p1

p1

tiff('v2_figures/4_1_cellcount_bar.tiff', width = 9, height = 4, res = 300, units = 'in')
print(p1)
dev.off()

# Gene positive cell density
genes <- c("PAX3", "PAX1", "PAX9", "TWIST1", "NKX3-2", "EBF2", "KDR", "LOX", 
           "MYF5", "MYF6", "SCX", "MEF2C", "PAX6", "WNT1", "DBX1", "DBX2", 
           "NKX6-1", "NKX6-2", "FOXA1", "FOXA2", "SHH", "OLIG2", "NKX2-2")

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

embs <- Embeddings(srt, reduction = 'umap.harmony')
gene_expression <- FetchData(srt, vars = genes)

embs <- cbind(embs, srt@meta.data[row.names(embs), c('orig.ident', 'celltype')])
embs <- cbind(embs, gene_expression[rownames(embs), ])

number <- 1
for (gene in genes){
    label <- paste0(gene,'+')

    plot_embs <- embs[c('UMAP_1','UMAP_2','orig.ident','celltype',gene)]
    colnames(plot_embs) <- c('UMAP_1','UMAP_2','orig.ident','celltype','expr')
    plot_embs['expressing'] <- "NA"
    plot_embs[(plot_embs$expr > 0), 'expressing'] <- label

    plot.list <- c()
    n <- 1

    color_min <- 10000
    color_max <- 0
    for (s in unique(plot_embs$orig.ident)){
        target_embs <- plot_embs %>% filter(expressing == label, orig.ident == s)
        if (nrow(target_embs) > 0){
            s <- get_density(target_embs$UMAP_1, target_embs$UMAP_2, n = 100)

            color_min <- min(color_min, min(s))
            color_max <- max(color_max, max(s))
        }
    }

    for (s in unique(plot_embs$orig.ident)){
        target_embs <- plot_embs %>% filter(expressing == label, orig.ident == s)
        bg_embs <- plot_embs %>% filter(expressing != label, orig.ident == s)
        if (nrow(target_embs) > 0){
            target_embs$density <- get_density(target_embs$UMAP_1, target_embs$UMAP_2, n = 100)
            if (n == 1){
                target_embs %>% ggplot( aes(x = UMAP_1, y = UMAP_2, color=density)) + 
                    geom_point(aes(x = UMAP_1, y = UMAP_2), color='lightgray', data = bg_embs, size = .5) +
                    geom_point(size = .5) +
                    scale_color_gradient2(low='lightgray', mid='red', high='yellow', limits = c(color_min, color_max)) + 
                    theme_minimal() +
                    theme(
                        panel.background = element_blank(),  # Remove panel background
                        panel.grid.major = element_blank(),  # Remove major gridlines
                        panel.grid.minor = element_blank(),  # Remove minor gridlines
                        plot.background = element_blank(),   # Remove plot background
                        axis.ticks = element_blank(),        # Remove axis ticks
                        axis.text = element_blank(),
                        legend.position = "none"
                    ) +
                    labs(y = label, x = s, fill=',') -> p1
            } else if (n == 2){
                target_embs %>% ggplot( aes(x = UMAP_1, y = UMAP_2, color=density)) + 
                    geom_point(aes(x = UMAP_1, y = UMAP_2), color='lightgray', data = bg_embs, size = .5) +
                    geom_point(size = .5) +
                    scale_color_gradient2(low='lightgray', mid='red', high='yellow', limits = c(color_min, color_max)) + 
                    theme_minimal() +
                    theme(
                        panel.background = element_blank(),  # Remove panel background
                        panel.grid.major = element_blank(),  # Remove major gridlines
                        panel.grid.minor = element_blank(),  # Remove minor gridlines
                        plot.background = element_blank(),   # Remove plot background
                        axis.ticks = element_blank(),        # Remove axis ticks
                        axis.text = element_blank(),
                        axis.title.y = element_blank(),
                        legend.position = "none"
                    ) +
                    labs(y = label, x = s, fill=',') -> p1
            } else if (n == 3){
                target_embs %>% ggplot( aes(x = UMAP_1, y = UMAP_2, color=density)) + 
                    geom_point(aes(x = UMAP_1, y = UMAP_2), color='lightgray', data = bg_embs, size = .5) +
                    geom_point(size = .5) +
                    scale_color_gradient2(low='lightgray', mid='red', high='yellow', limits = c(color_min, color_max)) + 
                    theme_minimal() +
                    theme(
                        panel.background = element_blank(),  # Remove panel background
                        panel.grid.major = element_blank(),  # Remove major gridlines
                        panel.grid.minor = element_blank(),  # Remove minor gridlines
                        plot.background = element_blank(),   # Remove plot background
                        axis.ticks = element_blank(),        # Remove axis ticks
                        axis.text = element_blank(),
                        axis.title.y = element_blank(),
                    ) +
                    labs(y = label, x = s, color='Density') -> p1
            }
        }  else{
            if (n == 1){
                bg_embs %>% ggplot( aes(x = UMAP_1, y = UMAP_2, color=density)) + 
                    geom_point(aes(x = UMAP_1, y = UMAP_2), color='lightgray', data = bg_embs, size = .5) +
                    theme_minimal() +
                    theme(
                        panel.background = element_blank(),  # Remove panel background
                        panel.grid.major = element_blank(),  # Remove major gridlines
                        panel.grid.minor = element_blank(),  # Remove minor gridlines
                        plot.background = element_blank(),   # Remove plot background
                        axis.ticks = element_blank(),        # Remove axis ticks
                        axis.text = element_blank(),
                        legend.position = "none"
                    ) +
                    labs(y = label, x = s, fill=',') -> p1
            } else if (n == 2){
                bg_embs %>% ggplot( aes(x = UMAP_1, y = UMAP_2, color=density)) + 
                    geom_point(aes(x = UMAP_1, y = UMAP_2), color='lightgray', data = bg_embs, size = .5) +
                    theme_minimal() +
                    theme(
                        panel.background = element_blank(),  # Remove panel background
                        panel.grid.major = element_blank(),  # Remove major gridlines
                        panel.grid.minor = element_blank(),  # Remove minor gridlines
                        plot.background = element_blank(),   # Remove plot background
                        axis.ticks = element_blank(),        # Remove axis ticks
                        axis.text = element_blank(),
                        axis.title.y = element_blank(),
                        legend.position = "none"
                    ) +
                    labs(y = label, x = s, fill=',') -> p1
            } else if (n == 3){
                bg_embs %>% ggplot( aes(x = UMAP_1, y = UMAP_2, color=density)) + 
                    geom_point(aes(x = UMAP_1, y = UMAP_2), color='lightgray', data = bg_embs, size = .5) +
                    theme_minimal() +
                    theme(
                        panel.background = element_blank(),  # Remove panel background
                        panel.grid.major = element_blank(),  # Remove major gridlines
                        panel.grid.minor = element_blank(),  # Remove minor gridlines
                        plot.background = element_blank(),   # Remove plot background
                        axis.ticks = element_blank(),        # Remove axis ticks
                        axis.text = element_blank(),
                        axis.title.y = element_blank(),
                    ) +
                    labs(y = label, x = s, color='Density') -> p1
            }
        }
        plot.list[[n]] <- p1
        n <- n + 1
    }

    options(repr.plot.width = 10, repr.plot.height = 3)
    tiff(paste0('v2_figures/6_',sprintf("%02d", number),'_',gene,'_density.tiff'), width = 10, height = 3, res = 300, units = 'in')
    print(wrap_plots(plot.list, ncol = 3))
    dev.off()

    number <- number + 1



}

********************************************************************
********************************V3_plot*****************************
********************************************************************

********************************************
******************lineagegene***************
********************************************

library(Seurat)
library(ggplot2)
library(tidyr)

hTLO_v3 <- qs::qread(file = 'output/02_output/01_NTSM_integrated_decontX.qs')
##set colors
writeLines("gene_colors <- c('#E41A1C', '#377EB8', '#4DAF4A', '#984EA3','#FF7F00', '#FFFF33', '#A65628', '#F781BF', 
'#66C2A5', '#FC8D62', '#8DA0CB', '#E78AC3','#A6D854', '#FFD92F')",
           "resource/gene_colors.R")
source('resource/gene_colors.R')
##NT_AP
NT_AP=data.table::fread(file = 'output/02_output/06_NT_AP_palantir_result.csv',data.table = F)
colnames(NT_AP)[1] <- 'cell.id'
rownames(NT_AP) <- NT_AP[,1]
NT_GOI=FetchData(hTLO_v3,vars = c("SOX2", "TBXT", "CYP26A1", "CDX1", "WNT5B", "NKX1-2", 
                                  "HES1",  "NEFM", "PAX6", "SOX1", 
                                  "NKX6-1"))

NT_GOI$NMP <- ifelse(NT_GOI$TBXT >0 & NT_GOI$SOX2 >0,'TBXT+SOX2+','other')

GOI = cbind(rownames(NT_GOI),NT_GOI)
colnames(GOI)[1] <- 'cell.id'

exprSet <- merge(NT_AP,GOI, by = 'cell.id')
rownames(exprSet) <- exprSet[,1]
exprSet2  <- exprSet[,-c(1:10,12)]

write.csv(exprSet2,file = "output/04_task1/03_slide3/03_NT_AP.csv",row.names = F)
exprSet2<- data.table::fread(file = 'output/04_task1/03_slide3/03_NT_AP.csv',data.table = F)

data <- exprSet2 %>% 
  pivot_longer(col=2:12,
               names_to = 'gene',
               values_to = 'expression')

data$gene <- as.factor(data$gene)
data$gene <- factor(data$gene, levels = c("SOX2", "TBXT", "CYP26A1", "CDX1", "WNT5B", "NKX1-2", 
                                          "HES1", "NEFM", "PAX6", "SOX1","NKX6-1"))

##做图
library(grid)
library(ggplot2)
highlight_data <- data[data$NMP =='TBXT+SOX2+',]

ggplot(data = data, aes(x = pseudotime, y = expression, color = gene)) +
  geom_smooth(se = TRUE, alpha = 0.3,size=1.5) + 
  scale_color_manual(values = gene_colors) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(), 
    axis.line.x = element_blank(),
    axis.line.y = element_line(size = 0.8, color = "black"), 
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(size = 0.8, color = "black"), 
    axis.text = element_text(size = 12, face = "bold"), 
    axis.title = element_text(size = 14, face = "bold"), 
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold") 
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25), 
    labels = scales::number_format(accuracy = 0.01)
  ) +
  scale_y_continuous(
    expand = c(0, 0), 
    breaks = seq(0, max(data$expression, na.rm = TRUE), by = 0.5) 
  ) +
  geom_segment(
    data = highlight_data,  # 使用高亮数据
    aes(x = pseudotime, xend = pseudotime, y = -0.05, yend = -0.1),  # 每个细胞的条形线
    inherit.aes = FALSE,
    color = "black", size = 0.1  # 条形线颜色和宽度
  ) +
  labs(
    title = "Gene Expression Over Pseudotime",
    x = NULL, 
    y = "Normalized Gene Expression"
  )
ggsave(filename = 'output/04_task1/03_slide3/NT_GOI_optimized.tiff',height = 5, width = 10,dpi = 300,bg = 'white')




##SM_AP
SM_AP=data.table::fread(file = 'output/02_output/06_SM_AP_palantir_result.csv',data.table = F)
colnames(SM_AP)[1] <- 'cell.id'
rownames(SM_AP) <- SM_AP[,1]
SM_GOI=FetchData(hTLO_v3,vars = c('SOX2', 'TBXT', 'CYP26A1', 'CDX2', 'WNT8A', 'HES7', 'LFNG', 'MSGN1', 
                                  'MEOX1', 'TCF15', 'PAX3', 'PAX1', 'TWIST1', 'CTNNB1'))
SM_GOI$NMP <- ifelse(SM_GOI$TBXT >0 & SM_GOI$SOX2 >0,'TBXT+SOX2+','other')
GOI = cbind(rownames(SM_GOI),SM_GOI)
colnames(GOI)[1] <- 'cell.id'

exprSet <- merge(SM_AP,GOI, by = 'cell.id')
rownames(exprSet) <- exprSet[,1]
exprSet2  <- exprSet[,-c(1:10,12)]

write.csv(exprSet2,file = "output/04_task1/03_slide3/03_SM_AP.csv",row.names = F)

data <- exprSet2 %>% 
  pivot_longer(col=2:15,
               names_to = 'gene',
               values_to = 'expression')
data$gene <- as.factor(data$gene)
data$gene <- factor(data$gene, levels = c('SOX2', 'TBXT', 'CYP26A1', 'CDX2', 'WNT8A', 'HES7', 'LFNG', 'MSGN1', 
                                          'MEOX1', 'TCF15', 'PAX3', 'PAX1', 'TWIST1', 'CTNNB1'))

library(grid)
library(ggplot2)
highlight_data <- data[data$NMP =='TBXT+SOX2+',]

ggplot(data = data, aes(x = pseudotime, y = expression, color = gene)) +
  geom_smooth(se = TRUE, alpha = 0.3,size=1.5) +
  scale_color_manual(values = gene_colors) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    panel.grid = element_blank(),
    axis.line.x = element_blank(), 
    axis.line.y = element_line(size = 0.8, color = "black"), 
    axis.ticks.x = element_blank(), 
    axis.ticks.y = element_line(size = 0.8, color = "black"), 
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold") 
  ) +
  scale_x_continuous(
    limits = c(0, 1), 
    breaks = seq(0, 1, by = 0.25), 
    labels = scales::number_format(accuracy = 0.01)
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    breaks = seq(0, max(data$expression, na.rm = TRUE), by = 0.5) 
  ) +
  geom_segment(
    data = highlight_data,
    aes(x = pseudotime, xend = pseudotime, y = -0.05, yend = -0.1), 
    inherit.aes = FALSE,
    color = "black", size = 0.1 
  ) +
  labs(
    title = "Gene Expression Over Pseudotime",
    x = NULL, 
    y = "Normalized Gene Expression"
  )

ggsave(filename = 'output/04_task1/03_slide3/SM_GOI_optimized.tiff',height = 5, width = 10,dpi = 300,bg = 'white')


##Noto
Noto <- data.table::fread(file = 'output/02_output/06_Noto_palantir_result.csv',data.table = F)
Noto <- Noto[,c(1:11)]
colnames(Noto)[1] <- 'cell.id'
rownames(Noto) <- Noto[,1]
Noto_GOI=FetchData(hTLO_v3,vars = c(
  'SOX2','TBXT','NOTO', 'NODAL','SHH','TPPP3', 'CYP26A1','FOXA2','FOXJ1','WNT3A', 'C1orf189','CHRD'))
Noto_GOI$NMP <- ifelse(Noto_GOI$TBXT >0 & Noto_GOI$SOX2 >0,'TBXT+SOX2+','other')
GOI = cbind(rownames(Noto_GOI),Noto_GOI)
colnames(GOI)[1] <- 'cell.id'

exprSet <- merge(Noto,GOI, by = 'cell.id')
rownames(exprSet) <- exprSet[,1]
exprSet2  <- exprSet[,-c(1:9,11)]

write.csv(exprSet2,file = "output/04_task1/03_slide3/03_Noto.csv",row.names = F)
#main figfure
exprSet3 <- exprSet2[,-c(2,3,9:13)]
#sup figgfure
exprSet4 <- exprSet2[,-c(4:8)]


#main figure
data3 <- exprSet3 %>% 
  pivot_longer(col=2:6,
               names_to = 'gene',
               values_to = 'expression')
data3 <- as.data.frame(data3)
data3$gene <- as.factor(data3$gene)
data3$gene <- factor(data3$gene, levels = c('NODAL', 'NOTO', 'SHH', 'TPPP3','CYP26A1'))

##supp
data4 <- exprSet4 %>% 
  pivot_longer(col=2:8,
               names_to = 'gene',
               values_to = 'expression')
data4 <- as.data.frame(data4)
data4$gene <- as.factor(data4$gene)
data4$gene <- factor(data4$gene, levels = c('SOX2','TBXT','FOXA2','FOXJ1','WNT3A','C1orf189','CHRD'))

library(grid)
library(ggplot2)
highlight_data3 <- data3[data3$NMP =='TBXT+SOX2+',]

ggplot(data = data3, aes(x = pseudotime, y = expression, color = gene)) +
  geom_smooth(se = TRUE, alpha = 0.3,size=1.5) + 
  scale_color_manual(values = gene_colors) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA), 
    panel.grid = element_blank(), 
    axis.line.x = element_blank(), 
    axis.line.y = element_line(size = 0.8, color = "black"), 
    axis.ticks.x = element_blank(), 
    axis.ticks.y = element_line(size = 0.8, color = "black"), 
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"), 
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold") 
  ) +
  scale_x_continuous(
    limits = c(0, 1), 
    breaks = seq(0, 1, by = 0.25),
    labels = scales::number_format(accuracy = 0.01) 
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    breaks = seq(0, max(data$expression, na.rm = TRUE), by = 0.5)
  ) +
  geom_segment(
    data = highlight_data3, 
    aes(x = pseudotime, xend = pseudotime, y = -0.05, yend = -0.1),  
    inherit.aes = FALSE,
    color = "black", size = 0.1 
  ) +
  labs(
    title = "Gene Expression Over Pseudotime",
    x = NULL, 
    y = "Normalized Gene Expression"
  )

ggsave(filename = 'output/04_task1/03_slide3/NOTO_GOI_supplefigure_optimized.tiff',height = 5, width = 10,dpi = 300,bg = 'white')


********************************************
**************cellcount_piechart************
********************************************

library(Seurat)
library(ggplot2)
library(dplyr)

Noto <- qs::qread(file = 'output/02_output/01_Noto.qs')

markers <- c("NOTO","TBXT","FOXA2","FOXA1","SHH",'CHRD','LEFTY1','BAMBI','NOG','LEFTY2')

Expression <- FetchData(Noto, vars = markers)

NOTO <- sum(Expression$NOTO>0)
TN <- sum(Expression$NOTO>0 & Expression$TBXT >0)
TF2 <- sum(Expression$TBXT>0 & Expression$FOXA2 >0)
TF1 <- sum(Expression$TBXT>0 & Expression$FOXA1 >0)
TS <- sum(Expression$TBXT>0 & Expression$SHH>0)
NS <- sum(Expression$NOTO>0 & Expression$SHH>0)
NF2 <- sum(Expression$NOTO>0 & Expression$FOXA2>0)
NF1 <- sum(Expression$NOTO>0 & Expression$FOXA1>0)
F2S <- sum(Expression$FOXA2>0 & Expression$SHH>0)
F2F1 <- sum(Expression$FOXA2>0 & Expression$FOXA1>0)
F1S <- sum(Expression$FOXA1>0 & Expression$SHH>0)
CHRD <- sum(Expression$CHRD>0)
CL1 <- sum(Expression$CHRD>0 & Expression$LEFTY1>0)
CB <- sum(Expression$CHRD>0 & Expression$BAMBI>0)
CN <- sum(Expression$CHRD>0 & Expression$NOG>0)
CL1B <- sum(Expression$CHRD>0 & Expression$LEFTY1>0 & Expression$BAMBI >0)
CL2B <- sum(Expression$CHRD>0 & Expression$LEFTY2>0 & Expression$BAMBI >0)
CL1N <- sum(Expression$CHRD>0 & Expression$LEFTY1>0 & Expression$NOG >0)
CL2N <- sum(Expression$CHRD>0 & Expression$LEFTY2>0 & Expression$NOG >0)
CBN <- sum(Expression$CHRD>0 & Expression$BAMBI>0 & Expression$NOG >0)


categories <- c('CHRD','CHRD+LEFTY1+','CHRD+BAMBI+','CHRD+NOG+','CHRD+LEFTY+BAMBI+','CHRD+LEFTY+NOG+',
                'CHRD+BAMBI+NOG+','NOTO+','TBXT+NOTO+','TBXT+FOXA2+','TBXT+FOXA1+','TBXT+SHH+',
                'NOTO+SHH+','NOTO+FOXA2+','NOTO+FOXA1+','FOXA2+SHH+','FOXA2+FOXA1+',
                'FOXA1+SHH+')

counts <- c(474,52, 98, 64, 11, 11, 19,
            836, 511, 542, 110, 419, 368, 523, 63, 406, 86, 107)

data <- data.frame(
  Category = categories,
  Count = counts
)

data$percent <- (data$Count/1599) *100
write.csv(data,file = 'output/09_task5/01_NOTO_cellcount.csv',row.names = F)


###all cell
hTLOV3 <- qs::qread(file = 'output/02_output/01_NTSM_integrated_decontX.qs')
markers <- c("NOTO","TBXT","FOXA2","FOXA1","SHH")
Expression <- FetchData(hTLOV3, vars = markers)

NOTO <- sum(Expression$NOTO>0)
TN <- sum(Expression$NOTO>0 & Expression$TBXT >0)
TF2 <- sum(Expression$TBXT>0 & Expression$FOXA2 >0)
TF1 <- sum(Expression$TBXT>0 & Expression$FOXA1 >0)
TS <- sum(Expression$TBXT>0 & Expression$SHH>0)
NS <- sum(Expression$NOTO>0 & Expression$SHH>0)
NF2 <- sum(Expression$NOTO>0 & Expression$FOXA2>0)
NF1 <- sum(Expression$NOTO>0 & Expression$FOXA1>0)
F2S <- sum(Expression$FOXA2>0 & Expression$SHH>0)
F2F1 <- sum(Expression$FOXA2>0 & Expression$FOXA1>0)
F1S <- sum(Expression$FOXA1>0 & Expression$SHH>0)

categories <- c("F1S", "F2F1", "F2S", "NF1", "NF2", "NOTO", "NS", "TF1", "TF2", "TN", "TS")
counts <- c(196, 352, 594, 66, 539, 872, 378, 115, 591, 519, 435)

data <- data.frame(
  Category = categories,
  Count = counts
)

##day3
D3 <- subset(hTLOV3, subset = orig.ident == "hTEM_v3_day3")
Expression <- FetchData(D3, vars = markers)

NOTO <- sum(Expression$NOTO>0)
TN <- sum(Expression$NOTO>0 & Expression$TBXT >0)
TF2 <- sum(Expression$TBXT>0 & Expression$FOXA2 >0)
TF1 <- sum(Expression$TBXT>0 & Expression$FOXA1 >0)
TS <- sum(Expression$TBXT>0 & Expression$SHH>0)
NS <- sum(Expression$NOTO>0 & Expression$SHH>0)
NF2 <- sum(Expression$NOTO>0 & Expression$FOXA2>0)
NF1 <- sum(Expression$NOTO>0 & Expression$FOXA1>0)
F2S <- sum(Expression$FOXA2>0 & Expression$SHH>0)
F2F1 <- sum(Expression$FOXA2>0 & Expression$FOXA1>0)
F1S <- sum(Expression$FOXA1>0 & Expression$SHH>0)

data$D3 <- c('18','30','189','29','297','501','208','30','272','283','185')

#day6
D6 <- subset(hTLOV3, subset = orig.ident == "hTEM_v3_day6")
Expression <- FetchData(D6, vars = markers)

NOTO <- sum(Expression$NOTO>0)
TN <- sum(Expression$NOTO>0 & Expression$TBXT >0)
TF2 <- sum(Expression$TBXT>0 & Expression$FOXA2 >0)
TF1 <- sum(Expression$TBXT>0 & Expression$FOXA1 >0)
TS <- sum(Expression$TBXT>0 & Expression$SHH>0)
NS <- sum(Expression$NOTO>0 & Expression$SHH>0)
NF2 <- sum(Expression$NOTO>0 & Expression$FOXA2>0)
NF1 <- sum(Expression$NOTO>0 & Expression$FOXA1>0)
F2S <- sum(Expression$FOXA2>0 & Expression$SHH>0)
F2F1 <- sum(Expression$FOXA2>0 & Expression$FOXA1>0)
F1S <- sum(Expression$FOXA1>0 & Expression$SHH>0)

data$D6 <- c('23','26','42','8','15','34','21','16','35','22','31')

data$Category <- factor(data$Category)
data$Category <- factor(data$Category, levels = c('NOTO','TN','TF2','TF1','TS','NS',
                                                  'NF2','NF1','F2S','F2F1','F1S'))
data$Category <- factor(data$Category,labels = c('NOTO+','TBXT+NOTO+','TBXT+FOXA2+','TBXT+FOXA1+','TBXT+SHH+',
                                                 'NOTO+SHH+','NOTO+FOXA2+','NOTO+FOXA1+','FOXA2+SHH+',
                                                 'FOXA2+FOXA1+','FOXA1+SHH+'))
colnames(data) <- c('celltype','all','D3','D6')

write.csv(data,file = 'output/09_task5/01_Allcells_cellcount.csv',row.names = F)



####cell number for each clusters
df <- data.frame(table(Noto$seurat_clusters))
colnames(df)<-c('cluster','cell_number')
write.csv(df,file = 'output/09_task5/01_cluster_cellnumbers',row.names = F)

##cell number for each day
Noto$group <- Noto$orig.ident
Noto$group <- factor(Noto$group,labels = c('hTLO_v3_D3', 'hTLO_V3_D4','hTLO_V3_D5',"hTLO_V3_D6",'hTLO_V3_D7'))
df <-data.frame(table(Noto$group))
colnames(df) <- c('Time','cell_number')
write.csv(df,file = 'output/09_task5/01_time_cellnumbers',row.names = F)



###subset of TBXT+ cells
TBXT <- WhichCells(Noto,expression = TBXT >0)
TBXT <- subset(Noto,subset = TBXT>0)
markers <- c('WNT3A','NODAL','CHRD','SHH','NOTO','FOXA2','FOXA1','NOG','RFX2','FOXJ1')

Expression <- FetchData(TBXT, vars = markers)
WN <-sum(Expression$WNT3A>0 & Expression$NODAL >0)
CS <- sum(Expression$CHRD>0 & Expression$SHH >0)
NF <- sum(Expression$NOTO>0 & Expression$FOXA2>0)
FN <- sum(Expression$FOXA1>0 & Expression$NOG>0)
C<- sum(Expression$CHRD>0)
S <- sum(Expression$SHH>0)
N<- sum(Expression$NOTO>0)
F<-sum(Expression$FOXA2>0)
RF<-sum(Expression$RFX2 >0 & Expression$FOXJ1>0)
R<- sum(Expression$RFX2>0)
FOXJ<-sum(Expression$FOXJ1>0)

data <- data.frame(
  Category = c('TBXT+','WNT3A+NODAL+',"CHRD+SHH+",'NOTO+FOXA2+','FOXA1+NOG+','RFX2+FOXJ1+','others','CHRD','SHH','NOTO','FOXA2','RFX2','FOXJ1' ),
  Count = c(ncol(TBXT),WN,CS,NF,FN,RF,ncol(TBXT)-WN-CS-NF-FN-RF,C,S,N,F,R,FOXJ),
  percent=c(ncol(TBXT)/ncol(Noto),WN/ncol(TBXT),CS/ncol(TBXT),NF/ncol(TBXT),FN/ncol(TBXT),RF/ncol(TBXT),(ncol(TBXT)-WN-CS-NF-FN-RF)/ncol(TBXT),
            C/ncol(TBXT),S/ncol(TBXT),N/ncol(TBXT),F/ncol(TBXT),R/ncol(TBXT),FOXJ/ncol(TBXT))*100)

df <- data[c(1:7),]
inner_data <- df[2:7,]
inner_data$percent_int<- paste0(round(inner_data$percent), "%")
library(RColorBrewer)
library(plotrix)
inner_data$Category<-factor(inner_data$Category,levels  = c('NOTO+FOXA2+','CHRD+SHH+','RFX2+FOXJ1+','WNT3A+NODAL+','FOXA1+NOG+','others'))
inner_data<- inner_data[order(inner_data$Category),]

pdf("output/09_task5/01_NOTO_cellcount_allcluster.pdf", width = 10, height = 6) 
pie(inner_data$percent , labels = paste0(inner_data$Category," : ",inner_data$percent_int),border='white',col=brewer.pal(6, "Set2"),main = "TBXT+ Cells")
draw.circle(0, 0, radius = 0.8, border = "red", lwd = 16)
dev.off()



##cluster 1
Noto_1 <- subset(Noto, subset = celltype2 =='Noto.1')
TBXT_1 <- subset(Noto_1,subset = TBXT>0)
markers <- c('WNT3A','NODAL','CHRD','SHH','NOTO','FOXA2','FOXA1','NOG','RFX2','FOXJ1')

Expression <- FetchData(TBXT_1, vars = markers)
WN <-sum(Expression$WNT3A>0 & Expression$NODAL >0)
CS <- sum(Expression$CHRD>0 & Expression$SHH >0)
NF <- sum(Expression$NOTO>0 & Expression$FOXA2>0)
FN <- sum(Expression$FOXA1>0 & Expression$NOG>0)
C<- sum(Expression$CHRD>0)
S <- sum(Expression$SHH>0)
N<- sum(Expression$NOTO>0)
F<-sum(Expression$FOXA2>0)
RF<-sum(Expression$RFX2 >0 & Expression$FOXJ1>0)
R<- sum(Expression$RFX2>0)
FOXJ<-sum(Expression$FOXJ1>0)
# 创建 DataFrame
data <- data.frame(
  Category = c('TBXT+','WNT3A+NODAL+',"CHRD+SHH+",'NOTO+FOXA2+','FOXA1+NOG+','RFX2+FOXJ1+','others','CHRD','SHH','NOTO','FOXA2','RFX2','FOXJ1' ),
  Count = c(ncol(TBXT_1),WN,CS,NF,FN,RF,ncol(TBXT_1)-WN-CS-NF-FN-RF,C,S,N,F,R,FOXJ),
  percent=c(ncol(TBXT_1)/ncol(Noto_1),WN/ncol(TBXT_1),CS/ncol(TBXT_1),NF/ncol(TBXT_1),FN/ncol(TBXT_1),RF/ncol(TBXT_1),(ncol(TBXT_1)-WN-CS-NF-FN-RF)/ncol(TBXT_1),
            C/ncol(TBXT_1),S/ncol(TBXT_1),N/ncol(TBXT_1),F/ncol(TBXT_1),R/ncol(TBXT_1),FOXJ/ncol(TBXT_1))*100)

df <- data[c(1:7),]
inner_data <- df[2:7,]
inner_data$percent_int<- paste0(round(inner_data$percent), "%")
library(RColorBrewer)
library(plotrix)
inner_data$Category<-factor(inner_data$Category,levels  = c('NOTO+FOXA2+','CHRD+SHH+','RFX2+FOXJ1+','WNT3A+NODAL+','FOXA1+NOG+','others'))
inner_data<- inner_data[order(inner_data$Category),]

pdf("output/09_task5/01_NOTO_cellcount_cluster_1.pdf", width = 10, height = 6) 
pie(inner_data$percent , labels = paste0(inner_data$Category," : ",inner_data$percent_int),border='white',col=brewer.pal(6, "Set2"),main = "TBXT+ Cells")
draw.circle(0, 0, radius = 0.8, border = "red", lwd = 16)
dev.off()


##cluster 2
Noto_2 <- subset(Noto, subset = celltype2 =='Noto.2')
TBXT_2 <- subset(Noto_2,subset = TBXT>0)
markers <- c('WNT3A','NODAL','CHRD','SHH','NOTO','FOXA2','FOXA1','NOG','RFX2','FOXJ1')

Expression <- FetchData(TBXT_2, vars = markers)
WN <-sum(Expression$WNT3A>0 & Expression$NODAL >0)
CS <- sum(Expression$CHRD>0 & Expression$SHH >0)
NF <- sum(Expression$NOTO>0 & Expression$FOXA2>0)
FN <- sum(Expression$FOXA1>0 & Expression$NOG>0)
C<- sum(Expression$CHRD>0)
S <- sum(Expression$SHH>0)
N<- sum(Expression$NOTO>0)
F<-sum(Expression$FOXA2>0)
RF<-sum(Expression$RFX2 >0 & Expression$FOXJ1>0)
R<- sum(Expression$RFX2>0)
FOXJ<-sum(Expression$FOXJ1>0)
data <- data.frame(
  Category = c('TBXT+','WNT3A+NODAL+',"CHRD+SHH+",'NOTO+FOXA2+','FOXA1+NOG+','RFX2+FOXJ1+','others','CHRD','SHH','NOTO','FOXA2','RFX2','FOXJ1' ),
  Count = c(ncol(TBXT_2),WN,CS,NF,FN,RF,ncol(TBXT_2)-WN-CS-NF-FN-RF,C,S,N,F,R,FOXJ),
  percent=c(ncol(TBXT_2)/ncol(Noto_2),WN/ncol(TBXT_2),CS/ncol(TBXT_2),NF/ncol(TBXT_2),FN/ncol(TBXT_2),RF/ncol(TBXT_2),(ncol(TBXT_2)-WN-CS-NF-FN-RF)/ncol(TBXT_2),
            C/ncol(TBXT_2),S/ncol(TBXT_2),N/ncol(TBXT_2),F/ncol(TBXT_2),R/ncol(TBXT_2),FOXJ/ncol(TBXT_2))*100)

df <- data[c(1:7),]
inner_data <- df[2:7,]
inner_data$percent_int<- paste0(round(inner_data$percent), "%")
library(RColorBrewer)
library(plotrix)
inner_data$Category<-factor(inner_data$Category,levels  = c('NOTO+FOXA2+','CHRD+SHH+','RFX2+FOXJ1+','WNT3A+NODAL+','FOXA1+NOG+','others'))
inner_data<- inner_data[order(inner_data$Category),]

pdf("output/09_task5/01_NOTO_cellcount_cluster_2.pdf", width = 10, height = 6) 
pie(inner_data$percent , labels = paste0(inner_data$Category," : ",inner_data$percent_int),border='white',col=brewer.pal(6, "Set2"),main = "TBXT+ Cells")
draw.circle(0, 0, radius = 0.8, border = "red", lwd = 16)
dev.off()



##cluster 3
Noto_3 <- subset(Noto, subset = celltype2 =='Noto.3')
TBXT_3 <- subset(Noto_3,subset = TBXT>0)
markers <- c('WNT3A','NODAL','CHRD','SHH','NOTO','FOXA2','FOXA1','NOG','RFX2','FOXJ1')

Expression <- FetchData(TBXT_3, vars = markers)
WN <-sum(Expression$WNT3A>0 & Expression$NODAL >0)
CS <- sum(Expression$CHRD>0 & Expression$SHH >0)
NF <- sum(Expression$NOTO>0 & Expression$FOXA2>0)
FN <- sum(Expression$FOXA1>0 & Expression$NOG>0)
C<- sum(Expression$CHRD>0)
S <- sum(Expression$SHH>0)
N<- sum(Expression$NOTO>0)
F<-sum(Expression$FOXA2>0)
RF<-sum(Expression$RFX2 >0 & Expression$FOXJ1>0)
R<- sum(Expression$RFX2>0)
FOXJ<-sum(Expression$FOXJ1>0)

data <- data.frame(
  Category = c('TBXT+','WNT3A+NODAL+',"CHRD+SHH+",'NOTO+FOXA2+','FOXA1+NOG+','RFX2+FOXJ1+','others','CHRD','SHH','NOTO','FOXA2','RFX2','FOXJ1' ),
  Count = c(ncol(TBXT_3),WN,CS,NF,FN,RF,ncol(TBXT_3)-WN-CS-NF-FN-RF,C,S,N,F,R,FOXJ),
  percent=c(ncol(TBXT_3)/ncol(Noto_3),WN/ncol(TBXT_3),CS/ncol(TBXT_3),NF/ncol(TBXT_3),FN/ncol(TBXT_3),RF/ncol(TBXT_3),(ncol(TBXT_3)-WN-CS-NF-FN-RF)/ncol(TBXT_3),
            C/ncol(TBXT_3),S/ncol(TBXT_3),N/ncol(TBXT_3),F/ncol(TBXT_3),R/ncol(TBXT_3),FOXJ/ncol(TBXT_3))*100)

df <- data[c(1:6),]
inner_data <- df[2:6,]
inner_data$percent_int<- paste0(round(inner_data$percent), "%")
library(RColorBrewer)
library(plotrix)
inner_data$Category<-factor(inner_data$Category,levels  = c('NOTO+FOXA2+','CHRD+SHH+','RFX2+FOXJ1+','WNT3A+NODAL+','FOXA1+NOG+','others'))
inner_data<- inner_data[order(inner_data$Category),]

pdf("output/09_task5/01_NOTO_cellcount_cluster_3.pdf", width = 10, height = 6) 
pie(inner_data$percent , labels = paste0(inner_data$Category," : ",inner_data$percent_int),border='white',col=brewer.pal(5, "Set2"),main = "TBXT+ Cells")
draw.circle(0, 0, radius = 0.8, border = "red", lwd = 16)
dev.off()


##cluster 4
Noto_4<- subset(Noto, subset = celltype2 =='Noto.4')
TBXT_4 <- subset(Noto_4,subset = TBXT>0)
markers <- c('WNT3A','NODAL','CHRD','SHH','NOTO','FOXA2','FOXA1','NOG','RFX2','FOXJ1')

Expression <- FetchData(TBXT_4, vars = markers)
WN <-sum(Expression$WNT3A>0 & Expression$NODAL >0)
CS <- sum(Expression$CHRD>0 & Expression$SHH >0)
NF <- sum(Expression$NOTO>0 & Expression$FOXA2>0)
FN <- sum(Expression$FOXA1>0 & Expression$NOG>0)
C<- sum(Expression$CHRD>0)
S <- sum(Expression$SHH>0)
N<- sum(Expression$NOTO>0)
F<-sum(Expression$FOXA2>0)
RF<-sum(Expression$RFX2 >0 & Expression$FOXJ1>0)
R<- sum(Expression$RFX2>0)
FOXJ<-sum(Expression$FOXJ1>0)
# 创建 DataFrame
data <- data.frame(
  Category = c('TBXT+','WNT3A+NODAL+',"CHRD+SHH+",'NOTO+FOXA2+','FOXA1+NOG+','RFX2+FOXJ1+','others','CHRD','SHH','NOTO','FOXA2','RFX2','FOXJ1' ),
  Count = c(ncol(TBXT_4),WN,CS,NF,FN,RF,ncol(TBXT_4)-WN-CS-NF-FN-RF,C,S,N,F,R,FOXJ),
  percent=c(ncol(TBXT_4)/ncol(Noto_4),WN/ncol(TBXT_4),CS/ncol(TBXT_4),NF/ncol(TBXT_4),FN/ncol(TBXT_4),RF/ncol(TBXT_4),(ncol(TBXT_4)-WN-CS-NF-FN-RF)/ncol(TBXT_4),
            C/ncol(TBXT_4),S/ncol(TBXT_4),N/ncol(TBXT_4),F/ncol(TBXT_4),R/ncol(TBXT_4),FOXJ/ncol(TBXT_4))*100)

df <- data[c(1:7),]
inner_data <- df[2:7,]
inner_data$percent_int<- paste0(round(inner_data$percent), "%")
library(RColorBrewer)
library(plotrix)
inner_data$Category<-factor(inner_data$Category,levels  = c('NOTO+FOXA2+','CHRD+SHH+','RFX2+FOXJ1+','WNT3A+NODAL+','FOXA1+NOG+','others'))
inner_data<- inner_data[order(inner_data$Category),]

pdf("output/09_task5/01_NOTO_cellcount_cluster_4.pdf", width = 10, height = 6) 
pie(inner_data$percent , labels = paste0(inner_data$Category," : ",inner_data$percent_int),border='white',col=brewer.pal(6, "Set2"),main = "TBXT+ Cells")
draw.circle(0, 0, radius = 0.8, border = "red", lwd = 16)
dev.off()


combined_data <- bind_cols(data,data_1,data_2,data_3,data_4)

colnames(combined_data) <- c('cluster_all','category_all','percent_all',
                             'cluster_1','category_1','percent_1',
                             'cluster_2','category_2','percent_2',
                             'cluster_3','category_3','percent_3',
                             'cluster_4','category_4','percent_4')
write.csv(combined_data, file = "output/09_task5/01_cluster_cells.csv", row.names = FALSE)


####pie chart

data <- read.csv(file = 'output/09_task5/01_cluster_cells.csv')
TBXT <- WhichCells(Noto,expression = TBXT >0)

Noto_1<-subset(Noto,subset = celltype2 == 'Noto.1')
TBXT_1 <- WhichCells(Noto_1, expression = TBXT> 0)

Noto_2<-subset(Noto,subset = celltype2 == 'Noto.2')
TBXT_2<- WhichCells(Noto_2,expression = TBXT >0)

Noto_3<-subset(Noto,subset = celltype2 == 'Noto.3')
TBXT_3<- WhichCells(Noto_3,expression = TBXT >0)

Noto_4<-subset(Noto,subset = celltype2 == 'Noto.4')
TBXT_4<- WhichCells(Noto_4,expression = TBXT >0)

df <- data.frame(
  Category <- c('All','Cluster_1','Cluster_2','Cluster_3','Cluster_4'),
  Counts <- c(length(TBXT),length(TBXT_1),length(TBXT_2),length(TBXT_3),length(TBXT_4)),
  percent <- c(length(TBXT)/ncol(Noto)*100,length(TBXT_1)/ncol(Noto_1)*100,length(TBXT_2)/ncol(Noto_2)*100,length(TBXT_3)/ncol(Noto_3)*100,length(TBXT_4)/ncol(Noto_4)*100))
colnames(df) <- c('cluster','Category','TBXT+')
df$'TBXT-'<- 100-df$`TBXT+`

write.csv(df, file = "output/09_task5/01_Noto_TBXT_cells.csv", row.names = FALSE)

data <- data.frame(
  cluster = c("All", "Cluster_1", "Cluster_2", "Cluster_3", "Cluster_4"),
  TBXT_plus = c(55.53471, 32.82675, 52.64484, 58.60349, 71.18644),
  TBXT_minus = c(44.46529, 67.17325, 47.35516, 41.39651, 28.81356)
)


data_long <- melt(data, id.vars = "cluster", 
                  variable.name = "Condition", 
                  value.name = "Percentage")

data_long <- melt(data, id.vars = "cluster", 
                  variable.name = "Condition", 
                  value.name = "Percentage")

ggplot(data_long, aes(x = cluster, y = Percentage, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +  
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), 
            position = position_dodge(0.9), vjust = -0.5, size = 3.5) +
  labs(title = "TBXT+ and TBXT- Percentages by Cluster", 
       x = "Cluster", 
       y = "Percentage") +
  scale_fill_manual(values = c("TBXT_plus" = "#1f77b4", "TBXT_minus" = "#ff7f0e"),
                    labels = c("TBXT+", "TBXT-")) +
  theme_minimal() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1), 
    panel.grid = element_blank(), 
    panel.background = element_blank(),  
    axis.line = element_line(color = "black") 
  )

********************************************
*********************dotplot****************
********************************************
library(Seurat)
library(ggplot2)
library(scales)
library(reshape2)
load('resource/tmp/colors.RData')
v3<- qs::qread(file = 'output/02_output/01_NTSM_integrated_decontX.qs')

SM <- subset(v3, detailed_celltype_2nd %in% c('E-SM','M-SM','Sclerotome','Syndetome','Endotome',
                                              'Dermomyotome','Myogenic Progenitors'))
NT <- subset(v3, detailed_celltype_2nd %in% c('Caudal Neural Plate',"Early NT",
                                              'Neural Tube','Early FP','Floor Plate'))

SM_GOI  <- rev(unique(c( "MEOX1", "LEF1", "TCF15", "IGFBP5", "COL5A2", "PDGFRA",
                         "PAX3", "SIX1", "SNAI2", "RDH10",
                         "FN1", "MYF5", "MYF6", "NEB", "KLHL41",
                         "PAX7", "MYOD1", "MYOG", "RARA", "JAG1", "NOTCH1",
                         "SNAI1", "DHRS3", "SIM1", "CDH7", "FZD2",
                         "EBF2", "EBF3", "ALCAM", "PKDCC", "BMP5", "NEFM", "CXCL12",
                         "TGFBI", "CXCR4", "FGF18", "LOX", "SCX", "FOXP1", "COL1A1",
                         "ALDH1A2", "TWIST1", "COL3A1", "PTCH1", "GLI1", "SOX9", 
                         "PAX1", "PAX9", "NKX3-2", "COL3A1", "FRZB")))            
NT_GOI <- rev(unique(c("SOX2", "SOX3", "RFX4", "NKX1-2", 'CDX1',"CDX2", "CYP26A1","FGF8",'FGF17',
                       "RARA", "RARB","WNT1", "WNT3", "WNT3A", "WNT5A", "WNT5B",'FRZB',"NRG1", "PDGFA", "IGFBP6",
                       "LMX1A", "MSX1","BMP4",'BMP7', "OLIG3", "HES1", "SFRP2", 'TTYH1', "CAMK2N1", "SOX21","PAX3", "IRX3", "IRX5", "PAX6", "GBX2", "PAX7",
                       "DBX2","NKX6-2", "SP8","HES4", "HES5", "KRT8", "EPCAM", "SAT1", 'NKX6-1',"OLIG2", "OLIG1", "NKX2-8", "NKX2-2",
                       "FOXA2", "ARX", "FOXA1", "FOXJ1","SHH", "ONECUT1", "ONECUT2", "CRABP1")))


dotplot_SM <- DotPlot(SM, features =SM_GOI, group.by = 'detailed_celltype_2nd',cols = c('lightblue','brown'),dot.scale = 10) +coord_flip()+ RotatedAxis() + ggtitle('SM')+
  theme(plot.title = element_text(hjust = 0.5))

dotplot_data_SM <- dotplot_SM$data


ggplot(dotplot_data_SM, aes(y = id, x = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_size_continuous(name = "Percent Expressed", range = c(1, 8),
                        limits = c(5,100),
                        breaks = c(5,25,50,75)) +
  scale_color_gradient(low = "lightblue", high = "brown") +
  labs(y = "Identity", x = "Features", title = "SM") +
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
ggsave(filename="output/09_task5/05_SM_dotplot.pdf",width=16,height=3.7,dpi=300,bg='white')


##############NT #########
NT$detailed_celltype_2nd = factor(NT$detailed_celltype_2nd,levels =  c('Caudal Neural Plate',"Early NT",
                                                                       'Neural Tube','Early FP','Floor Plate'))
dotplot_NT <- DotPlot(NT, features =NT_GOI, group.by = 'detailed_celltype_2nd',cols = c('lightblue','brown'),dot.scale = 10) +coord_flip()+ RotatedAxis() + ggtitle('NT')+
  theme(plot.title = element_text(hjust = 0.5))
dotplot_data_NT <- dotplot_NT$data

ggplot(dotplot_data_NT, aes(y = id, x = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_size_continuous(name = "Percent Expressed", range = c(1, 8),
                        limits = c(5,100),
                        breaks =  c(5,25,50,75)) +
  scale_color_gradient(low = "lightblue", high = "brown") +
  labs(y = "Identity", x = "Features", title = "NT") +
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
ggsave(filename="output/09_task5/05_NT_dotplot.pdf",width=16,height=3.7,dpi=300,bg='white')


********************************************
***************mountain_range***************
********************************************

library(Seurat)
library(ggplot2)
library(ggridges)
library(tidyr)

color_data <- data.frame(
  ColorCode = c("#FF1677", "#FAA31C", "#CE471C", "#D8EA60", "#F9DCAE", "#AEF700", "#D58A81", "#00FB0D", "#009285", "#F96ADB", "#0D8AFF", "#2A40FE", "#5335C1", "#D194EE", "#2E7300", "#9A8116", "#B6CAFA", "#8C8F6E", "#732E8F"),
  Name = c("Noto", "NMP-Neural", "Caud.NP", "E-NT", "E-FP", "NT", "FP", "NMP-Meso", "pPSM", "aPSM", "E-SM", "M-SM", "L-SM", "Caud.Meso", "Int.Meso", "EC", "NE", "SE", "Neuron")
)

color_values <- setNames(color_data$ColorCode,color_data$Name)
save(color_values, file = "resource/color_values.RData")
##NT_AP,
#celltype
NT_AP=data.table::fread(file = 'output/02_output/06_NT_AP_palantir_result.csv',data.table = F)
NT_AP$group = factor(NT_AP$orig.ident)
NT_AP$group =factor(NT_AP$group,labels = c('hTLO_v3_D3','hTLO_v3_D4','hTLO_v3_D5','hTLO_v3_D6','hTLO_v3_D7'))

NT_AP$celltype <- as.factor(NT_AP$celltype)
NT_AP$celltype <- factor(NT_AP$celltype, levels = c('NMP-Neural','Caud.NP','E-FP','FP','E-NT','NT'))
load('resource/color_values.RData')
ggplot(NT_AP, aes(x=pseudotime,y=celltype)) +
  geom_density_ridges(aes(fill= celltype, scale = 0.75))+
  scale_x_continuous(limits = c(0,1),expand = c(0, 0))+
  scale_fill_manual(values = color_values)+
  theme_ridges(font_size = 13, grid = FALSE)

ggsave(filename = 'output/06_task2/01_NP_celltype.pdf',height = 4, width = 8,dpi = 300,bg = "white")

#day
load('resource/day_color_values.Rdata')
ggplot(NT_AP, aes(x=pseudotime,y=group)) +
  geom_density_ridges(aes(fill= group, scale = 0.75))+
  scale_x_continuous(limits = c(0,1),expand = c(0, 0))+
  scale_fill_manual(values = day_color_values)+
  theme_ridges(font_size = 13, grid = FALSE)
ggsave(filename = 'output/06_task2/01_NP_day.tiff',height = 4, width = 8,dpi = 300,bg = "white")

####SM_AP
SM_AP=data.table::fread(file = 'output/02_output/06_SM_AP_palantir_result.csv',data.table = F)
SM_AP$group = factor(SM_AP$orig.ident)
SM_AP$group =factor(SM_AP$group,labels = c('hTLO_v3_D3','hTLO_v3_D4','hTLO_v3_D5','hTLO_v3_D6','hTLO_v3_D7'))

SM_AP$celltype <- as.factor(SM_AP$celltype)
SM_AP$celltype <- factor(SM_AP$celltype, levels = c('NMP-Meso','pPSM','aPSM','E-SM','M-SM','L-SM'))

ggplot(SM_AP, aes(x=pseudotime,y=celltype)) +
  geom_density_ridges(aes(fill= celltype, scale = 0.75))+
  scale_x_continuous(limits = c(0,1),expand = c(0, 0))+
  scale_fill_manual(values = color_values)+
  theme_ridges(font_size = 13, grid = FALSE)
ggsave(filename = 'output/06_task2/01_SM_celltype.tiff',height = 4, width = 8,dpi = 300,bg = "white")

#day
ggplot(SM_AP, aes(x=pseudotime,y=group)) +
  geom_density_ridges(aes(fill= group, scale = 0.75))+
  scale_x_continuous(limits = c(0,1),expand = c(0, 0))+
  scale_fill_manual(values = colors)+
  theme_ridges(font_size = 13, grid = FALSE)
ggsave(filename = 'output/06_task2/01_SM_day.pdf',height = 4, width = 8,dpi = 300,bg = "white")


********************************************
******************heatmap*******************
********************************************

*****01.moudle_heatmap



library(Seurat)
library(tidyverse)
library(clusterProfiler)

hTLO <- qs::qread("output/02_output/01_NTSM_integrated_decontX.qs")
Module.score <- data.table::fread(file = 'output/06_task2/v3_timeseries_with_modulescore_meta.csv',data.table = F)
rownames(Module.score)<- Module.score[,1]
Module.score <- Module.score[,-c(1:27)]
hTLO <- AddMetaData(hTLO,metadata = Module.score)
qs::qsave(hTLO,file = 'output/06_task2/03_1_hTLO.qs')
hTLO <- qs::qread(file = 'output/06_task2/03_1_hTLO.qs')
seu <- subset(hTLO,subset = celltype %in% c('Noto','L-SM','NT','NMP-Neural','NMP-Meso'))

exprSet <- FetchData(seu,vars = c('celltype',paste0("ModuleScore", 1:12)))

library(tidyr)
library(scales)
exprSet2 <- exprSet %>% 
  pivot_longer(cols = 2:13,
               names_to = 'Module',
               values_to = 'Score')
var=paste0("ModuleScore", 1:12)
exprSet2$Module <- as.factor(exprSet2$Module)
exprSet2$Module <- factor(exprSet2$Module,levels = c('ModuleScore12','ModuleScore11','ModuleScore10','ModuleScore9','ModuleScore8',"ModuleScore7",
                                                     "ModuleScore6","ModuleScore5","ModuleScore4","ModuleScore3","ModuleScore2","ModuleScore1"))
exprSet3 <- exprSet2 %>%
  group_by(celltype, Module) %>%
  mutate(
    z_score = (scale(Score) - min(scale(Score))) / (max(scale(Score)) - min(scale(Score))) * (1 - (-0.8)) +(-0.8)
  ) %>%
  ungroup() 



exprSet3$celltype <- as.factor(exprSet3$celltype)
exprSet3$celltype <- factor(exprSet3$celltype,levels = c('Noto','NMP-Meso','NMP-Neural','L-SM','NT')) 

library(dplyr)

exprSet4<- exprSet2 %>%
  group_by(celltype, Module) %>%
  mutate(
    avg_score = mean(Score)) %>% 
  ungroup() %>% 
  mutate(
    z_score = (scale(avg_score) - min(scale(avg_score))) / (max(scale(avg_score)) - min(scale(avg_score))) * 2 - 1 # 将 Z-score 缩放到 [-1, 1]
  )



library(ggplot2)
ggplot(data = exprSet3, aes(x=celltype, y=Module, fill=z_score)) + 
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "#2b83ba", high = "red",mid = 'white',
                       midpoint=0,
                       limits=c(-0.6,1))+
  theme_minimal() +  # 设置基础主题
  theme(
    panel.grid = element_blank(),      # 移除网格线
    panel.background = element_blank(), # 移除背景
    axis.ticks = element_blank(),      # 移除轴刻度
    axis.text = element_text(size = 12), # 调整轴文字大小（可选）
    axis.title = element_blank()       # 可选：移除轴标题
  )
ggsave(filename = 'output/06_task2/03_1_five_celltype_heatmap.tiff',width = 7,height = 7,dpi = 300,bg='white')
ggsave(filename = 'output/06_task2/03_1_five_celltype_heatmap.pdf',width = 7,height = 7,dpi = 300)


******0.2 signaling_heatmap



library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(Cairo)
NT <- qs::qread(file = 'output/02_output/01_NT_AP.qs')
SM <- qs::qread(file = 'output/02_output/01_SM_AP.qs')
Noto<-qs::qread(file = 'output/02_output/01_Noto.qs')
exprSet_NT <-FetchData(NT,vars = c('celltype','pseudotime'))
exprSet_SM <-FetchData(SM,vars = c('celltype','pseudotime'))
exprSet_Noto<-FetchData(Noto,vars = c('celltype2','pseudotime'))
RA <-c('ALDH1A1', 'ALDH1A2', 'ALDH1A3', 'RDH10', 'CYP26A1', 'DHRS3', 'CRABP2', 'RBP1', 'RBP2', 'RBP3', 'RBP4', 'RBP5', 'RBP7', 'RARA', 'RARB', 'RARG', 'RXRA', 'RXRB', 'RXRG', 'NCOR1', 'NCOR2', 'DHRS7B')
Notch <- c('NOTCH1', 'NOTCH2', 'NOTCH3','PRKCA', 'ST3GAL6', 'SAP30', 'DTX2', 'CCND1', 'RBX1', 'MAML2', 'LFNG', 'DLL1', 'DLL3', 'JAG1', 'DTX4', 'FZD7',
           'PPARD', 'KAT2A', 'PSENEN', 'SKP1', 'HEYL','HES1','CUL1', 'FBXW11','APH1A', 'PSEN2','HES4','HES5','HES7')
WNT <- c('WNT1', 'WNT2', 'WNT2B', 'WNT3', 'WNT3A', 'WNT4', 'WNT5A', 'WNT5B', 'WNT6', 'WNT7A', 'WNT7B', 'WNT8A', 'WNT8B', 'WNT9A', 'WNT9B', 'WNT10A', 'WNT10B', 
         'WNT11', 'WNT16', 'RSPO1', 'RSPO2', 'RSPO3', 'RSPO4', 'SFRP1', 'SFRP2', 'SFRP4', 'SFRP5', 'DKK1', 'DKK2', 'DKK3', 'DKK4', 'LEF1', 
         'CTNNB1', 'DACT1', 'DACT2', 'DACT3')
FGF <- c('FGF1', 'FGF2', 'FGF3', 'FGF4', 'FGF5', 'FGF6', 'FGF7', 'FGF8', 'FGF9', 'FGF10', 
         'FGF11', 'FGF12', 'FGF13', 'FGF14', 'FGF16', 'FGF17', 'FGF18', 'FGF19', 'FGF20', 
         'FGF21', 'FGF22', 'FGF23', 'FGFBP1', 'FGFBP2', 'FGFBP3', 'FGFR1', 'FGFR2', 
         'FGFR3', 'FGFR4', 'DUSP6', 'SPRY1', 'SPRY2', 'SPRY3', 'SPRY4', 'SHISA2', 'SHISA3', 
         'SHISA4', 'SHISA5', 'SHISA6', 'SHISA7', 'SHISA8', 'SHISA9')
HOXA <- paste0("HOXA", c(1:7,9:13))
HOXB <- paste0("HOXB", c(1:9,13))
HOXC <- paste0("HOXC", c(4:6,8:13))
HOXD <- paste0("HOXD", c(1,3,4,8:13))
##########################################################RA##################################################
exprSet_NT_RA <- FetchData(NT,vars = unique(RA))
exprSet_NT_RA <- as.data.frame(exprSet_NT_RA[,c(1:21)])
exprSet_NT_RA <- cbind(exprSet_NT,exprSet_NT_RA) %>%  .[,-1]

exprSet_SM_RA <- FetchData(SM,vars = unique(RA))
exprSet_SM_RA <- as.data.frame(exprSet_SM_RA[,c(1:21)])
exprSet_SM_RA <- cbind(exprSet_SM,exprSet_SM_RA) %>%  .[,-1]

exprSet_Noto_RA <- FetchData(Noto,vars = unique(RA))
exprSet_Noto_RA <- as.data.frame(exprSet_Noto_RA[,c(1:19)])
exprSet_Noto_RA[c("ALDH1A1", "RXRG")] <- 0
exprSet_Noto_RA <- cbind(exprSet_Noto,exprSet_Noto_RA) %>%  .[,-1]
########NT
exprSet_NT_RA1 <- exprSet_NT_RA %>% 
  pivot_longer(cols = 2:22,
               names_to = 'Gene',
               values_to = 'GeneExpression')
exprSet_NT_RA1 <- exprSet_NT_RA1 %>%
  mutate(pseudotime_binned = cut(pseudotime,breaks = seq(0, 1, length.out = 1001), labels = FALSE)) %>%
  group_by(Gene, pseudotime_binned) %>%
  summarize(GeneExpression_adjusted = mean(GeneExpression, na.rm = TRUE), .groups = "drop")

gene_matrix <- exprSet_NT_RA1 %>%
  select(Gene, GeneExpression_adjusted, pseudotime_binned) %>%
  pivot_wider(names_from = pseudotime_binned, values_from = GeneExpression_adjusted, values_fill = 0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

n_cols <- ncol(gene_matrix) 
label_positions <- seq(1, n_cols, length.out = 5)
labels <- rep("", n_cols)  # 初始化所有标记为空字符串
labels[label_positions] <- c(0, 0.25, 0.5, 0.75, 1) 

col_fun <- colorRamp2(
  breaks = c(0, 0.5, 1, 1.5,2),  # 4 个数值
  colors = c("#164788", "#E2E3E1", "#FAA31C", "#FF1677",'#D33880')  # 4 个颜色
)

tiff("output/07_task3/03_heatmapwithcluster/01_3_NT_RA_pseudotime_withcluster.tiff", width = 10, height = 7, units = "in", res = 300)
p1=Heatmap(
  gene_matrix,
  name = "Gene Expression",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "NT_RA",
  row_title = "Genes",
  row_dend_side = "left",
  column_labels = labels,# 设置横坐标为 0, 0.25, 0.5, 0.75, 1
  column_names_rot = 0,
  heatmap_legend_param = list(
    title = "Expression",
    at = c(0, 0.5, 1, 1.5, 2),
    labels = c("0", "0.5", "1", "1.5", "2")
  ),
  width = unit(10,'cm')
)
dev.off()

#######SM
exprSet_SM_RA1 <- exprSet_SM_RA %>% 
  pivot_longer(cols = 2:22,
               names_to = 'Gene',
               values_to = 'GeneExpression')
exprSet_SM_RA1 <- exprSet_SM_RA1 %>%
  mutate(pseudotime_binned = cut(pseudotime,breaks = seq(0, 1, length.out = 1001), labels = FALSE))%>%
  group_by(Gene, pseudotime_binned) %>%
  summarize(GeneExpression_adjusted = mean(GeneExpression, na.rm = TRUE), .groups = "drop")

gene_matrix <- exprSet_SM_RA1 %>%
  select(Gene, GeneExpression_adjusted, pseudotime_binned) %>%
  pivot_wider(names_from = pseudotime_binned, values_from = GeneExpression_adjusted, values_fill = 0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

n_cols <- ncol(gene_matrix) 
label_positions <- seq(1, n_cols, length.out = 5)
labels <- rep("", n_cols)  # 初始化所有标记为空字符串
labels[label_positions] <- c(0, 0.25, 0.5, 0.75, 1) 

col_fun <- colorRamp2(
  breaks = c(0, 0.5, 1, 1.5,2),  # 4 个数值
  colors = c("#164788", "#E2E3E1", "#FAA31C", "#FF1677",'#D33880')  # 4 个颜色
)

# 绘制热图
tiff("output/07_task3/03_heatmapwithcluster/01_4_SM_RA_pseudotime_withcluster.tiff", width = 10, height = 7, units = "in", res = 300)
p2=Heatmap(
  gene_matrix,
  name = "Gene Expression",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "SM_RA",
  row_title = "Genes",
  row_dend_side = "left",
  column_labels = labels,# 设置横坐标为 0, 0.25, 0.5, 0.75, 1
  column_names_rot = 0,
  heatmap_legend_param = list(
    title = "Expression",
    at = c(0, 0.5, 1, 1.5, 2),
    labels = c("0", "0.5", "1", "1.5", "2")
  ),
  width = unit(10,'cm')
)
dev.off()

#########Noto
exprSet_Noto_RA1 <- exprSet_Noto_RA %>% 
  pivot_longer(cols = 2:22,
               names_to = 'Gene',
               values_to = 'GeneExpression')
exprSet_Noto_RA1 <- exprSet_Noto_RA1 %>%
  mutate(pseudotime_binned = cut(pseudotime,breaks = seq(0, 1, length.out = 1001), labels = FALSE))%>%
  group_by(Gene, pseudotime_binned) %>%
  summarize(GeneExpression_adjusted = mean(GeneExpression, na.rm = TRUE), .groups = "drop")

gene_matrix <- exprSet_Noto_RA1 %>%
  select(Gene, GeneExpression_adjusted, pseudotime_binned) %>%
  pivot_wider(names_from = pseudotime_binned, values_from = GeneExpression_adjusted, values_fill = 0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

n_cols <- ncol(gene_matrix) 
label_positions <- seq(1, n_cols, length.out = 5)
labels <- rep("", n_cols)  # 初始化所有标记为空字符串
labels[label_positions] <- c(0, 0.25, 0.5, 0.75, 1) 

col_fun <- colorRamp2(
  breaks = c(0, 0.5, 1, 1.5,2),  # 4 个数值
  colors = c("#164788", "#E2E3E1", "#FAA31C", "#FF1677",'#D33880')  # 4 个颜色
)

# 绘制热图
tiff("output/07_task3/03_heatmapwithcluster/01_5_Noto_RA_pseudotime_withcluster.tiff", width = 10, height = 7, units = "in", res = 300)
p3=Heatmap(
  gene_matrix,
  name = "Gene Expression",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "NOTO_RA",
  row_title = "Genes",
  row_dend_side = "left",
  column_labels = labels,# 设置横坐标为 0, 0.25, 0.5, 0.75, 1
  column_names_rot = 0,
  heatmap_legend_param = list(
    title = "Expression",
    at = c(0, 0.5, 1, 1.5, 2),
    labels = c("0", "0.5", "1", "1.5", "2")
  ),
  width = unit(10,'cm')
)
dev.off()

##paired
CairoTIFF("output/07_task3/03_heatmapwithcluster/RA_paired.tiff", width = 15, height = 5, units = "in", res = 300)
p1+p2+p3
dev.off()

##########################################################WNT######################################################################
exprSet_NT_WNT <- FetchData(NT,vars = unique(WNT))
exprSet_NT_WNT <- as.data.frame(exprSet_NT_WNT)
missing_genes <- setdiff(unique(WNT),colnames(exprSet_NT_WNT))
exprSet_NT_WNT[missing_genes]<-0
exprSet_NT_WNT <- cbind(exprSet_NT,exprSet_NT_WNT) %>%  .[,-1]

exprSet_SM_WNT <- FetchData(SM,vars = unique(WNT))
exprSet_SM_WNT <- as.data.frame(exprSet_SM_WNT)
missing_genes <- setdiff(unique(WNT),colnames(exprSet_SM_WNT))
exprSet_SM_WNT[missing_genes]<-0
exprSet_SM_WNT <- cbind(exprSet_SM,exprSet_SM_WNT) %>%  .[,-1]

exprSet_Noto_WNT <- FetchData(Noto,vars = unique(WNT))
exprSet_Noto_WNT <- as.data.frame(exprSet_Noto_WNT)
missing_genes <- setdiff(unique(WNT),colnames(exprSet_Noto_WNT))
exprSet_Noto_WNT[missing_genes]<-0
exprSet_Noto_WNT <- cbind(exprSet_Noto,exprSet_Noto_WNT) %>%  .[,-1]
#########NT
exprSet_NT_WNT1 <- exprSet_NT_WNT %>% 
  pivot_longer(cols = 2:37,
               names_to = 'Gene',
               values_to = 'GeneExpression')
exprSet_NT_WNT1 <- exprSet_NT_WNT1 %>%
  mutate(pseudotime_binned = cut(pseudotime,breaks = seq(0, 1, length.out = 1001),labels = FALSE)) %>%
  group_by(Gene, pseudotime_binned) %>%
  summarize(GeneExpression_adjusted = mean(GeneExpression, na.rm = TRUE), .groups = "drop")

gene_matrix <- exprSet_NT_WNT1 %>%
  select(Gene, GeneExpression_adjusted, pseudotime_binned) %>%
  pivot_wider(names_from = pseudotime_binned, values_from = GeneExpression_adjusted, values_fill = 0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

n_cols <- ncol(gene_matrix)
label_positions <- seq(1, n_cols, length.out = 5)
labels <- rep("", n_cols)  # 初始化所有标记为空字符串
labels[label_positions] <- c(0, 0.25, 0.5, 0.75, 1)

col_fun <- colorRamp2(
  breaks = c(0, 0.5, 1, 1.5,2),  # 4 个数值
  colors = c("#164788", "#E2E3E1", "#FAA31C", "#FF1677",'#D33880')  # 4 个颜色
)

# 绘制热图
tiff("output/07_task3/03_heatmapwithcluster/01_3_NT_WNT_pseudotime_withcluster.tiff", 
     width = 10, height = 7, units = "in", res = 300)
Heatmap(
  gene_matrix,
  name = "Gene Expression",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "NT_WNT",
  row_title = "Genes",
  row_dend_side = "left",
  column_labels = labels, # 设置横坐标为 0, 0.25, 0.5, 0.75, 1
  column_names_rot = 0,
  heatmap_legend_param = list(
    title = "Expression",
    at = c(0, 0.5, 1, 1.5, 2),
    labels = c("0", "0.5", "1", "1.5", "2")
  ),
  width = unit(10,'cm')
)
dev.off()
################SM
exprSet_SM_WNT1 <- exprSet_SM_WNT %>% 
  pivot_longer(cols = 2:37,
               names_to = 'Gene',
               values_to = 'GeneExpression')
exprSet_SM_WNT1 <- exprSet_SM_WNT1 %>%
  mutate(pseudotime_binned = cut(pseudotime,breaks = seq(0, 1, length.out = 1001),abels = FALSE))%>%
  group_by(Gene, pseudotime_binned) %>%
  summarize(GeneExpression_adjusted = mean(GeneExpression, na.rm = TRUE), .groups = "drop")

gene_matrix <- exprSet_SM_WNT1 %>%
  select(Gene, GeneExpression_adjusted, pseudotime_binned) %>%
  pivot_wider(names_from = pseudotime_binned, values_from = GeneExpression_adjusted, values_fill = 0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

n_cols <- ncol(gene_matrix)
label_positions <- seq(1, n_cols, length.out = 5)
labels <- rep("", n_cols)
labels[label_positions] <- c(0, 0.25, 0.5, 0.75, 1)

col_fun <- colorRamp2(
  breaks = c(0, 0.5, 1, 1.5,2),  # 4 个数值
  colors = c("#164788", "#E2E3E1", "#FAA31C", "#FF1677",'#D33880')  # 4 个颜色
)

# 绘制热图
tiff("output/07_task3/03_heatmapwithcluster/01_4_SM_WNT_pseudotime_withcluster.tiff", 
     width = 10, height = 7,units = "in", res = 300)
Heatmap(
  gene_matrix,
  name = "Gene Expression",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "SM_WNT",
  row_title = "Genes",
  row_dend_side = "left",
  column_labels = labels, # 设置横坐标为 0, 0.25, 0.5, 0.75, 1
  column_names_rot = 0,
  heatmap_legend_param = list(
    title = "Expression",
    at = c(0, 0.5, 1, 1.5, 2),
    labels = c("0", "0.5", "1", "1.5", "2")
  ),
  width = unit(10,'cm')
)
dev.off()

#########Noto
exprSet_Noto_WNT1 <- exprSet_Noto_WNT %>% 
  pivot_longer(cols = 2:37,
               names_to = 'Gene',
               values_to = 'GeneExpression')
exprSet_Noto_WNT1 <- exprSet_Noto_WNT1 %>%
  mutate(pseudotime_binned = cut(pseudotime,breaks = seq(0, 1, length.out = 1001), labels = FALSE))%>%
  group_by(Gene, pseudotime_binned) %>%
  summarize(GeneExpression_adjusted = mean(GeneExpression, na.rm = TRUE), .groups = "drop")

gene_matrix <- exprSet_Noto_WNT1 %>%
  select(Gene, GeneExpression_adjusted, pseudotime_binned) %>%
  pivot_wider(names_from = pseudotime_binned, values_from = GeneExpression_adjusted, values_fill = 0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

n_cols <- ncol(gene_matrix) 
label_positions <- seq(1, n_cols, length.out = 5)
labels <- rep("", n_cols)  # 初始化所有标记为空字符串
labels[label_positions] <- c(0, 0.25, 0.5, 0.75, 1) 

col_fun <- colorRamp2(
  breaks = c(0, 0.5, 1, 1.5,2),  # 4 个数值
  colors = c("#164788", "#E2E3E1", "#FAA31C", "#FF1677",'#D33880')  # 4 个颜色
)

# 绘制热图
tiff("output/07_task3/03_heatmapwithcluster/01_5_Noto_WNT_pseudotime_withcluster.tiff", width = 10, height = 7, units = "in", res = 300)
p3=Heatmap(
  gene_matrix,
  name = "Gene Expression",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "NOTO_WNT",
  row_title = "Genes",
  row_dend_side = "left",
  column_labels = labels,# 设置横坐标为 0, 0.25, 0.5, 0.75, 1
  column_names_rot = 0,
  heatmap_legend_param = list(
    title = "Expression",
    at = c(0, 0.5, 1, 1.5, 2),
    labels = c("0", "0.5", "1", "1.5", "2")
  ),
  width = unit(10,'cm')
)
dev.off()

##paired
CairoTIFF("output/07_task3/03_heatmapwithcluster/WNT_paired.tiff", width = 15, height = 6, units = "in", res = 300)
p1+p2+p3
dev.off()


##########################################################FGF########################################################################
exprSet_NT_FGF <- FetchData(NT, vars = unique(FGF))
exprSet_NT_FGF <- as.data.frame(exprSet_NT_FGF)
exprSet_NT_FGF <- exprSet_NT_FGF[, !(colnames(exprSet_NT_FGF) %in% c("FGFBP2"))]
exprSet_NT_FGF <- cbind(exprSet_NT, exprSet_NT_FGF) %>% .[, -1]

exprSet_SM_FGF <- FetchData(SM, vars = unique(FGF))
exprSet_SM_FGF <- as.data.frame(exprSet_SM_FGF)
exprSet_SM_FGF <- exprSet_SM_FGF[, !(colnames(exprSet_SM_FGF) %in% c("FGF5",'FGF6','FGF16'))]
exprSet_SM_FGF <- cbind(exprSet_SM, exprSet_SM_FGF) %>% .[, -1]


exprSet_Noto_FGF <- FetchData(Noto, vars = unique(FGF))
exprSet_Noto_FGF <- as.data.frame(exprSet_Noto_FGF)
missing_genes <- setdiff(unique(FGF),colnames(exprSet_Noto_FGF))
exprSet_Noto_FGF[c('FGF3','FGF7','FGF18','FGF19','FGF22','FGF23','SHISA3','SHISA7','SHISA8')] <-0
exprSet_Noto_FGF <- cbind(exprSet_Noto, exprSet_Noto_FGF) %>% .[, -1]

######### NT
exprSet_NT_FGF1 <- exprSet_NT_FGF %>% 
  pivot_longer(cols = 2:38,
               names_to = 'Gene',
               values_to = 'GeneExpression')
exprSet_NT_FGF1 <- exprSet_NT_FGF1 %>%
  mutate(pseudotime_binned = cut(pseudotime,breaks = seq(0, 1, length.out = 1001),labels = FALSE))%>%
  group_by(Gene, pseudotime_binned) %>%
  summarize(GeneExpression_adjusted = mean(GeneExpression, na.rm = TRUE), .groups = "drop")

gene_matrix <- exprSet_NT_FGF1 %>%
  select(Gene, GeneExpression_adjusted, pseudotime_binned) %>%
  pivot_wider(names_from = pseudotime_binned, values_from = GeneExpression_adjusted, values_fill = 0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

n_cols <- ncol(gene_matrix)
label_positions <- seq(1, n_cols, length.out = 5)
labels <- rep("", n_cols)  # 初始化所有标记为空字符串
labels[label_positions] <- c(0, 0.25, 0.5, 0.75, 1)

col_fun <- colorRamp2(
  breaks = c(0, 0.5, 1, 1.5,2),  # 4 个数值
  colors = c("#164788", "#E2E3E1", "#FAA31C", "#FF1677",'#D33880')  # 4 个颜色
)

# 绘制热图
tiff("output/07_task3/03_heatmapwithcluster/01_3_NT_FGF_pseudotime_withcluster.tiff", 
     width = 10, height = 7, units = "in", res = 300)
Heatmap(
  gene_matrix,
  name = "Gene Expression",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "NT_FGF",
  row_title = "Genes",
  row_dend_side = "left",
  column_labels = labels, # 设置横坐标为 0, 0.25, 0.5, 0.75, 1
  column_names_rot = 0,
  heatmap_legend_param = list(
    title = "Expression",
    at = c(0, 0.5, 1, 1.5, 2),
    labels = c("0", "0.5", "1", "1.5", "2")
  ),
  width = unit(10,'cm')
)
dev.off()
######### SM
exprSet_SM_FGF1 <- exprSet_SM_FGF %>% 
  pivot_longer(cols = 2:38,
               names_to = 'Gene',
               values_to = 'GeneExpression')
exprSet_SM_FGF1 <- exprSet_SM_FGF1 %>%
  mutate(pseudotime_binned = cut(pseudotime,breaks = seq(0, 1, length.out = 1001),labels = FALSE))%>%
  group_by(Gene, pseudotime_binned) %>%
  summarize(GeneExpression_adjusted = mean(GeneExpression, na.rm = TRUE), .groups = "drop")

gene_matrix <- exprSet_SM_FGF1 %>%
  select(Gene, GeneExpression_adjusted, pseudotime_binned) %>%
  pivot_wider(names_from = pseudotime_binned, values_from = GeneExpression_adjusted, values_fill = 0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

n_cols <- ncol(gene_matrix)
label_positions <- seq(1, n_cols, length.out = 5)
labels <- rep("", n_cols)
labels[label_positions] <- c(0, 0.25, 0.5, 0.75, 1)

col_fun <- colorRamp2(
  breaks = c(0, 0.5, 1, 1.5,2),  # 4 个数值
  colors = c("#164788", "#E2E3E1", "#FAA31C", "#FF1677",'#D33880')  # 4 个颜色
)

# 绘制热图
tiff("output/07_task3/03_heatmapwithcluster/01_4_SM_FGF_pseudotime_withcluster.tiff", 
     width = 10, height = 7, units = "in", res = 300)
Heatmap(
  gene_matrix,
  name = "Gene Expression",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "SM_FGF",
  row_title = "Genes",
  row_dend_side = "left",
  column_labels = labels, # 设置横坐标为 0, 0.25, 0.5, 0.75, 1
  column_names_rot = 0,
  heatmap_legend_param = list(
    title = "Expression",
    at = c(0, 0.5, 1, 1.5, 2),
    labels = c("0", "0.5", "1", "1.5", "2")
  ),
  width = unit(10,'cm')
)
dev.off()

#########Noto
exprSet_Noto_FGF1 <- exprSet_Noto_FGF %>% 
  pivot_longer(cols = 2:38,
               names_to = 'Gene',
               values_to = 'GeneExpression')
exprSet_Noto_FGF1 <- exprSet_Noto_FGF1 %>%
  mutate(pseudotime_binned = cut(pseudotime,breaks = seq(0, 1, length.out = 1001), labels = FALSE))%>%
  group_by(Gene, pseudotime_binned) %>%
  summarize(GeneExpression_adjusted = mean(GeneExpression, na.rm = TRUE), .groups = "drop")

gene_matrix <- exprSet_Noto_FGF1 %>%
  select(Gene, GeneExpression_adjusted, pseudotime_binned) %>%
  pivot_wider(names_from = pseudotime_binned, values_from = GeneExpression_adjusted, values_fill = 0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

n_cols <- ncol(gene_matrix) 
label_positions <- seq(1, n_cols, length.out = 5)
labels <- rep("", n_cols)  # 初始化所有标记为空字符串
labels[label_positions] <- c(0, 0.25, 0.5, 0.75, 1) 

col_fun <- colorRamp2(
  breaks = c(0, 0.5, 1, 1.5,2),  # 4 个数值
  colors = c("#164788", "#E2E3E1", "#FAA31C", "#FF1677",'#D33880')  # 4 个颜色
)

# 绘制热图
tiff("output/07_task3/03_heatmapwithcluster/01_5_Noto_FGF_pseudotime_withcluster.tiff", width = 10, height = 7, units = "in", res = 300)
Heatmap(
  gene_matrix,
  name = "Gene Expression",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "NOTO_FGF",
  row_title = "Genes",
  row_dend_side = "left",
  column_labels = labels,# 设置横坐标为 0, 0.25, 0.5, 0.75, 1
  column_names_rot = 0,
  heatmap_legend_param = list(
    title = "Expression",
    at = c(0, 0.5, 1, 1.5, 2),
    labels = c("0", "0.5", "1", "1.5", "2")
  ),
  width = unit(10,'cm')
)
dev.off()

##paired
CairoTIFF("output/07_task3/03_heatmapwithcluster/FGF_paired.tiff", width = 15, height = 6, units = "in", res = 300)
p1+p2+p3
dev.off()
#######################Notch################
exprSet_NT_Notch <- FetchData(NT,vars = unique(Notch))
exprSet_NT_Notch <- as.data.frame(exprSet_NT_Notch)
exprSet_NT_Notch <- cbind(exprSet_NT,exprSet_NT_Notch) %>%  .[,-1]

exprSet_SM_Notch <- FetchData(SM,vars = unique(Notch))
exprSet_SM_Notch <- as.data.frame(exprSet_SM_Notch)
exprSet_SM_Notch <- cbind(exprSet_SM,exprSet_SM_Notch) %>%  .[,-1]

exprSet_Noto_Notch <- FetchData(Noto,vars = unique(Notch))
exprSet_Noto_Notch <- as.data.frame(exprSet_Noto_Notch)
exprSet_Noto_Notch <- cbind(exprSet_Noto,exprSet_Noto_Notch) %>%  .[,-1]
######NT
exprSet_NT_Notch1 <- exprSet_NT_Notch %>% 
  pivot_longer(cols = 2:30,
               names_to = 'Gene',
               values_to = 'GeneExpression')
exprSet_NT_Notch1 <- exprSet_NT_Notch1 %>%
  mutate(pseudotime_binned = cut(pseudotime, breaks = seq(0, 1, length.out = 1001), labels = FALSE))%>%
  group_by(Gene, pseudotime_binned) %>%
  summarize(GeneExpression_adjusted = mean(GeneExpression, na.rm = TRUE), .groups = "drop")

gene_matrix <- exprSet_NT_Notch1 %>%
  select(Gene, GeneExpression_adjusted, pseudotime_binned) %>%
  pivot_wider(names_from = pseudotime_binned, values_from = GeneExpression_adjusted, values_fill = 0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

n_cols <- ncol(gene_matrix)
label_positions <- seq(1, n_cols, length.out = 5)
labels <- rep("", n_cols)  # 初始化所有标记为空字符串
labels[label_positions] <- c(0, 0.25, 0.5, 0.75, 1)

col_fun <- colorRamp2(
  breaks = c(0, 0.5, 1, 1.5,2),  # 4 个数值
  colors = c("#164788", "#E2E3E1", "#FAA31C", "#FF1677",'#D33880')  # 4 个颜色
)

# 绘制热图
tiff("output/07_task3/03_heatmapwithcluster/01_3_NT_Notch_pseudotime_withcluster.tiff", 
     width = 10, height = 7, units = "in", res = 300)
Heatmap(
  gene_matrix,
  name = "Gene Expression",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "NT_Notch",
  row_title = "Genes",
  row_dend_side = "left",
  column_labels = labels, # 设置横坐标为 0, 0.25, 0.5, 0.75, 1
  column_names_rot = 0,
  heatmap_legend_param = list(
    title = "Expression",
    at = c(0, 0.5, 1, 1.5, 2),
    labels = c("0", "0.5", "1", "1.5", "2")
  ),
  width = unit(10, 'cm')
)
dev.off()

############SM
exprSet_SM_Notch1 <- exprSet_SM_Notch %>% 
  pivot_longer(cols = 2:30,
               names_to = 'Gene',
               values_to = 'GeneExpression')
exprSet_SM_Notch1 <- exprSet_SM_Notch1 %>%
  mutate(pseudotime_binned = cut(pseudotime, breaks = seq(0, 1, length.out = 1001), labels = FALSE))%>%
  group_by(Gene, pseudotime_binned) %>%
  summarize(GeneExpression_adjusted = mean(GeneExpression, na.rm = TRUE), .groups = "drop")

gene_matrix <- exprSet_SM_Notch1 %>%
  select(Gene, GeneExpression_adjusted, pseudotime_binned) %>%
  pivot_wider(names_from = pseudotime_binned, values_from = GeneExpression_adjusted, values_fill = 0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

n_cols <- ncol(gene_matrix)
label_positions <- seq(1, n_cols, length.out = 5)
labels <- rep("", n_cols)
labels[label_positions] <- c(0, 0.25, 0.5, 0.75, 1)

col_fun <- colorRamp2(
  breaks = c(0, 0.5, 1, 1.5,2),  # 4 个数值
  colors = c("#164788", "#E2E3E1", "#FAA31C", "#FF1677",'#D33880')  # 4 个颜色
)

# 绘制热图
tiff("output/07_task3/03_heatmapwithcluster/01_4_SM_Notch_pseudotime_withcluster.tiff", 
     width = 13, height = 7, units = "in", res = 300)
Heatmap(
  gene_matrix,
  name = "Gene Expression",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "SM_Notch",
  row_title = "Genes",
  row_dend_side = "left",
  column_labels = labels, # 设置横坐标为 0, 0.25, 0.5, 0.75, 1
  column_names_rot = 0,
  heatmap_legend_param = list(
    title = "Expression",
    at = c(0, 0.5, 1, 1.5, 2),
    labels = c("0", "0.5", "1", "1.5", "2")
  ),
  width = unit(10,'cm')
)
dev.off()

#########Noto
exprSet_Noto_Notch1 <- exprSet_Noto_Notch %>% 
  pivot_longer(cols = 2:30,
               names_to = 'Gene',
               values_to = 'GeneExpression')
exprSet_Noto_Notch1 <- exprSet_Noto_Notch1 %>%
  mutate(pseudotime_binned = cut(pseudotime,breaks = seq(0, 1, length.out = 1001), labels = FALSE))%>%
  group_by(Gene, pseudotime_binned) %>%
  summarize(GeneExpression_adjusted = mean(GeneExpression, na.rm = TRUE), .groups = "drop")

gene_matrix <- exprSet_Noto_Notch1 %>%
  select(Gene, GeneExpression_adjusted, pseudotime_binned) %>%
  pivot_wider(names_from = pseudotime_binned, values_from = GeneExpression_adjusted, values_fill = 0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

n_cols <- ncol(gene_matrix) 
label_positions <- seq(1, n_cols, length.out = 5)
labels <- rep("", n_cols)  # 初始化所有标记为空字符串
labels[label_positions] <- c(0, 0.25, 0.5, 0.75, 1) 

col_fun <- colorRamp2(
  breaks = c(0, 0.5, 1, 1.5,2),  # 4 个数值
  colors = c("#164788", "#E2E3E1", "#FAA31C", "#FF1677",'#D33880')  # 4 个颜色
)

# 绘制热图
tiff("output/07_task3/03_heatmapwithcluster/01_5_Noto_Notch_pseudotime_withcluster.tiff", width = 10, height = 7, units = "in", res = 300)
Heatmap(
  gene_matrix,
  name = "Gene Expression",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "NOTO_Notch",
  row_title = "Genes",
  row_dend_side = "left",
  column_labels = labels,# 设置横坐标为 0, 0.25, 0.5, 0.75, 1
  column_names_rot = 0,
  heatmap_legend_param = list(
    title = "Expression",
    at = c(0, 0.5, 1, 1.5, 2),
    labels = c("0", "0.5", "1", "1.5", "2")
  ),
  width = unit(10,'cm')
)
dev.off()

##paired
CairoTIFF("output/07_task3/03_heatmapwithcluster/Notch_paired.tiff", width = 15, height = 6, units = "in", res = 300)
p1+p2+p3
dev.off()


###HOXA########
exprSet_NT_A <- FetchData(NT,vars = unique(HOXA))
exprSet_NT_A <- as.data.frame(exprSet_NT_A[,c(1:10)])
exprSet_NT_A <- cbind(exprSet_NT,exprSet_NT_A) %>%  .[,-1]
exprSet_NT_A[, 2:11] <- exprSet_NT_A[, 2:11] + 0.001

exprSet_SM_A <- FetchData(SM,vars = unique(HOXA))
exprSet_SM_A <- as.data.frame(exprSet_SM_A[,c(1:9)])
exprSet_SM_A <- cbind(exprSet_SM,exprSet_SM_A) %>%  .[,-1]
exprSet_SM_A$HOXA11 <- 0
exprSet_SM_A[, 2:11] <- exprSet_SM_A[, 2:11] + 0.001


exprSet_Noto_A <- FetchData(Noto,vars = unique(HOXA))
exprSet_Noto_A <- as.data.frame(exprSet_Noto_A[,c(1:9)])
exprSet_Noto_A <- cbind(exprSet_Noto,exprSet_Noto_A) %>%  .[,-1]
exprSet_Noto_A$HOXA11 <- 0
exprSet_Noto_A[, 2:11] <- exprSet_Noto_A[, 2:11] + 0.001

############ NT
exprSet_NT_A1 <- exprSet_NT_A %>% 
  pivot_longer(cols = 2:11,
               names_to = 'Gene',
               values_to = 'GeneExpression')
exprSet_NT_A1 <- exprSet_NT_A1 %>%
  mutate(pseudotime_binned = cut(pseudotime,breaks = seq(0, 1, length.out = 1001), labels = FALSE)) %>%
  group_by(Gene, pseudotime_binned) %>%
  summarize(GeneExpression_adjusted = mean(GeneExpression, na.rm = TRUE), .groups = "drop")

gene_matrix <- exprSet_NT_A1 %>%
  select(Gene, GeneExpression_adjusted, pseudotime_binned) %>%
  pivot_wider(names_from = pseudotime_binned, values_from = GeneExpression_adjusted, values_fill = 0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

order <- c("HOXA1", "HOXA2", "HOXA3", "HOXA4", "HOXA5", "HOXA6", "HOXA7", "HOXA9", "HOXA10", "HOXA11")
gene_matrix <- gene_matrix[order, , drop = FALSE]

n_cols <- ncol(gene_matrix) 
label_positions <- seq(1, n_cols, length.out = 5)
labels <- rep("", n_cols)  # 初始化所有标记为空字符串
labels[label_positions] <- c(0, 0.25, 0.5, 0.75, 1) 

col_fun <- colorRamp2(
  breaks = c(0, 0.5, 1, 1.5,2),  # 4 个数值
  colors = c("#164788", "#E2E3E1", "#FAA31C", "#FF1677",'#D33880')  # 4 个颜色
)

tiff("output/07_task3/01_1_HOXA_pseudotime.tiff", width = 14, height = 5, units = "in", res = 300)
p1=Heatmap(
  gene_matrix,
  name = "Gene Expression",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "NT_HOXA",
  row_title = NULL,
  row_names_side = "left",
  column_labels = labels,# 设置横坐标为 0, 0.25, 0.5, 0.75, 1
  column_names_rot = 0,
  heatmap_legend_param = list(
    title = "Expression",
    at = c(0, 0.5, 1, 1.5, 2),
    labels = c("0", "0.5", "1", "1.5", "2")
  ),
  width = unit(20,'cm')
)
dev.off()

#####SM
exprSet_SM_A1 <- exprSet_SM_A %>% 
  pivot_longer(cols = 2:11,
               names_to = 'Gene',
               values_to = 'GeneExpression')
exprSet_SM_A1 <- exprSet_SM_A1 %>%
  mutate(pseudotime_binned = cut(pseudotime,breaks = seq(0, 1, length.out = 1001), labels = FALSE)) %>%
  group_by(Gene, pseudotime_binned) %>%
  summarize(GeneExpression_adjusted = mean(GeneExpression, na.rm = TRUE), .groups = "drop")

gene_matrix <- exprSet_SM_A1 %>%
  select(Gene, GeneExpression_adjusted, pseudotime_binned) %>%
  pivot_wider(names_from = pseudotime_binned, values_from = GeneExpression_adjusted, values_fill = 0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

order <- c("HOXA1", "HOXA2", "HOXA3", "HOXA4", "HOXA5", "HOXA6", "HOXA7", "HOXA9", "HOXA10", "HOXA11")
gene_matrix <- gene_matrix[order, , drop = FALSE]

n_cols <- ncol(gene_matrix) 
label_positions <- seq(1, n_cols, length.out = 5)
labels <- rep("", n_cols)  # 初始化所有标记为空字符串
labels[label_positions] <- c(0, 0.25, 0.5, 0.75, 1) 

col_fun <- colorRamp2(
  breaks = c(0, 0.5, 1, 1.5,2),  # 4 个数值
  colors = c("#164788", "#E2E3E1", "#FAA31C", "#FF1677",'#D33880')  # 4 个颜色
)

tiff("output/07_task3/01_2_HOXA_pseudotime.tiff", width = 14, height = 5, units = "in", res = 300)
p2=Heatmap(
  gene_matrix,
  name = "Gene Expression",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "SM_HOXA",
  row_title = NULL,
  row_names_side = "left",
  column_labels = labels,# 设置横坐标为 0, 0.25, 0.5, 0.75, 1
  column_names_rot = 0,
  heatmap_legend_param = list(
    title = "Expression",
    at = c(0, 0.5, 1, 1.5, 2),
    labels = c("0", "0.5", "1", "1.5", "2")
  ),
  width = unit(20,'cm')
)
dev.off()

###NOto

exprSet_Noto_A1 <- exprSet_Noto_A %>% 
  pivot_longer(cols = 2:11,
               names_to = 'Gene',
               values_to = 'GeneExpression')
exprSet_Noto_A1 <- exprSet_Noto_A1 %>%
  mutate(pseudotime_binned = cut(pseudotime,breaks = seq(0, 1, length.out = 1001), labels = FALSE)) %>%
  group_by(Gene, pseudotime_binned) %>%
  summarize(GeneExpression_adjusted = mean(GeneExpression, na.rm = TRUE), .groups = "drop")

gene_matrix <- exprSet_Noto_A1 %>%
  select(Gene, GeneExpression_adjusted, pseudotime_binned) %>%
  pivot_wider(names_from = pseudotime_binned, values_from = GeneExpression_adjusted, values_fill = 0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

order <- c("HOXA1", "HOXA2", "HOXA3", "HOXA4", "HOXA5", "HOXA6", "HOXA7", "HOXA9", "HOXA10", "HOXA11")
gene_matrix <- gene_matrix[order, , drop = FALSE]

n_cols <- ncol(gene_matrix) 
label_positions <- seq(1, n_cols, length.out = 5)
labels <- rep("", n_cols)  # 初始化所有标记为空字符串
labels[label_positions] <- c(0, 0.25, 0.5, 0.75, 1) 

col_fun <- colorRamp2(
  breaks = c(0, 0.5, 1, 1.5,2),  # 4 个数值
  colors = c("#164788", "#E2E3E1", "#FAA31C", "#FF1677",'#D33880')  # 4 个颜色
)

tiff("output/07_task3/01_3_HOXA_pseudotime.tiff", width = 14, height = 5, units = "in", res = 300)
p3=Heatmap(
  gene_matrix,
  name = "Gene Expression",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "NOTO_HOXA",
  row_title = NULL,
  row_names_side = "left",
  column_labels = labels,# 设置横坐标为 0, 0.25, 0.5, 0.75, 1
  column_names_rot = 0,
  heatmap_legend_param = list(
    title = "Expression",
    at = c(0, 0.5, 1, 1.5, 2),
    labels = c("0", "0.5", "1", "1.5", "2")
  ),
  width = unit(20,'cm')
)
dev.off()


##paired
CairoTIFF("output/07_task3/03_heatmapwithcluster/HOXA_paired.tiff", width = 30, height = 5, units = "in", res = 300)
p1+p2+p3
dev.off()

###########HOXB####################
exprSet_NT_B <- FetchData(NT,vars = unique(HOXB))
exprSet_NT_B <- as.data.frame(exprSet_NT_B[,c(1:9)])
exprSet_NT_B <- cbind(exprSet_NT,exprSet_NT_B) %>%  .[,-1]

exprSet_SM_B <- FetchData(SM,vars = unique(HOXB))
exprSet_SM_B <- as.data.frame(exprSet_SM_B[,c(1:9)])
exprSet_SM_B <- cbind(exprSet_SM,exprSet_SM_B) %>%  .[,-1]

exprSet_Noto_B <- FetchData(Noto,vars = unique(HOXB))
exprSet_Noto_B <- as.data.frame(exprSet_Noto_B[,c(1:9)])
exprSet_Noto_B <- cbind(exprSet_Noto,exprSet_Noto_B) %>%  .[,-1]

############ NT
exprSet_NT_B1 <- exprSet_NT_B %>% 
  pivot_longer(cols = 2:10,
               names_to = 'Gene',
               values_to = 'GeneExpression')
exprSet_NT_B1 <- exprSet_NT_B1 %>%
  mutate(pseudotime_binned = cut(pseudotime,breaks = seq(0, 1, length.out = 1001), labels = FALSE)) %>%
  group_by(Gene, pseudotime_binned) %>%
  summarize(GeneExpression_adjusted = mean(GeneExpression, na.rm = TRUE), .groups = "drop")

gene_matrix <- exprSet_NT_B1 %>%
  select(Gene, GeneExpression_adjusted, pseudotime_binned) %>%
  pivot_wider(names_from = pseudotime_binned, values_from = GeneExpression_adjusted, values_fill = 0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

order <-  rev(c('HOXB9','HOXB8','HOXB7','HOXB6','HOXB5','HOXB4','HOXB3','HOXB2','HOXB1'))
gene_matrix <- gene_matrix[order, , drop = FALSE]
gene_matrix <- gene_matrix[order, , drop = FALSE]

n_cols <- ncol(gene_matrix) 
label_positions <- seq(1, n_cols, length.out = 5)
labels <- rep("", n_cols)  # 初始化所有标记为空字符串
labels[label_positions] <- c(0, 0.25, 0.5, 0.75, 1) 

col_fun <- colorRamp2(
  breaks = c(0, 0.5, 1, 1.5,2),  # 4 个数值
  colors = c("#164788", "#E2E3E1", "#FAA31C", "#FF1677",'#D33880')  # 4 个颜色
)

tiff("output/07_task3/01_1_HOXB_pseudotime.tiff", width = 14, height = 5, units = "in", res = 300)
Heatmap(
  gene_matrix,
  name = "Gene Expression",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "NT_HOXB",
  row_title = NULL,
  row_names_side = "left",
  column_labels = labels,# 设置横坐标为 0, 0.25, 0.5, 0.75, 1
  column_names_rot = 0,
  heatmap_legend_param = list(
    title = "Expression",
    at = c(0, 0.5, 1, 1.5, 2),
    labels = c("0", "0.5", "1", "1.5", "2")
  ),
  width = unit(20,'cm')
)
dev.off()

#####SM
exprSet_SM_B1 <- exprSet_SM_B %>% 
  pivot_longer(cols = 2:10,
               names_to = 'Gene',
               values_to = 'GeneExpression')
exprSet_SM_B1 <- exprSet_SM_B1 %>%
  mutate(pseudotime_binned = cut(pseudotime,breaks = seq(0, 1, length.out = 1001), labels = FALSE)) %>%
  group_by(Gene, pseudotime_binned) %>%
  summarize(GeneExpression_adjusted = mean(GeneExpression, na.rm = TRUE), .groups = "drop")

gene_matrix <- exprSet_SM_B1 %>%
  select(Gene, GeneExpression_adjusted, pseudotime_binned) %>%
  pivot_wider(names_from = pseudotime_binned, values_from = GeneExpression_adjusted, values_fill = 0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

order <-  rev(c('HOXB9','HOXB8','HOXB7','HOXB6','HOXB5','HOXB4','HOXB3','HOXB2','HOXB1'))
gene_matrix <- gene_matrix[order, , drop = FALSE]

n_cols <- ncol(gene_matrix) 
label_positions <- seq(1, n_cols, length.out = 5)
labels <- rep("", n_cols)  # 初始化所有标记为空字符串
labels[label_positions] <- c(0, 0.25, 0.5, 0.75, 1) 

col_fun <- colorRamp2(
  breaks = c(0, 0.5, 1, 1.5,2),  # 4 个数值
  colors = c("#164788", "#E2E3E1", "#FAA31C", "#FF1677",'#D33880')  # 4 个颜色
)

tiff("output/07_task3/01_2_HOXB_pseudotime.tiff", width = 14, height = 5, units = "in", res = 300)
Heatmap(
  gene_matrix,
  name = "Gene Expression",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "SM_HOXB",
  row_title = NULL,
  row_names_side = "left",
  column_labels = labels,# 设置横坐标为 0, 0.25, 0.5, 0.75, 1
  column_names_rot = 0,
  heatmap_legend_param = list(
    title = "Expression",
    at = c(0, 0.5, 1, 1.5, 2),
    labels = c("0", "0.5", "1", "1.5", "2")
  ),
  width = unit(20,'cm')
)
dev.off()

###NOto

exprSet_Noto_B1 <- exprSet_Noto_B %>% 
  pivot_longer(cols = 2:10,
               names_to = 'Gene',
               values_to = 'GeneExpression')
exprSet_Noto_B1 <- exprSet_Noto_B1 %>%
  mutate(pseudotime_binned = cut(pseudotime,breaks = seq(0, 1, length.out = 1001), labels = FALSE)) %>%
  group_by(Gene, pseudotime_binned) %>%
  summarize(GeneExpression_adjusted = mean(GeneExpression, na.rm = TRUE), .groups = "drop")

gene_matrix <- exprSet_Noto_B1 %>%
  select(Gene, GeneExpression_adjusted, pseudotime_binned) %>%
  pivot_wider(names_from = pseudotime_binned, values_from = GeneExpression_adjusted, values_fill = 0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

order <-  rev(c('HOXB9','HOXB8','HOXB7','HOXB6','HOXB5','HOXB4','HOXB3','HOXB2','HOXB1'))
gene_matrix <- gene_matrix[order, , drop = FALSE]

n_cols <- ncol(gene_matrix) 
label_positions <- seq(1, n_cols, length.out = 5)
labels <- rep("", n_cols)  # 初始化所有标记为空字符串
labels[label_positions] <- c(0, 0.25, 0.5, 0.75, 1) 

col_fun <- colorRamp2(
  breaks = c(0, 0.5, 1, 1.5,2),  # 4 个数值
  colors = c("#164788", "#E2E3E1", "#FAA31C", "#FF1677",'#D33880')  # 4 个颜色
)

tiff("output/07_task3/01_3_HOXB_pseudotime.tiff", width = 14, height = 5, units = "in", res = 300)
Heatmap(
  gene_matrix,
  name = "Gene Expression",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "NOTO_HOXB",
  row_title = NULL,
  row_names_side = "left",
  column_labels = labels,# 设置横坐标为 0, 0.25, 0.5, 0.75, 1
  column_names_rot = 0,
  heatmap_legend_param = list(
    title = "Expression",
    at = c(0, 0.5, 1, 1.5, 2),
    labels = c("0", "0.5", "1", "1.5", "2")
  ),
  width = unit(20,'cm')
)
dev.off()


##paired
CairoTIFF("output/07_task3/03_heatmapwithcluster/HOXB_paired.tiff", width = 30, height = 5, units = "in", res = 300)
P1+p2+p3
dev.off()

######HOXC##############
exprSet_NT_C <- FetchData(NT,vars = unique(HOXC))
exprSet_NT_C <- as.data.frame(exprSet_NT_C[,c(1:6)])
exprSet_NT_C <- cbind(exprSet_NT,exprSet_NT_C) %>%  .[,-1]

exprSet_SM_C <- FetchData(SM,vars = unique(HOXC))
exprSet_SM_C <- as.data.frame(exprSet_SM_C[,c(1:6)])
exprSet_SM_C <- cbind(exprSet_SM,exprSet_SM_C) %>%  .[,-1]

exprSet_Noto_C <- FetchData(Noto,vars = unique(HOXC))
exprSet_Noto_C$HOXC10 <- 0
exprSet_Noto_C <- as.data.frame(exprSet_Noto_C[,c(1:6)])
exprSet_Noto_C <- cbind(exprSet_Noto,exprSet_Noto_C) %>%  .[,-1]

############ NT
exprSet_NT_C1 <- exprSet_NT_C %>% 
  pivot_longer(cols = 2:7,
               names_to = 'Gene',
               values_to = 'GeneExpression')
exprSet_NT_C1 <- exprSet_NT_C1 %>%
  mutate(pseudotime_binned = cut(pseudotime,breaks = seq(0, 1, length.out = 1001), labels = FALSE)) %>%
  group_by(Gene, pseudotime_binned) %>%
  summarize(GeneExpression_adjusted = mean(GeneExpression, na.rm = TRUE), .groups = "drop")

gene_matrix <- exprSet_NT_C1 %>%
  select(Gene, GeneExpression_adjusted, pseudotime_binned) %>%
  pivot_wider(names_from = pseudotime_binned, values_from = GeneExpression_adjusted, values_fill = 0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

order <-  rev(c('HOXC10','HOXC9','HOXC8','HOXC6','HOXC5','HOXC4'))
gene_matrix <- gene_matrix[order, , drop = FALSE]

n_cols <- ncol(gene_matrix) 
label_positions <- seq(1, n_cols, length.out = 5)
labels <- rep("", n_cols)  # 初始化所有标记为空字符串
labels[label_positions] <- c(0, 0.25, 0.5, 0.75, 1) 

col_fun <- colorRamp2(
  breaks = c(0, 0.5, 1, 1.5,2),  # 4 个数值
  colors = c("#164788", "#E2E3E1", "#FAA31C", "#FF1677",'#D33880')  # 4 个颜色
)

tiff("output/07_task3/01_1_HOXC_pseudotime.tiff", width = 14, height = 5, units = "in", res = 300)
P1=Heatmap(
  gene_matrix,
  name = "Gene Expression",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "NT_HOXC",
  row_title = NULL,
  row_names_side = "left",
  column_labels = labels,# 设置横坐标为 0, 0.25, 0.5, 0.75, 1
  column_names_rot = 0,
  heatmap_legend_param = list(
    title = "Expression",
    at = c(0, 0.5, 1, 1.5, 2),
    labels = c("0", "0.5", "1", "1.5", "2")
  ),
  width = unit(20,'cm')
)
dev.off()

#####SM
exprSet_SM_C1 <- exprSet_SM_C %>% 
  pivot_longer(cols = 2:7,
               names_to = 'Gene',
               values_to = 'GeneExpression')
exprSet_SM_C1 <- exprSet_SM_C1 %>%
  mutate(pseudotime_binned = cut(pseudotime,breaks = seq(0, 1, length.out = 1001), labels = FALSE)) %>%
  group_by(Gene, pseudotime_binned) %>%
  summarize(GeneExpression_adjusted = mean(GeneExpression, na.rm = TRUE), .groups = "drop")

gene_matrix <- exprSet_SM_C1 %>%
  select(Gene, GeneExpression_adjusted, pseudotime_binned) %>%
  pivot_wider(names_from = pseudotime_binned, values_from = GeneExpression_adjusted, values_fill = 0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

order <-  rev(c('HOXC10','HOXC9','HOXC8','HOXC6','HOXC5','HOXC4'))
gene_matrix <- gene_matrix[order, , drop = FALSE]

n_cols <- ncol(gene_matrix) 
label_positions <- seq(1, n_cols, length.out = 5)
labels <- rep("", n_cols)  # 初始化所有标记为空字符串
labels[label_positions] <- c(0, 0.25, 0.5, 0.75, 1) 

col_fun <- colorRamp2(
  breaks = c(0, 0.5, 1, 1.5,2),  # 4 个数值
  colors = c("#164788", "#E2E3E1", "#FAA31C", "#FF1677",'#D33880')  # 4 个颜色
)

tiff("output/07_task3/01_2_HOXC_pseudotime.tiff", width = 14, height = 5, units = "in", res = 300)
p2=Heatmap(
  gene_matrix,
  name = "Gene Expression",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "SM_HOXC",
  row_title = NULL,
  row_names_side = "left",
  column_labels = labels,# 设置横坐标为 0, 0.25, 0.5, 0.75, 1
  column_names_rot = 0,
  heatmap_legend_param = list(
    title = "Expression",
    at = c(0, 0.5, 1, 1.5, 2),
    labels = c("0", "0.5", "1", "1.5", "2")
  ),
  width = unit(20,'cm')
)
dev.off()

###NOto

exprSet_Noto_C1 <- exprSet_Noto_C %>% 
  pivot_longer(cols = 2:7,
               names_to = 'Gene',
               values_to = 'GeneExpression')
exprSet_Noto_C1 <- exprSet_Noto_C1 %>%
  mutate(pseudotime_binned = cut(pseudotime,breaks = seq(0, 1, length.out = 1001), labels = FALSE)) %>%
  group_by(Gene, pseudotime_binned) %>%
  summarize(GeneExpression_adjusted = mean(GeneExpression, na.rm = TRUE), .groups = "drop")

gene_matrix <- exprSet_Noto_C1 %>%
  select(Gene, GeneExpression_adjusted, pseudotime_binned) %>%
  pivot_wider(names_from = pseudotime_binned, values_from = GeneExpression_adjusted, values_fill = 0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

order <-  rev(c('HOXC10','HOXC9','HOXC8','HOXC6','HOXC5','HOXC4'))
gene_matrix <- gene_matrix[order, , drop = FALSE]

n_cols <- ncol(gene_matrix) 
label_positions <- seq(1, n_cols, length.out = 5)
labels <- rep("", n_cols)  # 初始化所有标记为空字符串
labels[label_positions] <- c(0, 0.25, 0.5, 0.75, 1) 

col_fun <- colorRamp2(
  breaks = c(0, 0.5, 1, 1.5,2),  # 4 个数值
  colors = c("#164788", "#E2E3E1", "#FAA31C", "#FF1677",'#D33880')  # 4 个颜色
)

tiff("output/07_task3/01_3_HOXC_pseudotime.tiff", width = 14, height = 5, units = "in", res = 300)
Heatmap(
  gene_matrix,
  name = "Gene Expression",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "NOTO_HOXC",
  row_title = NULL,
  row_names_side = "left",
  column_labels = labels,# 设置横坐标为 0, 0.25, 0.5, 0.75, 1
  column_names_rot = 0,
  heatmap_legend_param = list(
    title = "Expression",
    at = c(0, 0.5, 1, 1.5, 2),
    labels = c("0", "0.5", "1", "1.5", "2")
  ),
  width = unit(20,'cm')
)
dev.off()


##paired
CairoTIFF("output/07_task3/03_heatmapwithcluster/HOXC_paired.tiff", width = 30, height = 5, units = "in", res = 300)
P1+p2+p3
dev.off()

###############HOXD################
exprSet_NT_D <- FetchData(NT,vars = unique(HOXD))
exprSet_NT_D <- as.data.frame(exprSet_NT_D[,c(1:7)])
exprSet_NT_D <- cbind(exprSet_NT,exprSet_NT_D) %>%  .[,-1]
exprSet_NT_D[, 2:8] <- exprSet_NT_D[, 2:8] + 0.001

exprSet_SM_D <- FetchData(SM,vars = unique(HOXD))
exprSet_SM_D$HOXD10 <- 0
exprSet_SM_D$HOXD11 <- 0
exprSet_SM_D <- as.data.frame(exprSet_SM_D[,c(1:7)])
exprSet_SM_D <- cbind(exprSet_SM,exprSet_SM_D) %>%  .[,-1]
exprSet_SM_D[, 2:8] <- exprSet_SM_D[, 2:8] + 0.001

exprSet_Noto_D <- FetchData(Noto,vars = unique(HOXD))
exprSet_Noto_D$HOXD4 <- 0
exprSet_Noto_D <- as.data.frame(exprSet_Noto_D[,c(1:7)])
exprSet_Noto_D <- cbind(exprSet_Noto,exprSet_Noto_D) %>%  .[,-1]
exprSet_Noto_D[, 2:8] <- exprSet_Noto_D[, 2:8] + 0.001
############ NT
exprSet_NT_D1 <- exprSet_NT_D %>% 
  pivot_longer(cols = 2:8,
               names_to = 'Gene',
               values_to = 'GeneExpression')
exprSet_NT_D1 <- exprSet_NT_D1 %>%
  mutate(pseudotime_binned = cut(pseudotime,breaks = seq(0, 1, length.out = 1001), labels = FALSE)) %>%
  group_by(Gene, pseudotime_binned) %>%
  summarize(GeneExpression_adjusted = mean(GeneExpression, na.rm = TRUE), .groups = "drop")

gene_matrix <- exprSet_NT_D1 %>%
  select(Gene, GeneExpression_adjusted, pseudotime_binned) %>%
  pivot_wider(names_from = pseudotime_binned, values_from = GeneExpression_adjusted, values_fill = 0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

order <-  rev(c('HOXD11','HOXD10','HOXD9','HOXD8','HOXD4','HOXD3','HOXD1'))
gene_matrix <- gene_matrix[order, , drop = FALSE]

n_cols <- ncol(gene_matrix) 
label_positions <- seq(1, n_cols, length.out = 5)
labels <- rep("", n_cols)  # 初始化所有标记为空字符串
labels[label_positions] <- c(0, 0.25, 0.5, 0.75, 1) 

col_fun <- colorRamp2(
  breaks = c(0, 0.5, 1, 1.5,2),  # 4 个数值
  colors = c("#164788", "#E2E3E1", "#FAA31C", "#FF1677",'#D33880')  # 4 个颜色
)

tiff("output/07_task3/01_1_HOXD_pseudotime.tiff", width = 14, height = 5, units = "in", res = 300)
P1=Heatmap(
  gene_matrix,
  name = "Gene Expression",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "NT_HOXD",
  row_title = NULL,
  row_names_side = "left",
  column_labels = labels,# 设置横坐标为 0, 0.25, 0.5, 0.75, 1
  column_names_rot = 0,
  heatmap_legend_param = list(
    title = "Expression",
    at = c(0, 0.5, 1, 1.5, 2),
    labels = c("0", "0.5", "1", "1.5", "2")
  ),
  width = unit(20,'cm')
)
dev.off()

#####SM
exprSet_SM_D1 <- exprSet_SM_D %>% 
  pivot_longer(cols = 2:8,
               names_to = 'Gene',
               values_to = 'GeneExpression')
exprSet_SM_D1 <- exprSet_SM_D1 %>%
  mutate(pseudotime_binned = cut(pseudotime,breaks = seq(0, 1, length.out = 1001), labels = FALSE)) %>%
  group_by(Gene, pseudotime_binned) %>%
  summarize(GeneExpression_adjusted = mean(GeneExpression, na.rm = TRUE), .groups = "drop")

gene_matrix <- exprSet_SM_D1 %>%
  select(Gene, GeneExpression_adjusted, pseudotime_binned) %>%
  pivot_wider(names_from = pseudotime_binned, values_from = GeneExpression_adjusted, values_fill = 0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

order <-  rev(c('HOXD11','HOXD10','HOXD9','HOXD8','HOXD4','HOXD3','HOXD1'))
gene_matrix <- gene_matrix[order, , drop = FALSE]

n_cols <- ncol(gene_matrix) 
label_positions <- seq(1, n_cols, length.out = 5)
labels <- rep("", n_cols)  # 初始化所有标记为空字符串
labels[label_positions] <- c(0, 0.25, 0.5, 0.75, 1) 

col_fun <- colorRamp2(
  breaks = c(0, 0.5, 1, 1.5,2),  # 4 个数值
  colors = c("#164788", "#E2E3E1", "#FAA31C", "#FF1677",'#D33880')  # 4 个颜色
)

tiff("output/07_task3/01_2_HOXD_pseudotime.tiff", width = 14, height = 5, units = "in", res = 300)
p2=Heatmap(
  gene_matrix,
  name = "Gene Expression",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "SM_HOXD",
  row_title = NULL,
  row_names_side = "left",
  column_labels = labels,# 设置横坐标为 0, 0.25, 0.5, 0.75, 1
  column_names_rot = 0,
  heatmap_legend_param = list(
    title = "Expression",
    at = c(0, 0.5, 1, 1.5, 2),
    labels = c("0", "0.5", "1", "1.5", "2")
  ),
  width = unit(20,'cm')
)
dev.off()

###NOto

exprSet_Noto_D1 <- exprSet_Noto_D %>% 
  pivot_longer(cols = 2:8,
               names_to = 'Gene',
               values_to = 'GeneExpression')
exprSet_Noto_D1 <- exprSet_Noto_D1 %>%
  mutate(pseudotime_binned = cut(pseudotime,breaks = seq(0, 1, length.out = 1001), labels = FALSE)) %>%
  group_by(Gene, pseudotime_binned) %>%
  summarize(GeneExpression_adjusted = mean(GeneExpression, na.rm = TRUE), .groups = "drop")

gene_matrix <- exprSet_Noto_D1 %>%
  select(Gene, GeneExpression_adjusted, pseudotime_binned) %>%
  pivot_wider(names_from = pseudotime_binned, values_from = GeneExpression_adjusted, values_fill = 0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

order <-  rev(c('HOXD11','HOXD10','HOXD9','HOXD8','HOXD4','HOXD3','HOXD1'))
gene_matrix <- gene_matrix[order, , drop = FALSE]

n_cols <- ncol(gene_matrix) 
label_positions <- seq(1, n_cols, length.out = 5)
labels <- rep("", n_cols)  # 初始化所有标记为空字符串
labels[label_positions] <- c(0, 0.25, 0.5, 0.75, 1) 

col_fun <- colorRamp2(
  breaks = c(0, 0.5, 1, 1.5,2),  # 4 个数值
  colors = c("#164788", "#E2E3E1", "#FAA31C", "#FF1677",'#D33880')  # 4 个颜色
)

tiff("output/07_task3/01_3_HOXD_pseudotime.tiff", width = 14, height = 5, units = "in", res = 300)
p3=Heatmap(
  gene_matrix,
  name = "Gene Expression",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "NOTO_HOXD",
  row_title = NULL,
  row_names_side = "left",
  column_labels = labels,# 设置横坐标为 0, 0.25, 0.5, 0.75, 1
  column_names_rot = 0,
  heatmap_legend_param = list(
    title = "Expression",
    at = c(0, 0.5, 1, 1.5, 2),
    labels = c("0", "0.5", "1", "1.5", "2")
  ),
  width = unit(20,'cm')
)
dev.off()


##paired
CairoTIFF("output/07_task3/03_heatmapwithcluster/HOXD_paired.tiff", width = 30, height = 5, units = "in", res = 300)
P1+p2+p3
dev.off()


********************************************
******************violin_plot***************
********************************************


library(Seurat)
library(patchwork)
library(ggplot2)
######
source('resource/color_map.R')
NTSM <- qs::qread(file = 'output/02_output/hTEM_v3_int_obj.qs')


c <- subset(NTSM,detailed_celltype_2nd %in% c('Notochord','NMP-Neural','NMP-Meso','pPSM','aPSM',
                                              'E-SM','M-SM','Sclerotome','Syndetome','Endotome',
                                              'Myogenic Progenitors','Dermomyotome','Caudal Neural Plate',
                                              'Early NT','Early FP','Neural Tube','Floor Plate'))


vangl1 <- subset(c, subset = VANGL1 >0)
vangl2 <- subset(c, subset = VANGL2 >0)

VlnPlot(vangl1, features =c("VANGL1"), split.by = 'detailed_celltype_2nd',split.plot = TRUE, pt.size = 1)+
  scale_fill_manual(values=v3_deCell2)+
  stat_summary(fun = "mean", geom = "crossbar", size = 0.5,color = "black",position = position_dodge(width = 0.9))+
  guides(fill = guide_legend(override.aes = list(color = NA)))
ggsave(filename = 'output/10_task11/01_vangle1.pdf',width = 15.4,height = 5,dpi = 300)

VlnPlot(vangl2, features =c("VANGL2"), split.by = 'detailed_celltype_2nd',split.plot = TRUE, pt.size = 1)+
  scale_fill_manual(values=v3_deCell2)+
  stat_summary(fun = "mean", geom = "crossbar", size = 0.5,color = "black",position = position_dodge(width = 0.9))+
  guides(fill = guide_legend(override.aes = list(color = NA)))
ggsave(filename = 'output/10_task11/02_vangle2.pdf',width = 15.4,height = 5,dpi = 300)


********************************************
******************feature_plot**************
********************************************

library(Seurat)
library(tidyverse)
library(ggplot2)
library(Nebulosa)

#######nerual tube##############
NT <- qs::qread(file = 'output/02_output/01_NT_DV.qs')

Markers <- c('SOX2', 'TBXT', 'NKX1-2', 'NEFM', 'FOXJ1', 'CDX1', 'CDX2', 'MSX1', 'MIXL1',
             'LEFTY1', 'BMP7', 'KRT8', 'KRT19', 'CDH1', 'TIMP3', 'SEMA3C', 'FGF8', 'FGF17', 
             'FGF18', 'FGF19', 'FGFBP3', 'FGFR2', 'WNT3A', 'WNT5A', 'WNT5B', 'SFRP2', 
             'LEF1', 'CTNNB1', 'DKK1', 'DKK3', 'DLL3', 'CCND1', 'DTX4', 'BMP2', 
             'BMP4', 'BMP7', 'ID4', 'FST', 'FZD7', 'JAG1', 'CITED4', 'CITED2', 'TGFBI', 
             'DCHS1', 'CELSR1', 'DBX1', 'DBX2', 'WNT1', 'PAX6', 'PAX7', 'IRX3', 'IRX5', 
             'TUBB3', 'GRHL2', 'GRHL3', 'NKX6-1', 'NKX6-2', 'OLIG1', 'OLIG2', 'OLIG3', 
             'FOXA1', 'FOXA2', 'SHH', 'CDH1', 'CDH2', 'CDH4', 'CDH7', 'SOX1', 'SOX21', 
             'ISL1', 'NEUROG2', 'CRABP1')

library(Nebulosa)
for (marker in Markers) {
  filename <- paste0("output/08_task4/", marker, ".tiff")
  tiff(filename = filename, width = 5.5, height = 3.3, units = "in", res = 300)
  print(plot_density(NT, features = marker, reduction = 'umap.bbknn'))
  dev.off()
}

##a whole one
tiff(filename = 'output/08_task4/03_NT_Nebulosa.tiff', width = 5500, height = 3300, units = 'px')
plot_density(NT, features = Markers, reduction = 'umap.bbknn')+ plot_layout(ncol = 10)
dev.off()

################################################SM#######################################################
SM <- qs::qread(file = 'output/02_output/01_SM_DV.qs')

Markers <- c(
  'MSX1', 'RIPPLY2', 'TBX6', 'LFNG', 'MESP2', 'MESP1', 'TBX18', 'UNCX',
  'CER1', 'MESP2', 'SEMA3C', 'NEFM', 'DLL1', 'DLL3',
  'BMP4', 'BMP5', 'BMP7', 'FST',
  'KRT18', 'KRT8', 'EPHA1', 'EPHA4', 'WNT5A', 'WNT5B', 'FRZB', 'SPRY2', 
  'DACT2', 'RSPO3', 'EYA1', 'CRABP2', 'DKK1', 'DKK3', 'SOX2', 
  'PAX3', 'MEOX1', 'MEIS1', 'SIX1', 'LEF1', 'TCF15', 'FOXC1', 'FOXC2',
  'DCHS1', 'CELSR1', 'ALDH1A2', 'RARA', 'RXRG', 'DHRS3', 
  'MYF5', 'MYF6', 'PAX7', 'SCX', 'LOX', 'NEB', 'TTN', 'TWIST1', 'PAX1',
  'PAX9', 'NKX3-2', 'GLI1', 'TGFBI', 'SNAI1',
  'KDR', 'ETV2', 'CLDN5', 'CDH5'
)


library(Nebulosa)
for (marker in Markers) {
  filename <- paste0("output/08_task4/03_SM_Nebulosa/", marker, ".tiff")
  tiff(filename = filename, width = 5.5, height = 3.3, units = "in", res = 300)
  print(plot_density(SM, features = marker, reduction = 'umap.bbknn'))
  dev.off()
  message("Saved: ", filename)  # 打印保存成功信息
}

tiff(filename = 'output/08_task4/03_SM_Nebulosa/EBF2.tiff', width = 5.5, height = 3.3, units = "in", res = 300)
print(plot_density(SM, features = 'EBF2', reduction = 'umap.bbknn'))
dev.off()
##a whole one
tiff(filename = 'output/08_task4/03_SM_Nebulosa.tiff', width = 5500, height = 3300, units = 'px')
plot_density(SM, features = Markers, reduction = 'umap.bbknn')+ plot_layout(ncol = 10)
dev.off()

p1=plot_density(SM, features = 'MSX1', reduction = 'umap.bbknn',pal = 'viridis')
p2=plot_density(SM, features = 'MSX1', reduction = 'umap.bbknn',pal = 'magma')
p3=plot_density(SM, features = 'MSX1', reduction = 'umap.bbknn',pal = 'cividis')
p4=plot_density(SM, features = 'MSX1', reduction = 'umap.bbknn',pal = 'inferno')
p5=plot_density(SM, features = 'MSX1', reduction = 'umap.bbknn',pal = 'plasma')
p1+p2+p3+p4+p5
####Noto########
Noto <- qs::qread(file = 'output/02_output/01_Noto.qs')

GOI <- c("TBXT", "SOX2", "SCG3", "KRT18", "KRT19", "TIMP3", "TPPP3", "NOTO", "FOXA1", "FOXA2", "SHH",
         "NOG", "DKK1",'WNT3A', "WNT5A", "WNT5B", "WNT8A", "WNT11", "BMP2", "BMP4", "BMP7", "FST", "BAMBI", 
         "TGFB1", "TGFB2", "TGFB3", "NODAL", "EGF", "TGFBI", "SPRY2", "NEFM", "LEFTY1", "LEFTY2", 
         "SOX9", "FOXJ1", "CHRD", "CDX2", "HES7", "CDH1", "SEMA3C", "SEMA3E", "FGF8", "FGF17", 
         "CYP26A1", "CYP26C1", "CYP26B1", "VANGL1", "VANGL2", "NFIA", "RFX2", "SMAD3", "NFE2L3", 
         "MNX1", "SIX5", "FOS", "FOSL1", "LHX9", "SP5", "SERPINE2", "C1orf189", "LRRIQ1", "SPATA17", 
         "FLT4", "TFF3", "GLIS1", "GLI1", "SDR16C5", "NQO1", "TCTEX1D1",'GSC','EOMES','LEF1')

tiff(filename = "output/02_output/Figure/01_NOTO_GOI.tiff", 
     width = 30, height = 33, units = "in",res = 300)
FeaturePlot(Noto, features = unique(GOI), pt.size = 0.8,reduction = "umap.bbknn",
            ncol = 8, order = TRUE)
dev.off()

DimPlot(Noto, reduction = 'umap.bbknn',group.by = 'orig.ident',split.by = 'orig.ident')
ggsave(filename = 'output/02_output/Figure/01_NOTO_dimplot_split.tiff',height = 3.93, width = 13.51,dpi = 300)


********************************************
****************scenic_plot*****************
********************************************

******01.rss_heatmap

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

rss <- readRDS('scenic/savedobj/scenic_rss.rds')


tiff('RSS_heatmap.tiff', width = 10, height = 80, res = 300, units = 'in')
print(Heatmap(rss))
dev.off()

******02.circular_heatmap
# r-scenic
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
library(corrplot)
library(patchwork)

cluster_to_rctd <- c('28' = 'Notochord',
  '31' = 'NMP-Neural',
  '5' = 'NMP-Meso',
  '16' = 'Caudal Neural Plate',
  '12' = 'Neural Tube',
  '32' = 'Floor Plate',
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
  '3' = 'SM',
  '10' = 'SM',
  '15' = 'SM',
  '27' = 'SM',
  '8' = 'SM',
  '18' = 'SM',
  '21' = 'SM',
  '7' = 'SM',
  '22' = 'SM',
  '24' = 'Sclerotome',
  '4' = 'Syndetome',
  '6' = 'Endotome',
  '14' = 'Myogenic Progenitors',
  '11' = 'Myogenic Progenitors',
  '34' = 'Dermomyotome',
  '26' = 'Endothelial Cells',
  '36' = 'Neuroectoderm',
  '30' = 'Surface Ectoderm',
  '13' = 'Caudal Meso Progenitors',
  '9' = 'Intermediate Meso',
  '35' = 'Neurons')

srt <- readRDS("finalized_hTLO/01_NTSM_integrated_decontX_annotated_third_rctd.rds")

cellInfo <- readRDS("scenic/script/int/cellInfo.Rds")

cellInfo['detailed_rctd'] <- cluster_to_rctd[as.character(cellInfo$seuratCluster)]

exprMat_filtered <- readRDS("scenic/savedobj/checkpoint_2_exprMat_filtered.rds")
scenicOptions <- readRDS("scenic/script/int/scenicOptions.Rds")

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$celltype),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

# Figure 4 circular heatmap #
regulons <- c(
  'T', 'SOX2', 'NOTO', 'NFE2L3', 'ANXA1', 'ZBTB33', 'CDX2', 'CDX1', 'FOXJ1', 'RFX2', 
  'SALL1', 'SMAD3', 'EVX1', 'RARG', 'HOXB1', 'SP2', 
  'MIXL1', 'NKX1-2', 'RARG', 'POU5F1', 'MSGN1', 'LEF1', 'STAT1', 
  'MESP2', 'MYCN', 'SMARCA4', 'ZIC3', 'FOXC2', 'MEOX1', 'EZH2',
  'PAX9', 'PRRX2', 'FOXD1', 'NR2F6', 'TWIST1', 'TWIST2', 'NR5A2', 'ZNF124', 'STAT6',
  'MYF5', 'MYF6', 'CREB3L2', 'SMARCA4', 
  'PAX6', 'IRX3', 'NKX6-2', 'SOX3', 'POU3F1', 'NFAT2C', 'GRHL2', 'ERG1', 'THRB', 'SOX21'
)

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$detailed_rctd),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

regulon_full_name <- c()
for (r in regulons){
    for (j in rownames(regulonActivity_byCellType_Scaled)){
        if (grepl(r, j) == TRUE){
            regulon_full_name <- c(regulon_full_name, j)
        }
    }
}
regulon_full_name <- unique(regulon_full_name)

keep_clusters <- c(
  "Notochord", 
  "NMP-Meso", 
  "NMP-Neural", 
  "pPSM", 
  "aPSM", 
  "SM",
  'Sclerotome',
  'Syndetome',
  'Endotome',
  "Dermomyotome", 
  "Myogenic Progenitors", 
  "Neural Tube", 
  "Floor Plate"
)

circ_mat <- regulonActivity_byCellType_Scaled[regulon_full_name, keep_clusters]
cleaned_rownames <- sub("_.*", "(+)", sub(" .*", "", rownames(circ_mat)))
rownames(circ_mat) <- cleaned_rownames

regulons_to_remove <- c("ETV6", 'STAT3(+)', 'TFAP2A(+)', 'STAT5A(+)', 'ETS2', 'ETV2', 'ONECUT1', 'ONECUT2', 'TBP(+)', 'TAL1', 'GATA2', 'GATA6')
regulons_to_keep <- setdiff(rownames(circ_mat), regulons_to_remove)

circ_mat <- circ_mat[regulons_to_keep, ]

tiff("/home/datastore/brian/16_dalton/20250218/figure_plotting/figures/20250911_figures/revised_fig4_circle_heatmap.tiff", width=12, height=8, res=300, units='in')

interval <-  1.18
start <- 8.35

# Define the ranges and corresponding colors
rect_data <- data.frame(
  y_start = seq(from = start, by = interval, length.out = length(keep_clusters)),
  y_end = seq(from = start+interval, by = interval, length.out = length(keep_clusters)),
  col = rev(rctd_colors[keep_clusters])
)

circos.clear()

lower_val = min(circ_mat)
upper_val = max(circ_mat)

col_fun1 = colorRamp2(c(lower_val, upper_val), c("white", "firebrick"))

circos.par(start.degree = 0, gap.degree = 25)
circos.heatmap(circ_mat, dend.side = "inside", col = col_fun1, 
    dend.track.height = 0.2, rownames.side = 'outside', rownames.cex = 0.8,
    track.height = 0.4)
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  # Loop through each row of rect_data
  for (i in 1:nrow(rect_data)) {
    circos.rect(
      -6 + convert_x(0, "mm"), rect_data$y_start[i], 
      -6 + convert_x(1, "mm"), rect_data$y_end[i],
      col = rect_data$col[i], border = NA
    )
  }
}, bg.border = NA)

lgd_celltype = Legend(at = keep_clusters, 
    type = "grid", 
    legend_gp = gpar(fill = as.vector(rctd_colors[keep_clusters]), lwd = 2), 
    title_position = "leftcenter-rot", 
    title = "")
lgd_zscore = Legend(at = c(lower_val, upper_val), col_fun = col_fun1, labels = c("Low", "High"),
    title_position = "leftcenter-rot", title = "Z-score")

draw(lgd_celltype, x = unit(251, "mm"), y = unit(20, "mm"), just = c("left", "bottom"))
draw(lgd_zscore, x = unit(246, "mm"), y = unit(80, "mm"), just = c("left", "bottom"))
dev.off()

circos.clear()

*****03.scenic_regulon_correlation_analysis

# r-scenic
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
library(corrplot)
library(patchwork)

cluster_to_cell_type <- c('28' = 'Noto',
  '31' = 'NMP-Neural',
  '5' = 'NMP-Meso',
  '16' = 'Caud. NP',
  '12' = 'E-NT',
  '32' = 'E-FP',
  '1' = 'NT',
  '20' = 'FP',
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
  '24' = 'L-SM',
  '4' = 'M-SM',
  '6' = 'M-SM',
  '14' = 'L-SM',
  '11' = 'L-SM',
  '34' = 'L-SM',
  '26' = 'EC',
  '36' = 'NE',
  '30' = 'SE',
  '13' = 'Caud. Meso',
  '9' = 'Int. Meso',
  '35' = 'Neuron')

cellInfo <- readRDS("/home/datastore/brian/16_dalton/20250218/scenic/script/int/cellInfo.Rds")

cellInfo['celltype'] <- cluster_to_cell_type[as.character(cellInfo$seuratCluster)]

exprMat_filtered <- readRDS("/home/datastore/brian/16_dalton/20250218/scenic/savedobj/checkpoint_2_exprMat_filtered.rds")
scenicOptions <- readRDS("/home/datastore/brian/16_dalton/20250218/scenic/script/int/scenicOptions.Rds")

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$celltype),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

pearson_corr_scaled <- cor(t(regulonActivity_byCellType_Scaled), method = "pearson")

# Heatmap 
corrplot(pearson_corr_scaled, col = COL2('RdYlBu'), order = 'hclust', addrect = 12, method = 'color')


regulons <- readRDS('scenic/script/int/2.6_regulons_asGeneSet.Rds')

# Identify regulon groups
hc <- hclust(as.dist(1-pearson_corr_scaled), method = "complete")
k <- 12  # Number of clusters
clusters <- cutree(hc, k = k)

clusters <- clusters %>% data.frame()
colnames(clusters) <- 'cluster'

table(clusters$cluster)
clusters['recluster'] <- 0
clusters[clusters$cluster == 9, 'recluster'] <- 1
clusters[clusters$cluster == 11, 'recluster'] <- 2
clusters[clusters$cluster == 10, 'recluster'] <- 3
clusters[clusters$cluster == 5, 'recluster'] <- 4
clusters[clusters$cluster == 7, 'recluster'] <- 5
clusters[clusters$cluster == 4, 'recluster'] <- 6
clusters[clusters$cluster == 8, 'recluster'] <- 7
clusters[clusters$cluster == 6, 'recluster'] <- 8
clusters[clusters$cluster == 2, 'recluster'] <- 9
clusters[clusters$cluster == 3, 'recluster'] <- 10
clusters[clusters$cluster == 12, 'recluster'] <- 11
clusters[clusters$cluster == 1, 'recluster'] <- 12


gene_lists <- c()
for (c in sort(unique(clusters$recluster))){
    regulon_list <- clusters %>% filter(recluster == c)
    regulon_list <- rownames(regulon_list)
    regulon_list <- gsub("\\s*\\(.*\\)", "", regulon_list)
    gene_list <- c()
    for (r in regulon_list){
        gene_list <- c(gene_list, regulons[[r]])
    }
    gene_list <- unique(gene_list)
    gene_list <- gsub("^T$", "TBXT", gene_list)
    gene_lists[[c]] <- gene_list
}

srt <- readRDS("finalized_hTLO/01_NTSM_integrated_decontX_annotated_third_rctd.rds")
srt <- AddModuleScore(srt, gene_lists, name='ModuleScores')

for (c in 1:12){
    clusters %>% filter(recluster == c)
    regulon_list <- clusters %>% filter(recluster == c)
    regulon_list <- sort(rownames(regulon_list))
    ncol_use <- ceiling(sqrt(length(regulon_list)))
    tiff(paste0('module_',c,'_regulons.tiff'), width = ncol_use*5, height = ncol_use*5, res = 300, units = 'in')
    print(FeaturePlot(srt, features = regulon_list, reduction = 'umap.bbknn', order = T, ncol = ncol_use))
    dev.off()
}


*****0.4.igraph

# r-scenic
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
library(corrplot)
library(patchwork)
library(igraph)
library(ggrepel)

GRNBoost_output <- read.delim("/home/datastore/brian/16_dalton/20250218/scenic/savedobj/intermediate_for_grnboost/output/grn_output.tsv", header=FALSE)
colnames(GRNBoost_output) <- c("TF","Target","weight")
cleaned_rownames <- replace_T_with_TBXT(GRNBoost_output$Target)
GRNBoost_output$Target <- cleaned_rownames

regulons <- c('NKX1-2', 'NOTO', 'TBX6', 'SOX2', 'SIX1', 'ERG', 'FOXA1', 'FOXA2', 'NKX6-2', 'PAX6', 'MYF5')

plot_list <- c()
n <- 1
for (regulon_use in regulons){
    grn_df <- GRNBoost_output %>% filter(TF == regulon_use)

    shuffled_grn_df <- grn_df[sample(nrow(grn_df)), ]

    grn_graph <- graph_from_data_frame(d = shuffled_grn_df, directed = TRUE)
    genes_to_label <- grn_df %>% head(15)
    genes_to_label <- c(regulon_use, genes_to_label$Target)

    # Label only transcription factors (TFs) and hide other labels
    V(grn_graph)$label <- ifelse(V(grn_graph)$name %in% genes_to_label, V(grn_graph)$name, NA)

    node_df <- data.frame(
    name = V(grn_graph)$name,
    type = ifelse(V(grn_graph)$name %in% shuffled_grn_df$source, "TF", "Target")  # Label TFs vs. Targets
    )

    # Edge data
    edge_df <- as_data_frame(grn_graph, what = "edges")

    # Get layout (e.g., Fruchterman-Reingold layout)
    layout <- layout_on_sphere(grn_graph) %>% data.frame()
    colnames(layout) <- c('x2', 'y2')
    # Add x, y coordinates to node_df
    node_df$x <- layout[, 1]
    node_df$y <- layout[, 2]
    layout$x1 <- node_df[node_df$name == regulon_use, ]$x
    layout$y1 <- node_df[node_df$name == regulon_use, ]$y
    node_df['type'] <- ifelse(node_df$name == regulon_use, 'TF', 'Target')
    node_df['label'] <- ifelse(node_df$name %in% genes_to_label, node_df$name, NA)
    node_df['color'] <- ifelse(node_df$type == 'TF', '#009285', 'black')

    ggplot() +
    # Add edges (lines between nodes)
    geom_segment(data = layout, aes(x = x1, 
                                    y = y1, 
                                    xend = x2, 
                                    yend = y2),
                color = "#FAA31C", size = 0.5, arrow = arrow(length = unit(0, "cm"))) +
    # Add nodes (points)
    geom_point(data = node_df, aes(x = x, y = y), color = 'black', size = 2) +
    # Add labels with ggrepel
    geom_point(data = node_df[node_df$type == 'TF', ], aes(x = x, y = y), color = '#009285', size = 5) +
    geom_label_repel(data = node_df, aes(x = x, y = y, label = label, color=color), 
                    size = 6, box.padding = 0.5, max.overlaps = 50) +
    scale_color_identity()+
    # Customize theme
    theme_minimal() +
    labs(color = "Node Type") +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) -> p1
    plot_list[[n]] <- p1
    n <- n + 1
}

tiff('scenic/output/igraph_GRN/revised_graph_grn.tiff', width=40, height=30, res=300, units='in')
print(wrap_plots(plot_list, ncol = 4))
dev.off()


********************************************
*******v3_mouse_monkey_comparison***********
********************************************
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyr)              
  library(reticulate)
  library(Matrix)
  library(ComplexHeatmap)
  library(circlize)
  library(ClusterFoldSimilarity)
  library(BiocParallel)
  library(ggplot2)
})
use_assay_if_present <- function(obj, assay) {
  if (assay %in% Assays(obj)) DefaultAssay(obj) <- assay
  obj
}

subset_by_features <- function(obj, feats) {
  feats_use <- intersect(feats, rownames(obj))
  if (length(feats_use) == 0) stop("No overlapping features found.")
  subset(obj, features = feats_use)
}

prep_obj <- function(obj, idents_keep = NULL, feats_keep = NULL, prefer_assay = "decontX") {
  obj <- use_assay_if_present(obj, prefer_assay)
  if (!is.null(idents_keep)) obj <- subset(obj, idents = idents_keep)
  if (!is.null(feats_keep))  obj <- subset_by_features(obj, feats_keep)
  obj
}

normalize_variable_scale <- function(obj) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj, features = VariableFeatures(obj))
  obj
}

var_only <- function(obj) obj[VariableFeatures(obj), , drop = FALSE]

run_cfs <- function(obj_list, sample_names, topN = Inf, nSubsampling = 24, workers = 10) {
  BiocParallel::register(BiocParallel::MulticoreParam(workers = workers))
  clusterFoldSimilarity(
    scList = obj_list,
    sampleNames = sample_names,
    topN = topN,
    nSubsampling = nSubsampling,
    parallel = TRUE
  )
}

rescale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(x*0)  
  (x - rng[1]) / (rng[2] - rng[1])
}

rescale_blocks <- function(mat, blocks) {
  for (b in blocks) {
    r <- b$rows; c <- b$cols
    mat[r, c] <- t(apply(mat[r, c, drop = FALSE], 1, rescale01))
  }
  mat
}

combined.sct <- readRDS("01_NTSM_integrated_decontX_up2.rds")
Idents(combined.sct) <- "celltype"

features <- read.csv("features.csv", check.names = FALSE)
if (!"human" %in% colnames(features))
  stop('features.csv must contain a "human" column of gene names.')
human_feats <- unique(na.omit(as.character(features$human)))

subtype_human <- c(
  "Noto","NMP-Neural","NMP-Meso","Caud.NP","E-NT","E-FP","NT","FP","pPSM",
  "aPSM","E-SM","M-SM","L-SM","Caud.Meso","Neuron","NE","SE","Int.Meso","EC"
)

human <- prep_obj(combined.sct, idents_keep = subtype_human, feats_keep = human_feats)

mouse <- readRDS("MouseGastrulationData_E7_E8.5_subset_unique.rds")
decont_map <- read.csv("decontX_remap_mouse_feature.csv", check.names = FALSE)
mouse_feats <- unique(na.omit(as.character(decont_map[[1]])))

mouse_subtypes <- c(
  "Notochord","Caudal epiblast","Caudal Mesoderm","Somitic mesoderm","NMP",
  "Paraxial mesoderm","Forebrain/Midbrain/Hindbrain","Spinal cord","Primitive Streak",
  "Anterior Primitive Streak","Intermediate mesoderm","Surface ectoderm","Endothelium"
)

mouse <- prep_obj(mouse, idents_keep = mouse_subtypes, feats_keep = mouse_feats)

monkey <- readRDS("03_human_monkey_CS7_CS10.rds")
monkey_subtypes <- c(
  "Anterior Primitive Streak","anterior PSM","Endothelium","Epithelium","Hindbrain",
  "Forebrain/Midbrain","Intermediate Meso","Neurons","NMP-Neural","NMP-Meso","Notochord",
  "posterior PSM","Primitive Streak","Somite","Spinal Cord"
)

monkey <- prep_obj(monkey, idents_keep = monkey_subtypes, feats_keep = mouse_feats, prefer_assay = "RNA")

levs_monkey <- levels(Idents(monkey))
if ("Forebrain/Midbrain" %in% levs_monkey) monkey <- RenameIdents(monkey, "Forebrain/Midbrain" = "FB/MB/HB")
if ("Hindbrain"         %in% levs_monkey) monkey <- RenameIdents(monkey, "Hindbrain"         = "FB/MB/HB")

three_list <- list(human, mouse, monkey)
three_list_proc <- lapply(three_list, normalize_variable_scale)
three_list_var  <- lapply(three_list_proc, var_only)

sim_human_mouse_monkey <- run_cfs(three_list_var, sample_names = c("human","mouse","monkey"))
write.csv(sim_human_mouse_monkey, "similarity.table.human.mouse.monkey-up2.csv")

self_cfs <- function(obj, label_prefix, out_csv) {
  lst  <- list(obj, obj)                                 # two replicates
  lstp <- lapply(lst, normalize_variable_scale)
  lstv <- lapply(lstp, var_only)
  sim  <- run_cfs(lstv, sample_names = c(paste0(label_prefix,"1"),
                                         paste0(label_prefix,"2")))
  write.csv(sim, out_csv)
  invisible(sim)
}

sim_human_self <- self_cfs(human,  "human",  "similarity.table.human.self.csv")
sim_mouse_self <- self_cfs(mouse,  "mouse",  "similarity.table.mouse.self.csv")
sim_monkey_self<- self_cfs(monkey, "monkey", "similarity.table.monkey.self.csv")

heatmap_csv <- "corr_final_heatmap2.csv"
if (!file.exists(heatmap_csv)) {
  message("Heatmap input CSV not found: ", heatmap_csv)
} else {
  mat <- as.matrix(read.csv(heatmap_csv, row.names = 1, check.names = FALSE))

  blocks <- list(
    list(rows = 1:19,  cols = 1:19),
    list(rows = 20:46, cols = 1:19),
    list(rows = 1:19,  cols = 20:46),
    list(rows = 20:32, cols = 20:32),
    list(rows = 33:46, cols = 33:46),
    list(rows = 20:32, cols = 33:46),
    list(rows = 33:46, cols = 20:32)
  )
  blocks <- lapply(blocks, function(b) {
    b$rows <- intersect(b$rows, seq_len(nrow(mat)))
    b$cols <- intersect(b$cols, seq_len(ncol(mat)))
    b
  })
  mat <- rescale_blocks(mat, blocks)

  mat <- t(scale(t(mat), center = TRUE, scale = TRUE))

  col_fun <- colorRamp2(c(0, 0.5, 1), c("#98D8EF","white","#DE3163"))

  pdf("heatmap_corr_final_inid_scale_cluster_upname.pdf", width = 10, height = 10)
  print(Heatmap(
    mat,
    cluster_rows   = TRUE,
    cluster_columns= TRUE,
    col            = col_fun,
    use_raster     = FALSE,
    row_names_side = "left",
    row_names_gp   = gpar(fontsize = 6),
    column_names_gp= gpar(fontsize = 6),
    show_row_dend  = TRUE
  ))
  dev.off()
}

********************************************
*****************sankey_plot****************
********************************************
library(dplyr)
library(ggalluvial)
library(ggplot2)

#meta_monkey <- read.csv("03_v3_to_humanMonkey_metadata_edit.csv",
                        stringsAsFactors = FALSE) %>%
  rename(
    cell                   = cell.id,            # your barcode column
    human                  = celltype,        # human labels
    monkey_stage           = knn.pred.stage#  # monkey prediction
  )
#meta_monkey=meta_monkey[meta_monkey$knn.pred.labels!="FB/MB/HB",]

meta_mouse <- read.csv("02_v3_to_Mouse_metadata.csv",
                       stringsAsFactors = FALSE) %>%
  rename(
    cell                  = cell.id,
    human2                = celltype,       # same human labels
    mouse_stage           = knn.pred.stage#mouse prediction
  )

#meta_mouse=meta_mouse[meta_mouse$knn.pred.labels!="Forebrain/Midbrain/Hindbrain",]

meta_monkey$human[which(meta_monkey$human %in% c("E-SM","L-SM","M-SM"))]="SM"

meta_mouse$human2[which(meta_mouse$human2 %in% c("E-SM","L-SM","M-SM"))]="SM"

# 2. Inner‐join on cell → only keep barcodes present in both

combined <- inner_join(meta_monkey, meta_mouse, by = "cell")

# 3. (Optional) sanity check that human == human2
stopifnot(all(combined$human == combined$human2))

# 4. Filter to your three human types
combined <- combined %>%
  filter(human %in% c("Caud.NP","E-NT","E-FP","Neuron"))#)c("pPSM","Int.Meso","Caud.Meso","EC"))#c("SM", "aPSM", "NT","FP"))

# 5. Build the Sankey count table
sankey_df <- combined %>%
  count(monkey_stage, human, mouse_stage, name = "n") %>%
  rename(
    monkey = monkey_stage,
    mouse  = mouse_stage
  )

sankey_df <- combined %>%
  count(monkey_stage, human, mouse_stage, name = "n") %>%
  ungroup() %>%
  mutate(pct = n / sum(n) * 100)
  


# 6. (Optional) set factor levels to control vertical ordering
sankey_df <- sankey_df %>%
  mutate(
    monkey = factor(monkey_stage,levels = c("H-CS7","H-CS8", "H-CS10",
                               "M-CS8", "M-CS9", "M-CS11")),
    human  = factor(human,  levels = c("Caud.NP","E-NT","E-FP","Neuron")),#c("pPSM","Int.Meso","EC")),
    mouse  = factor(mouse_stage,  levels = c("E7.0",  "E7.25", "E7.5",  "E7.75", "E8.0",  "E8.25", "E8.5"))
  )


# 1. Prep & melt
sankey2 <- sankey_df %>%
  rename(
    monkey = monkey_stage,
    mouse  = mouse_stage
  ) %>%
  mutate(id = row_number())

long <- gather_set_data(sankey2, 1:3)

# 2. Build your fill map
monkey_lvls <- unique(long$y[long$x == 1])
mouse_lvls  <- unique(long$y[long$x == 3])
human_map   <- c(
  "SM"   = "#00879E",
  "aPSM" = "#AA60C8",
  "NT"   = "#EF9651",
  "FP"   = "#A6AEBF"
)

fill_map <- c(
  setNames(rep("#F4D793", length(monkey_lvls)), monkey_lvls),
  human_map,
  setNames(rep("#A0C878", length(mouse_lvls)), mouse_lvls)
)


monkey_order <- c("H-CS7","H-CS8","H-CS10","M-CS8","M-CS9","M-CS11")
human_order  <- c("FP","NT","SM","aPSM")
mouse_order  <- c("E7.0","E7.25","E7.5","E7.75","E8.0","E8.25","E8.5")

# 3. Re‐level the “y” factor so ggforce will respect those positions
long <- long %>%
  mutate(
    y = factor(
      y,
      levels = c(monkey_order, human_order, mouse_order)
    )
  )

# 3. Plot
p=ggplot(long,
       aes(x     = x,        # axis name: "monkey"/"human"/"mouse"
           id    = id,       # links flows
           split = y,        # <— must map your stratum column here
           value = pct)) +   # controls ribbon width
  # ribbons stay colored by human stage
  geom_parallel_sets(aes(fill = human),
                     axis.width = 0.2,
                     sep        = 0.05,alpha = 0.7)  +
  # axis blocks colored by the stratum in `y`
  geom_parallel_sets_axes(aes(fill = y),
                          axis.width = 0.2,
                          sep        = 0.05) +
  geom_parallel_sets_labels(angle = 0, size = 3) +
  # apply your unified palette
  scale_fill_manual(values = fill_map) +
  scale_x_discrete(limits = c("monkey","human","mouse"),
                   labels = c("Monkey","Human","Mouse")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x   = element_blank(),
        axis.text.x    = element_text(size = 10, face = "bold"))

ggsave("celltype2_knn.pred.stag-seperate1.pdf",plot=p,width=5,height=5)


********************************************
*****************cellchat_plot**************
********************************************
# Packages
library(dplyr)
library(purrr)
library(ggplot2)

cellchat <- readRDS("cellchat.rds")

cell_levels <- c(
  "Notochord","NMP-Neural","NMP-Meso","Early NT","Early FP","Neural Tube","Floor Plate",
  "pPSM","aPSM","Sclerotome","Syndetome","Endotome","Myogenic Progenitors","Dermomyotome",
  "Intermediate Meso","Endothelial Cells"
)
cellchat@idents <- factor(cellchat@idents, levels = cell_levels)

pairLR.use <- extractEnrichedLR(
  cellchat,
  signaling = c("BMP","TGFb","WNT","ncWNT","FGF","NOTCH","HH","FLRT","PCDH","RA",
                "PDGF","IGF","CLDN","VEGF","Desmosterol","IGFBP")
)

# Helper
subset_comm <- function(src, tgt, thr = 1) {
  subsetCommunication(
    object = cellchat, slot.name = "net",
    sources.use = src, targets.use = tgt,
    pairLR.use = pairLR.use, thresh = thr
  )
}

# Define all source–target requests once
requests <- list(
  list("Notochord",     c("NMP-Neural","NMP-Meso","Neural Tube","Floor Plate","Syndetome","Sclerotome"), 1),
  list("NMP-Meso",      c("aPSM","Neural Tube","Myogenic Progenitors","Dermomyotome"),                    1),
  list("NMP-Neural",    c("Neural Tube","Sclerotome","Syndetome","Myogenic Progenitors"),                 1),
  list(c("Floor Plate","Syndetome","Dermomyotome"), "Neural Tube",                                        1),
  list("Neural Tube",   c("Floor Plate","Syndetome","Dermomyotome"),                                      1),
  list("Floor Plate",   c("Sclerotome","Syndetome","Myogenic Progenitors"),                               1),
  list("Sclerotome",    c("Syndetome","Myogenic Progenitors"),                                            0.05),
  list("Syndetome",     c("Endotome","Myogenic Progenitors"),                                             0.05),
  list("Endotome",      c("Myogenic Progenitors","Endothelial Cells"),                                    0.05)
)

signaling <- purrr::pmap_dfr(requests, subset_comm)


genes <- unique(unlist(strsplit(signaling$interaction_name, "_")))

combined.sct <- readRDS("01_NTSM_integrated_decontX_up3.rds")
Idents(combined.sct) <- "detailed_celltype_2nd"

# AverageExpression returns a list; take the requested assay matrix
avg <- AverageExpression(
  combined.sct, assays = "decontX", slot = "data", features = genes
)$decontX  # genes x celltypes

# Split interaction into ligand + receptors
parts      <- strsplit(signaling$interaction_name, "_")
ligands    <- vapply(parts, `[`, character(1), 1)
receptors  <- lapply(parts, function(x) if (length(x) > 1) x[-1] else character(0))

# Sum ligand (in source) + receptors (in target), handling missing gracefully
sum_lig <- mapply(function(g, cl) {
  if (g %in% rownames(avg) && cl %in% colnames(avg)) avg[g, cl] else 0
}, ligands, signaling$source)

sum_rec <- mapply(function(gs, cl) {
  gs <- gs[gs %in% rownames(avg)]
  if (length(gs) == 0 || !(cl %in% colnames(avg))) return(0)
  sum(avg[gs, cl, drop = FALSE])
}, receptors, signaling$target)

data <- signaling %>%
  mutate(
    sum_expression   = sum_lig + sum_rec,
    ID               = paste(source, target, sep = " - "),
    ID               = factor(ID, levels = unique(ID)),
    interaction_name = factor(interaction_name, levels = rev(sort(unique(interaction_name))))
  )


p <- ggplot(data, aes(x = ID, y = interaction_name,
                      size = -log2(prob),
                      color = log2(sum_expression))) +
  geom_point() +
  scale_color_gradientn(
    colors = c("#B6CAFA", "#FAA31C", "#FF1777", "#A700C0"),
    values = scales::rescale(c(-3, -1, 0, 2)),
    name   = "log2 Sum (Ligand, Receptor)"
  ) +
  scale_size_continuous(range = c(2, 5), name = "-log2(prob)") +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10, color = "black"),
    axis.text.y      = element_text(size = 8, color = "black"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey90"),
    panel.border     = element_rect(color = "grey30", fill = NA, size = 1),
    axis.ticks.length = unit(0.2, "cm")
  )

ggsave("dotplot_cellchat_selected.pdf", plot = p, width = 10, height = 10)
