*********************************************************************
*****************************V2_Spatial******************************
*********************************************************************

*****************************************************************
*****************************01.V2_QC****************************
*****************************************************************


library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(arrow)
library(ggspatial)

#### 39L ####
filepath <- "./Aln_2/39L/outs"
srt <- Load10X_Spatial(data.dir=filepath, bin.size = c(8))

srt <- NormalizeData(srt)

srt <- FindVariableFeatures(srt)
srt <- ScaleData(srt)
srt <- RunPCA(srt)
srt <- FindNeighbors(srt, dims = 1:50)
srt <- FindClusters(srt, resolution = 1.6)
srt <- RunUMAP(srt, return.model = T, dims = 1:50)

saveRDS(srt, "v2_39L_filtered.rds")

#### 65L ####
filepath <- ".Aln_2/65T/outs"
srt <- Load10X_Spatial(data.dir=filepath, bin.size = c(8))

srt <- NormalizeData(srt)

srt <- FindVariableFeatures(srt)
srt <- ScaleData(srt)
srt <- RunPCA(srt)
srt <- FindNeighbors(srt, dims = 1:50)
srt <- FindClusters(srt, resolution = 1.8)
srt <- RunUMAP(srt, return.model = T, dims = 1:50)

saveRDS(srt, "v2_65T_filtered.rds")

*****************************************************************
*****************************02.V2_RCTD**************************
*****************************************************************

library(spacexr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(arrow)

#### 39L ####
# import object
visium.obj <- readRDS('v2_39L_filtered.rds')

query.counts <- GetAssayData(visium.obj, assay = "Spatial.008um", slot = "counts")[, Cells(visium.obj[['slice1.008um']])]
coords <- GetTissueCoordinates(visium.obj[['slice1.008um']], which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))

# import reference object
ref.obj <- readRDS("./V2_annotated.rds")
ref.obj <- UpdateSeuratObject(ref.obj)

Idents(ref.obj) <- "celltype"
# remove CR cells because there aren't enough of them for annotation
ref.obj <- subset(ref.obj, subset = celltype != "Neuron")
counts <- GetAssayData(ref.obj, assay = "decontX", slot = "counts")
cluster <- as.factor(ref.obj$celltype)
names(cluster) <- colnames(ref.obj)
nUMI <- ref.obj$nCount_decontX
names(nUMI) <- colnames(ref.obj)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
reference <- Reference(counts, cluster, nUMI)

# run RCTD with many cores
RCTD <- create.RCTD(query, reference, max_cores = 8)
RCTD <- run.RCTD(RCTD, doublet_mode = "full")

max_indices <- apply(RCTD@results$weights, 1, which.max)
celltypes <- colnames(RCTD@results$weights)
annotation_df = data.frame(cellid = rownames(RCTD@results$weights), celltype = celltypes[max_indices])
rownames(annotation_df) = annotation_df$cellid

visium.obj@meta.data$predicted.celltype <- annotation_df[rownames(visium.obj@meta.data), 'celltype']
keep.cells <- Cells(visium.obj)[!is.na(visium.obj$predicted.celltype)]
visium.obj <- subset(visium.obj, cells = keep.cells)

saveRDS(visium.obj, "v2_39L_rctd.rds")


#### 65T ####
# import object
visium.obj <- readRDS('v2_65T_filtered.rds')

query.counts <- GetAssayData(visium.obj, assay = "Spatial.008um", slot = "counts")[, Cells(visium.obj[['slice1.008um']])]
coords <- GetTissueCoordinates(visium.obj[['slice1.008um']], which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))

# import reference object
ref.obj <- readRDS("./V2_annotated.rds")
ref.obj <- UpdateSeuratObject(ref.obj)

Idents(ref.obj) <- "celltype"
# remove CR cells because there aren't enough of them for annotation
ref.obj <- subset(ref.obj, subset = celltype != "Neuron")
counts <- GetAssayData(ref.obj, assay = "decontX", slot = "counts")
cluster <- as.factor(ref.obj$celltype)
names(cluster) <- colnames(ref.obj)
nUMI <- ref.obj$nCount_decontX
names(nUMI) <- colnames(ref.obj)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
reference <- Reference(counts, cluster, nUMI)

# run RCTD with many cores
RCTD <- create.RCTD(query, reference, max_cores = 8)
RCTD <- run.RCTD(RCTD, doublet_mode = "full")

max_indices <- apply(RCTD@results$weights, 1, which.max)
celltypes <- colnames(RCTD@results$weights)
annotation_df = data.frame(cellid = rownames(RCTD@results$weights), celltype = celltypes[max_indices])
rownames(annotation_df) = annotation_df$cellid

visium.obj@meta.data$predicted.celltype <- annotation_df[rownames(visium.obj@meta.data), 'celltype']
keep.cells <- Cells(visium.obj)[!is.na(visium.obj$predicted.celltype)]
visium.obj <- subset(visium.obj, cells = keep.cells)

saveRDS(visium.obj, "v2_65T_rctd.rds")

*****************************************************************
*************************03.V2_annotation************************
*****************************************************************

# r-seurat5

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(arrow)
library(future)
library(svglite)

# Annotate 39L
srt <- readRDS('v2_39L_rctd.rds')

srt <- FindClusters(srt, resolution = 0.8)

cluster_to_celltype <- c(
  "10" = "Tail Bud",
  "7" = "Neural Lineage",
  "2" = "Neural Lineage",
  "14" = "Neural Lineage",
  "9" = "Neural Lineage",
  "13" = "Neural Lineage",
  "11" = "Somitic Lineage",
  "4" = "Somitic Lineage",
  "5" = "Somitic Lineage",
  "1" = "Somitic Lineage",
  "0" = "Somitic Lineage",
  "6" = "Somitic Lineage",
  "8" = "Somitic Lineage",
  "12" = "Endothelial Cells",
  "3" = "Neural Lineage"
)

cluster_to_detailed_celltype <- c(
  "10" = "NMP-Neural / Tail bud progenitors",
  "7" = "Caud. NP",
  "2" = "NT",
  "14" = "NT",
  "9" = "FP-like",
  "13" = "Neuron",
  "11" = "PSM",
  "4" = "Somite",
  "5" = "Sclerotome",
  "1" = "Sclerotome",
  "0" = "Sclerotome",
  "6" = "Syndetome",
  "8" = "Endotome",
  "12" = "Endothelial Cells",
  "3" = "Neuron"
)

srt@meta.data['celltype'] <- cluster_to_celltype[as.character(srt@meta.data$seurat_clusters)]
srt@meta.data['detailed_celltype'] <- cluster_to_detailed_celltype[as.character(srt@meta.data$seurat_clusters)]

# Annotate 65T
srt <- readRDS('v2_65T_rctd.rds')

srt <- FindClusters(srt, resolution = 1.4)

cluster_to_celltype <- c(
  "2" = "Tail bud",
  "10" = "Tail bud",
  "4" = "Tail bud",
  "9" = "Tail bud",
  "7" = "Neural Lineage",
  "5" = "Neural Lineage",
  "8" = "Neural Lineage",
  "3" = "Neural Lineage",
  "6" = "Somitic Lineage",
  "0" = "Somitic Lineage",
  "1" = "Somitic Lineage",
  "12" = "EC",
  "11" = "Mixed Meso"
)

cluster_to_detailed_celltype <- c(
  "2" = "Tail bud progenitors",
  "10" = "Tail bud progenitors",
  "4" = "NMP-Neural",
  "9" = "Node-like",
  "7" = "NT",
  "5" = "NT",
  "8" = "Doral NT",
  "3" = "Ventral NT",
  "6" = "Dorsal SM",
  "0" = "Ventral SM",
  "1" = "Somitic Meso",
  "12" = "Endothelial Cells",
  "11" = "Mixed Meso"
)

srt@meta.data['celltype'] <- cluster_to_celltype[as.character(srt@meta.data$seurat_clusters)]
srt@meta.data['detailed_celltype'] <- cluster_to_detailed_celltype[as.character(srt@meta.data$seurat_clusters)]


*********************************************************************
*****************************V3_Spatial******************************
*********************************************************************

*****************************************************************
*****************************01.V3_QC****************************
*****************************************************************
# Seurat 5

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(arrow)
library(future)


#### Preprocess v3 Day4 ####

filepath <- "./hTLOv3_Day4_slide9/outs"
srt <- Load10X_Spatial(data.dir=filepath, bin.size = c(8))
srt@meta.data['orig.ident'] <- 'v3_Day4'

# Filter low quality cells
srt <- subset(srt, subset = nCount_Spatial.008um > 700)

# Normalize
srt <- NormalizeData(srt)

# Comupte CC genes
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
srt <- CellCycleScoring(srt, s.features = s.genes, g2m.features = g2m.genes)

# Dimension reduction and clustering
DefaultAssay(srt) <- "Spatial.008um"
srt <- FindVariableFeatures(srt)
srt <- ScaleData(srt)
srt <- RunPCA(srt)
srt <- FindNeighbors(srt, dims = 1:50)
srt <- FindClusters(srt, resolution = 1.6)
srt <- RunUMAP(srt, return.model = T, dims = 1:50)

# Save object
saveRDS(srt, 'v3_day4_spatial.rds')



#### Preprocess v3 Day6.5 ####

filepath <- "./hTLOv3_Day6_5_slide22/outs"
srt <- Load10X_Spatial(data.dir=filepath, bin.size = c(8))
srt@meta.data['orig.ident'] <- 'v3_Day4'

# Filter low quality cells
srt.subset <- subset(srt, subset = nCount_Spatial.008um < 200 & percent.mt < 5)

# Normalize
srt <- NormalizeData(srt)

# Comupte CC genes
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
srt <- CellCycleScoring(srt, s.features = s.genes, g2m.features = g2m.genes)

# Dimension reduction and clustering
DefaultAssay(srt) <- "Spatial.008um"
srt <- FindVariableFeatures(srt)
srt <- ScaleData(srt)
srt <- RunPCA(srt)
srt <- FindNeighbors(srt, dims = 1:50)
srt <- FindClusters(srt, resolution = 1.6)
srt <- RunUMAP(srt, return.model = T, dims = 1:50)

# Save object
saveRDS(srt, 'v3_day6.5_spatial.rds')


*****************************************************************
*****************************02.V3_RCTD**************************
*****************************************************************

library(spacexr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(arrow)
# library(future)

# options(future.globals.maxSize= 16*1024^3)

#### Day 4 spatial ####
# Read Day 4 spatial
visium.obj <- readRDS('v3_day4_spatial.rds')

query.counts <- GetAssayData(visium.obj, assay = "Spatial.008um", slot = "counts")[, Cells(visium.obj[['slice1.008um']])]
coords <- GetTissueCoordinates(visium.obj[['slice1.008um']], which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))

# Read reference
ref.obj <- readRDS("finalized_hTLO/01_NTSM_integrated_decontX_annotated_third_rctd.rds")
ref.obj <- UpdateSeuratObject(ref.obj)

Idents(ref.obj) <- "rctd_selected_cell_type"
counts <- GetAssayData(ref.obj, assay = "decontX", slot = "counts")
cluster <- as.factor(ref.obj$rctd_selected_cell_type)
names(cluster) <- colnames(ref.obj)
nUMI <- ref.obj$nCount_decontX
names(nUMI) <- colnames(ref.obj)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
reference <- Reference(counts, cluster, nUMI)

# run RCTD with many cores
RCTD <- create.RCTD(query, reference, max_cores = 8)
RCTD <- run.RCTD(RCTD, doublet_mode = "full")

max_indices <- apply(RCTD@results$weights, 1, which.max)
celltypes <- colnames(RCTD@results$weights)
annotation_df = data.frame(cellid = rownames(RCTD@results$weights), celltype = celltypes[max_indices])
rownames(annotation_df) = annotation_df$cellid

visium.obj@meta.data$predicted.rctd_selected_cell_type <- annotation_df[rownames(visium.obj@meta.data), 'celltype']
keep.cells <- Cells(visium.obj)[!is.na(visium.obj$predicted.rctd_selected_cell_type)]
visium.obj <- subset(visium.obj, cells = keep.cells)

# Save object
saveRDS(visium.obj, "v3_day4_spatial_rctd.rds")


#### Day 6.5 spatial ####
# Read Day 6.5 spatial
visium.obj <- readRDS('v3_day6.5_spatial.rds')

query.counts <- GetAssayData(visium.obj, assay = "Spatial.008um", slot = "counts")[, Cells(visium.obj[['slice1.008um']])]
coords <- GetTissueCoordinates(visium.obj[['slice1.008um']], which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))

# Read reference
ref.obj <- readRDS("finalized_hTLO/01_NTSM_integrated_decontX_annotated_third_rctd.rds")
ref.obj <- UpdateSeuratObject(ref.obj)

Idents(ref.obj) <- "rctd_selected_cell_type"
# remove CR cells because there aren't enough of them for annotation
# ref.obj <- subset(ref.obj, subset = subclass != "CR")
counts <- GetAssayData(ref.obj, assay = "decontX", slot = "counts")
cluster <- as.factor(ref.obj$rctd_selected_cell_type)
names(cluster) <- colnames(ref.obj)
nUMI <- ref.obj$nCount_decontX
names(nUMI) <- colnames(ref.obj)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
reference <- Reference(counts, cluster, nUMI)

# run RCTD with many cores
RCTD <- create.RCTD(query, reference, max_cores = 8)
RCTD <- run.RCTD(RCTD, doublet_mode = "full")

max_indices <- apply(RCTD@results$weights, 1, which.max)
celltypes <- colnames(RCTD@results$weights)
annotation_df = data.frame(cellid = rownames(RCTD@results$weights), celltype = celltypes[max_indices])
rownames(annotation_df) = annotation_df$cellid

visium.obj@meta.data$predicted.rctd_selected_cell_type <- annotation_df[rownames(visium.obj@meta.data), 'celltype']
keep.cells <- Cells(visium.obj)[!is.na(visium.obj$predicted.rctd_selected_cell_type)]
visium.obj <- subset(visium.obj, cells = keep.cells)

# Save object
saveRDS(visium.obj, "v3_day6.5_spatial_rctd.rds")


*****************************************************************
*************************03.V3_annotation************************
*****************************************************************

# Seurat 5

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(arrow)
library(future)
library(scCustomize)
library(ComplexHeatmap)

#### Day 4 ####
srt <- readRDS("v3_day4_spatial_rctd.rds")

srt <- FindClusters(srt, resolution = 0.75, verbose=F)

cluster_annotations_1 <- c(
  "0" = "Intermediate Meso",
  "1" = "SM.1",
  "2" = "SM.2",
  "3" = "SM.3",
  "6" = "SM.4",
  "7" = "SM.5",
  "5" = "Neural Plate",
  "8" = "NMP-Meso / pPSM",
  "9" = "aPSM",
  "10" = "NMP-Neural",
  "11" = "Notochord",
  "12" = "Tail Bud Progenitors",
  "4" = "Ambiguous"
)

srt@meta.data$spatial_annot <- cluster_annotations_1[as.character(srt@meta.data$seurat_clusters)]

saveRDS(srt, "v3_day4_spatial_rctd_annotated.rds")

#### Day 6.5 ####
srt <- readRDS("v3_day6.5_spatial_rctd.rds")

srt <- FindClusters(srt, resolution = 1.0, verbose=F)

cluster_annotations <- c(
  "14" = "Notochord",
  "16" = "aPSM",
  "3" = "SM.1",
  "1" = "SM.2",
  "2" = "SM.3",
  "4" = "SM.4",
  "0" = "SM.5",
  "10" = "SM.6",
  "17" = "NMP-Neural",
  "8" = "Caud. NP",
  "15" = "Pre-NT",
  "6" = "Neural Tube.1",
  "13" = "Neural Tube.2",
  "9" = "Floor Plate",
  "12" = "Mixed Caud. Progenitors",
  "18" = "Int Meso",
  "5" = "Endothelial Cells",
  "19" = "Neuroectoderm",
  "7" = "Ambiguous",
  "11" = "Tail Bud Progenitors"
)

srt@meta.data$spatial_annot <- cluster_annotations[as.character(srt@meta.data$seurat_clusters)]

saveRDS(srt, "v3_day6.5_spatial_rctd_annotated.rds")

*************************************************
******************spatial plot*******************
*************************************************


********01.ap_axis_plot

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(arrow)
library(future)
library(svglite)
library(dbscan)

rotate_90_clockwise <- function(x, y) {
  x_new <- y
  y_new <- -x
  return(data.frame(x = x_new, y = y_new))
}

get_perpendicular_lines <- function(x1, x2, y1, y2, spacing = 100, length=3000){
    # Desired spacing and perpendicular line length
    line_length <- length

    # Calculate direction vector
    dx <- x2 - x1
    dy <- y2 - y1

    # Calculate length of the line segment
    length <- sqrt(dx^2 + dy^2)

    # Normalize the direction vector
    dx_unit <- dx / length
    dy_unit <- dy / length

    # Perpendicular vectors
    perp_dx <- -dy_unit
    perp_dy <- dx_unit

    # Generate points along the line at regular intervals
    t_values <- seq(0, length, by = spacing)
    perpendicular_lines <- data.frame()

    for (t in t_values) {
    # Position on the main line
    x_mid <- x1 + t * dx_unit
    y_mid <- y1 + t * dy_unit

    # Perpendicular line endpoints
    x1_perp <- x_mid + (line_length / 2) * perp_dx
    y1_perp <- y_mid + (line_length / 2) * perp_dy
    x2_perp <- x_mid - (line_length / 2) * perp_dx
    y2_perp <- y_mid - (line_length / 2) * perp_dy

    # Store the line segment
    perpendicular_lines <- rbind(perpendicular_lines, 
                                data.frame(x1 = x1_perp, y1 = y1_perp, 
                                            x2 = x2_perp, y2 = y2_perp))
    }
    return(perpendicular_lines)
}

# Function to calculate distance from a point to a line segment
point_to_line_distance <- function(px, py, x1, y1, x2, y2) {
  # Line direction vector
  dx <- x2 - x1
  dy <- y2 - y1
  
  # Handle vertical/horizontal lines
  if (dx == 0 && dy == 0) {
    # Line is a single point
    return(sqrt((px - x1)^2 + (py - y1)^2))
  }
  
  # Parameter t for projection onto the line
  t <- ((px - x1) * dx + (py - y1) * dy) / (dx^2 + dy^2)
  
  # Clamp t to [0,1] to stay within the segment
  t <- max(0, min(1, t))
  
  # Closest point on the line segment
  closest_x <- x1 + t * dx
  closest_y <- y1 + t * dy
  
  # Distance from point to closest point on the line
  return(sqrt((px - closest_x)^2 + (py - closest_y)^2))
}

get_closest_points <- function(point_df, perpendicular_lines, distance_threshold=100){
  # Initialize a dataframe to store results
  results <- data.frame(point_id = integer(), line_id = integer(), distance = numeric())

  # Loop through points and lines
  for (i in 1:nrow(point_df)) {
    px <- point_df$x[i]
    py <- point_df$y[i]
    point_id <- point_df$id[i]
    
    for (j in 1:nrow(perpendicular_lines)) {
      x1 <- perpendicular_lines$x1[j]
      y1 <- perpendicular_lines$y1[j]
      x2 <- perpendicular_lines$x2[j]
      y2 <- perpendicular_lines$y2[j]
      line_id <- perpendicular_lines$line_id[j]
      
      # Calculate distance from point to line
      dist <- point_to_line_distance(px, py, x1, y1, x2, y2)
      
      # If within the threshold, store the result
      if (dist <= distance_threshold) {
        results <- rbind(results, data.frame(point_id = point_id, line_id = line_id, distance = dist))
      }
    }
  }

  # Find the closest line for each point
  closest_lines <- results %>%
    group_by(point_id) %>%
    slice_min(order_by = distance, n = 1) %>%
    ungroup()

  closest_lines %>% data.frame() -> closest_lines
  # rownames(closest_lines) <- closest_lines$point_id
  return(closest_lines)

}

library(sf)
library(dplyr)

# Function to find and update intersecting lines
find_and_update_intersections <- function(data) {
  # Convert the data frame to an sf object with LINESTRING geometry
  lines_sf <- data %>%
    rowwise() %>%
    mutate(
      geometry = st_sfc(st_linestring(matrix(c(x1, y1, x2, y2), ncol = 2, byrow = TRUE)))
    ) %>%
    st_as_sf(crs = "+proj=longlat +datum=WGS84 +no_defs") # Set CRS explicitly
  
  # Find intersections
  intersections <- st_intersection(lines_sf)
  
  # Check if there are any intersections
  if (nrow(intersections) == 0) {
    message("No intersections found.")
    return(lines_sf)
  }
  
  # Extract intersection points (POINT geometries)
  intersection_points <- intersections %>%
    filter(st_geometry_type(geometry) == "POINT") %>%
    st_coordinates() %>%
    as.data.frame() %>%
    rename(x = X, y = Y)
  
  # Update the endpoints of intersecting lines to end at the intersection
  updated_lines <- lines_sf %>%
    rowwise() %>%
    mutate(
      geometry = {
        # Convert current line to geometry
        current_line <- st_cast(geometry, "LINESTRING")
        
        # Find intersection points for the current line
        intersecting_points <- st_intersection(current_line, st_as_sf(intersection_points, coords = c("x", "y"), crs = st_crs(lines_sf)))
        
        # If intersections exist, update the line's endpoint
        if (length(intersecting_points) > 0) {
          intersect_coords <- st_coordinates(intersecting_points)[1, ]
          st_sfc(st_linestring(matrix(c(x1, y1, intersect_coords[1], intersect_coords[2]), ncol = 2, byrow = TRUE)))
        } else {
          geometry # Keep original line if no intersection
        }
      }
    )
  
  return(updated_lines)
}

find_intersection <- function(x1, y1, x2, y2, x3, y3, x4, y4) {
  # Calculate the denominator of the intersection formula
  denominator <- (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
  
  # Check if the lines are parallel (denominator = 0)
  if (denominator == 0) {
    return(list(intersects = FALSE, point = NULL)) # No intersection
  }
  
  # Calculate the intersection point
  px <- ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) / denominator
  py <- ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) / denominator
  
  # Check if the intersection point lies on both line segments
  if (
    px >= min(x1, x2) && px <= max(x1, x2) &&
    py >= min(y1, y2) && py <= max(y1, y2) &&
    px >= min(x3, x4) && px <= max(x3, x4) &&
    py >= min(y3, y4) && py <= max(y3, y4)
  ) {
    return(list(intersects = TRUE, point = c(px, py))) # Intersection point
  } else {
    return(list(intersects = FALSE, point = NULL)) # No intersection on segments
  }
}

find_closer_point <- function(x1, y1, x2, y2, x3, y3) {
  # Calculate Euclidean distances
  dist1 <- sqrt((x1 - x3)^2 + (y1 - y3)^2) # Distance from (x1, y1) to (x3, y3)
  dist2 <- sqrt((x2 - x3)^2 + (y2 - y3)^2) # Distance from (x2, y2) to (x3, y3)
  
  # Compare distances and return 1 or 2
  if (dist1 < dist2) {
    return(1) # (x1, y1) is closer
  } else {
    return(2) # (x2, y2) is closer
  }
}

find_longer_line <- function(x1, y1, x2, y2, x3, y3) {
  # Calculate Euclidean distances (line lengths)
  length1 <- sqrt((x1 - x3)^2 + (y1 - y3)^2) # Length of line (x1, y1) to (x3, y3)
  length2 <- sqrt((x2 - x3)^2 + (y2 - y3)^2) # Length of line (x2, y2) to (x3, y3)
  
  # Compare lengths and return 1 or 2
  if (length1 > length2) {
    return(2) # Line (x1, y1) to (x3, y3) is longer
  } else {
    return(1) # Line (x2, y2) to (x3, y3) is longer
  }
}

find_middle_slope <- function(x1, y1, x2, y2, x3, y3, x4, y4) {
  # Calculate the slopes of the two lines
  slope1 <- ifelse((x2 - x1) != 0, (y2 - y1) / (x2 - x1), Inf) # Handle vertical line
  slope2 <- ifelse((x4 - x3) != 0, (y4 - y3) / (x4 - x3), Inf) # Handle vertical line
  
  # Handle case when one or both lines are vertical
  if (is.infinite(slope1) && is.infinite(slope2)) {
    stop("Both lines are vertical. No middle slope exists.")
  } else if (is.infinite(slope1)) {
    return(slope2) # Middle slope is the second slope
  } else if (is.infinite(slope2)) {
    return(slope1) # Middle slope is the first slope
  }
  
  # Calculate the middle slope (average of the two slopes)
  middle_slope <- (slope1 + slope2) / 2
  return(middle_slope)
}

find_endpoint <- function(x1, y1, slope, distance) {
  # Check if the slope is vertical
  if (is.infinite(slope)) {
    # Vertical line: x2 remains the same, y2 changes by distance
    x2 <- x1
    y2 <- y1 + distance # Positive distance moves up; negative distance moves down
  } else {
    # For non-vertical lines, calculate x2 and y2
    x2 <- x1 + distance / sqrt(1 + slope^2)
    y2 <- y1 + slope * (distance / sqrt(1 + slope^2))
  }
  
  return(c(x2 = x2, y2 = y2))
}


#### 39L ####

srt <- readRDS('v2_39L_rctd-1.rds')

coord <-  GetTissueCoordinates(srt[['slice1.008um']], which = "centroids")
coord <- rotate_90_clockwise(coord$x, coord$y)
rownames(point_df) <- rownames(coord)

# Use DBScan to identify individual organoids
result <- dbscan(point_df[c('x', 'y')], eps = 1000, minPts = 2)
point_df$cluster <- as.factor(result$cluster)

# Isolate organoid of interest
point_df <- point_df %>% filter(cluster == 1)
point_df['id'] <- rownames(point_df)

# Plot lines along the ap axis
spacing <- 200

x1 <- 17000
y1 <- -7450
x2 <- 15050
y2 <- -4950
length <- 4000
perpendicular_lines_1 <- get_perpendicular_lines(x1, x2, y1, y2, spacing, length)
perpendicular_lines_1['line_id'] <- 1:nrow(perpendicular_lines_1)

x1 <- 15050
y1 <- -4950
x2 <- 15250
y2 <- -3700
length <- 4000
perpendicular_lines_2 <- get_perpendicular_lines(x1, x2, y1, y2, spacing, length)
perpendicular_lines_2['line_id'] <- (max(perpendicular_lines_1['line_id'])+1):(max(perpendicular_lines_1['line_id'])+nrow(perpendicular_lines_2))

x1 <- 15250
y1 <- -3700
x2 <- 16500
y2 <- -1900
length <- 4000
perpendicular_lines_3 <- get_perpendicular_lines(x1, x2, y1, y2, spacing, length)
perpendicular_lines_3['line_id'] <- (max(perpendicular_lines_2['line_id'])+1):(max(perpendicular_lines_2['line_id'])+nrow(perpendicular_lines_3))

perpendicular_lines <- rbind(perpendicular_lines_1, perpendicular_lines_2, perpendicular_lines_3)

# update lines if there are any collisions
perpendicular_lines_updated <- perpendicular_lines
for (i in 1:nrow(perpendicular_lines_updated)){
    for (j in i:nrow(perpendicular_lines_updated)){
        if (i != j){
            x1 <- perpendicular_lines_updated$x1[i]
            y1 <- perpendicular_lines_updated$y1[i]
            x2 <- perpendicular_lines_updated$x2[i]
            y2 <- perpendicular_lines_updated$y2[i]
            x3 <- perpendicular_lines_updated$x1[j]
            y3 <- perpendicular_lines_updated$y1[j]
            x4 <- perpendicular_lines_updated$x2[j]
            y4 <- perpendicular_lines_updated$y2[j]

            intersection_result <- find_intersection(x1, y1, x2, y2, x3, y3, x4, y4)

            if (intersection_result$intersects == TRUE){
                pt <- intersection_result$point
                int_x1 <- pt[[1]]
                int_y1 <- pt[[2]]

                middle_slope <- find_middle_slope(x1, y1, x2, y2, x3, y3, x4, y4)
                end_point_1 <- find_endpoint(int_x1-0.1, int_y1-0.1, middle_slope, 500)
                end_point_2 <- find_endpoint(int_x1+0.1, int_y1+0.1, middle_slope, 500)

                for (line in c(i, j)){
                    x1 <- perpendicular_lines_updated$x1[line]
                    y1 <- perpendicular_lines_updated$y1[line]
                    x2 <- perpendicular_lines_updated$x2[line]
                    y2 <- perpendicular_lines_updated$y2[line]

                    closer <- find_longer_line(x1, y1, x2, y2, int_x1, int_y1)
                    if (closer == 1){
                        perpendicular_lines_updated[line, 'x1'] <- int_x1
                        perpendicular_lines_updated[line, 'y1'] <- int_y1
                    } else {
                        perpendicular_lines_updated[line, 'x2'] <- int_x1
                        perpendicular_lines_updated[line, 'y2'] <- int_y1
                    }

                    # perpendicular_lines_updated <- rbind(perpendicular_lines_updated, c(int_x1-0.1, int_y1-0.1, end_point_1[['x2']], end_point_1[['y2']], i))
                    # perpendicular_lines_updated <- rbind(perpendicular_lines_updated, c(int_x1+0.1, int_y1+0.1, end_point_2[['x2']], end_point_2[['y2']], j))
                }

            }

        }
    }
}

# Plot line segments
point_df %>% ggplot(aes(x = x, y = y)) + 
    geom_point() + 
    coord_fixed() +
    geom_segment(x = 17000, y = -7450, xend = 15050, yend = -4950, color='red') +
    geom_segment(x = 15050, y = -4950, xend = 15250, yend = -3700, color='red') +
    geom_segment(x = 15250, y = -3700, xend = 16500, yend = -1900, color='red') + 
    geom_segment(data = perpendicular_lines_updated, aes(x = x1, y = y1, xend = x2, yend = y2, color=as.factor(line_id)))

# Assign cells to their closest line
distance_threshold <- 500
closest_lines <- get_closest_points(point_df, perpendicular_lines_updated, distance_threshold)

max_per_group <- closest_lines %>%
  group_by(point_id) %>%
  summarize(line_id = max(line_id, na.rm = TRUE))
closest_lines <- max_per_group %>% data.frame()
rownames(closest_lines) <- closest_lines$point_id
point_df['closest_line'] <- closest_lines[rownames(point_df), 'line_id']

# Assign line identity to cells in the seurat object
srt@meta.data['band'] <- -1
srt@meta.data[rownames(point_df), 'band'] <- point_df$closest_line
srt@meta.data$band <- factor(srt@meta.data$band, levels = sort(unique(srt@meta.data$band)))

# Obtain average gene expression of target genes in each segment/band/line
target_genes <- c('TBX18', 'UNCX')

expr <- AggregateExpression(srt, features = target_genes, group.by = 'band')
expr$Spatial.008um %>% data.frame() %>% t() -> expr
rownames(expr) <- sort(unique(srt@meta.data$band))
expr <- expr %>% data.frame()
expr$band <- sort(unique(srt@meta.data$band))
expr$sample <- 's'
expr <- expr %>% filter(band != -1)

all_expr <- data.frame()
for (g in target_genes){
    expr_1 <- expr %>% select(band, sample)
    expr_1['expr_raw'] <- expr[,g]
    expr_1['expr_zscore'] <- scale(expr[,g])
    expr_1['expr_minmax'] <- expr[,g]/max(expr[,g])
    colnames(expr_1) <- c('band', 'sample', 'expr_raw', 'expr_zscore', 'expr_minmax')
    expr_1['gene'] <- g
    all_expr <- rbind(all_expr, expr_1)
}
all_expr$band <- factor(all_expr$band, levels = sort(unique(as.numeric(all_expr$band))))

# Plot expression along axis
pal <- c(
  "#FD1616",
  "#732E8F",
  "#FED426",
  "#00FB0D",
  "#FAA31C",
  "#4F474D"
)

all_expr %>% filter(band != -1) %>%
    ggplot(aes(x = band, y = expr_minmax, group = gene, color=gene)) + 
        geom_jitter(size = 0, alpha = 0) + 
        geom_smooth(method = "loess", se = FALSE, span = 0.2) + 
        labs(y = 'Min-max scaled expression', x = '', color='Gene', title='v2 39L - Min-max scaled expression') + 
        theme(
            # Remove all grid lines
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            
            axis.line.y = element_line(color = "black"),
            axis.line.x = element_blank(),
            
            # axis.title.x = element_blank(),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 12),
            plot.title = element_text(size = 12, hjust = 0.5),
            axis.ticks.x = element_blank()

        ) +
        scale_x_discrete(
            breaks = c('2', '34'),
            labels = c("P", "A")
        ) +
        scale_colour_manual(values = pal) -> p1

tiff('v2_39L_plot.tiff', width=7.5, height=4, units = 'in', res=300)
print(p1)
dev.off()


********02.spatial_two_gene_expression_plot

# r-seurat5

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(arrow)
library(future)
library(ComplexHeatmap)
library(scCustomize)
library(tidyr)
library(svglite)

# Double positive green: #02F321, blue: #0C76D3, yellow: #FFC30A. Negative #F4F4F4
joint_colour_gamma <- function(a, b,
                               gamma   = 0.75,
                               colA    = "#FF0000",   # blue
                               colB    = "#0000FF",   # red
                               colBoth = '#00FF00',   # will be filled in below
                               bg      = "#F4F4F4")
{
  # ---- helper --------------------------------------------------------------
  as_rgb <- function(col) col2rgb(col) / 255           # 0–1 matrix
  
  rgb_clip <- function(x) pmin(pmax(x, 0), 1)          # keep in [0,1]
  
  # ---- prepare expression values -------------------------------------------
  a <- rgb_clip(a)^gamma
  b <- rgb_clip(b)^gamma
  
  # ---- colours -------------------------------------------------------------
  rgbA <- as_rgb(colA)
  rgbB <- as_rgb(colB)
  
  if (is.null(colBoth)) {                              # additive default
    rgbBoth <- rgb_clip(rgbA + rgbB)
  } else {
    rgbBoth <- as_rgb(colBoth)
  }
  
  rgbBG <- as_rgb(bg)
  
  # ---- weights -------------------------------------------------------------
  w_exA   <- a * (1 - b)          # A only
  w_exB   <- b * (1 - a)          # B only
  w_both  <- pmin(a, b)           # overlap
  w_total <- pmax(a, b)           # overall intensity (0 … 1)
  
  # ---- build colour --------------------------------------------------------
  # pure colour (before blending with grey background)
  R <- w_exA * rgbA[1,] + w_exB * rgbB[1,] + w_both * rgbBoth[1,]
  G <- w_exA * rgbA[2,] + w_exB * rgbB[2,] + w_both * rgbBoth[2,]
  B <- w_exA * rgbA[3,] + w_exB * rgbB[3,] + w_both * rgbBoth[3,]
  
  # blend with background so that 0-expression → bg, full expression → pure
  R <- rgbBG[1,]*(1 - w_total) + R
  G <- rgbBG[2,]*(1 - w_total) + G
  B <- rgbBG[3,]*(1 - w_total) + B
  
  rgb(rgb_clip(R), rgb_clip(G), rgb_clip(B))
}

plot_two_genes_and_legend <- function(srt, gene1, gene2, pt.size=1, gamma=0.75, shape=21){
    expr <- FetchData(object = srt, vars = c(gene1, gene2))
    expr$normalized_g1 <- expr[[gene1]] / max(expr[[gene1]])
    expr$normalized_g2 <- expr[[gene2]] / max(expr[[gene2]])
    expr$final_color <- with(expr, joint_colour_gamma(normalized_g1, normalized_g2, gamma=gamma))

    srt@meta.data['merged_color'] <- expr[rownames(srt@meta.data), ]$final_color

    color_list <- c()
    for (i in unique(srt@meta.data$merged_color)){
        color_list[[i]] <- i
    }
    
    coord <-  GetTissueCoordinates(srt[['slice1.008um']], which = "centroids")
    myratio <- (max(coord$x) - min(coord$x)) / (max(coord$y) - min(coord$y))
    p1 <- SpatialDimPlot(srt, 'merged_color', cols = color_list, alpha = 0.8, image.alpha = 0.3, pt.size.factor = pt.size, shape=shape) +
        theme(legend.position = 'none', aspect.ratio = myratio, plot.title = element_text(hjust=0.5)) +
        labs(title = paste0(gene1,' (Red) and ', gene2, ' (Blue)'))

    scales_df <- data.frame(x = rep(seq(min(expr[[gene1]]), max(expr[[gene1]]), length.out = 11), 11),
        y = rep(seq(min(expr[[gene2]]), max(expr[[gene2]]), length.out = 11), each = 11),
        scale1 = rep(seq(0, 1, 0.1), 11),
        scale2 = rep(seq(0, 1, 0.1), each = 11))
    scales_df$final_color <- with(scales_df, joint_colour_gamma(scale1, scale2, gamma=gamma))

    p2 <- scales_df %>% ggplot(aes(x = x, y = y, fill = final_color)) +
        scale_fill_identity() +
        geom_tile() + 
        theme_minimal() +  # Use a minimal theme
        labs(x = gene1, y = gene2) +
        theme(
            axis.line = element_line(color = "black"),
            panel.grid = element_blank(),
            axis.line.x.top = element_blank(),
            axis.line.y.right = element_blank(),
            legend.position = 'none', aspect.ratio = 1
        ) + 
        coord_cartesian(xlim = c(min(scales_df$x), max(scales_df$x)), ylim = c(min(scales_df$y), max(scales_df$y)))
    
    plots <- c()
    plots[[1]] <- p1
    plots[[2]] <- p2

    return(plots)
}


##############################
# 39L
############################## 
srt <- readRDS("v2_39L_rctd.rds")
# 1. Define gene pairs
genes1 <- c("SOX2", "NOTO", "FOXA2", "EBF2", "CDH6", "FGF18", "EFNB2", "RIPPLY1", "WNT3A", "ALDH1A2", "HOXA1", "HOXC4", "CYP26A1", "CYP26A1")
genes2 <- c("SIX1", "FOXJ1", "SHH", "TGFBI", "CDH11", "FGFR3", "EPHA4", "LFNG", "FRZB", "RARB", "HOXA9", "HOXC9", "FGF17", "WNT3A")

# Combine the two vectors into a dataframe with specified column names
gene_pairs <- data.frame("gene1" = genes1, "gene2" = genes2, check.names = FALSE)

# Assert genes exists
existing_genes <- unique(colnames(FetchData(object = srt, vars = c(genes1, genes2))))

for (i in 1:nrow(gene_pairs)){
    gene1 <- gene_pairs$gene1[[i]]
    gene2 <- gene_pairs$gene2[[i]]
    if ((gene1 %in% existing_genes) && (gene2 %in% existing_genes)){
        plots <- plot_two_genes_and_legend(srt, gene1, gene2, gamma=0.5, shape=21, pt.size = 1.3)
        
        p1 <- plots[[1]]
        p2 <- plots[[2]]
        
        svglite(paste0("v2_39L_",i,"_",gene1,"_",gene2,".svg"),
            width  = 15,     # inches
            height = 15)
        print(p1)
        dev.off()

        svglite(paste0("v2_39L_",i,"_",gene1,"_",gene2,"_legend.svg"), width=6, height=6)
        print(p2)
        dev.off()
    }
}

##############################
# 65T
############################## 
srt <- readRDS('v2_65T_rctd.rds')
# 1. Define gene pairs
genes1 <- c("SOX2", "BMP4", "OLIG3", "WNT1", "WNT3A", "HES1", "ZIC2", "DLL3", "FOXA2", "EBF2", "CDH6", "FGF8", "FGF18", "ALDH1A2", "RDH10", "DHRS3")
genes2 <- c("SIX1", "NKX2-8", "OLIG1", "OLIG2", "FRZB", "HES5", "RFX4", "NOTCH1", "SHH", "SOX9", "CDH11", "FGFR1", "FGFR3", "RARB", "COL1A1", "PAX9")

# Combine the two vectors into a dataframe with specified column names
gene_pairs <- data.frame("gene1" = genes1, "gene2" = genes2, check.names = FALSE)

# Assert genes exists
existing_genes <- unique(colnames(FetchData(object = srt, vars = c(genes1, genes2))))

for (i in 1:nrow(gene_pairs)){
    gene1 <- gene_pairs$gene1[[i]]
    gene2 <- gene_pairs$gene2[[i]]
    if ((gene1 %in% existing_genes) && (gene2 %in% existing_genes)){
        plots <- plot_two_genes_and_legend(srt, gene1, gene2, gamma=0.5, shape=21, pt.size = 1.3)
        
        p1 <- plots[[1]]
        p2 <- plots[[2]]
        
        svglite(paste0("v2_65T_",i,"_",gene1,"_",gene2,".svg"),
            width  = 15,     # inches
            height = 15)
        print(p1)
        dev.off()

        svglite(paste0("v2_65T_",i,"_",gene1,"_",gene2,"_legend.svg"), width=6, height=6)
        print(p2)
        dev.off()
    }
}


********03.dotplot

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



