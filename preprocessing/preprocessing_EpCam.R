library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(RColorBrewer)
library(ggalluvial)
library(harmony)



EpCam <- readRDS('data/sc_raw/EpCam.rds')

goblet_4 <- WhichCells(EpCam, idents = 'Enterocytes_4')
EpCam <- SetIdent(EpCam, cells = goblet_4, value = 'Goblet')

EpCam@meta.data$celltypes_final[EpCam$celltypes_final == 'Enterocytes_4'] <- 'Goblet'
EpCam[['percent.mt']] <- PercentageFeatureSet(EpCam, pattern='^mt-')

EpCam <- subset(EpCam, subset = percent.mt < 25)
EpCam <- NormalizeData(EpCam)


harmonized_EpCam <- RunHarmony(EpCam, 
                               group.by.vars = c("experiment"), 
                               reduction = "pca", assay.use = "integrated", reduction.save = "harmony")

harmonized_EpCam <- RunUMAP(harmonized_EpCam, reduction = "harmony", assay = "integrated", dims = 1:20)

harmonized_EpCam <- FindNeighbors(object = harmonized_EpCam, reduction = "harmony")

harmonized_EpCam$RNA_snn_res.0.1 <- NULL
harmonized_EpCam$RNA_snn_res.0.2 <- NULL
harmonized_EpCam$RNA_snn_res.0.3 <- NULL
harmonized_EpCam$RNA_snn_res.0.4 <- NULL
harmonized_EpCam$RNA_snn_res.0.5 <- NULL
harmonized_EpCam$RNA_snn_res.0.6 <- NULL
harmonized_EpCam$RNA_snn_res.0.7 <- NULL
harmonized_EpCam$RNA_snn_res.0.8 <- NULL
harmonized_EpCam$RNA_snn_res.0.9 <- NULL

harmonized_EpCam <- FindClusters(harmonized_EpCam, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2))

Idents(harmonized_EpCam) <- harmonized_EpCam$RNA_snn_res.0.2

# Function to find markers for each cluster and store results in variables
find_cluster_markers <- function(seurat_obj) {
  cluster_markers <- list()  # Initialize a list to store markers for each cluster
  
  for (i in 0:9) {  # Iterate over cluster numbers 0 to 9
    markers <- FindMarkers(seurat_obj, ident.1 = i, min.pct = 0.25)
    cluster_markers[[paste0("cluster_", i)]] <- markers  # Store markers for the cluster
  }
  
  return(cluster_markers)
}

# Apply the function to your Seurat object
cluster_markers <- find_cluster_markers(harmonized_EpCam)

# Access markers for each cluster using variable names cluster_0, cluster_1, ..., cluster_9
cluster_0 <- cluster_markers[["cluster_0"]]
cluster_1 <- cluster_markers[["cluster_1"]]
cluster_2 <- cluster_markers[["cluster_2"]]
cluster_3 <- cluster_markers[["cluster_3"]]
cluster_4 <- cluster_markers[["cluster_4"]]
cluster_5 <- cluster_markers[["cluster_5"]]
cluster_6 <- cluster_markers[["cluster_6"]]
cluster_7 <- cluster_markers[["cluster_7"]]
cluster_8 <- cluster_markers[["cluster_8"]]
cluster_9 <- cluster_markers[["cluster_9"]]


harmonized_EpCam <- RenameIdents(harmonized_EpCam,
                                 `0` = 'Enterocytes_1',
                                 `1` = 'Enterocytes_2',
                                 `2` = 'Enterocytes_3',
                                 `3` = 'Enterocytes_4',
                                 `4` = 'Proliferating epithelial cells',
                                 `5` = 'TA',
                                 `6` = 'Goblet',
                                 `7` = 'Goblet',
                                 `8` = 'Tuft',
                                 `9` = 'Enteroendocrine')

harmonized_EpCam$celltypes <- Idents(harmonized_EpCam)

harmonized_EpCam <- FindSubCluster(harmonized_EpCam, cluster = 'Enterocytes_1', resolution = 0.2, graph.name = 'RNA_snn')

harmonized_EpCam <- SetIdent(harmonized_EpCam, value = 'sub.cluster')

enterocytes_1_0 <- WhichCells(harmonized_EpCam, idents = 'Enterocytes_1_0')
enterocytes_1_1 <- WhichCells(harmonized_EpCam, idents = c('Enterocytes_1_1'))
harmonized_EpCam <- SetIdent(harmonized_EpCam, cells = enterocytes_1_0, value = 'Enterocytes_1_0')
harmonized_EpCam <- SetIdent(harmonized_EpCam, cells = enterocytes_1_1, value = 'Enterocytes_1_1')
harmonized_EpCam$celltypes <- Idents(harmonized_EpCam)


saveRDS(harmonized_EpCam, "data/sc_final/EpCam_harmony.rds")