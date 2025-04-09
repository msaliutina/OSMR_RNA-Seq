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

load("data/sc_raw/SO_reclustered_GSEA_celltypeAnnotation_Stroma.rdata")
Stroma <- SO
rm(SO)

Stroma[['percent.mt']] <- PercentageFeatureSet(Stroma, pattern='^mt-')

Stroma <- subset(Stroma, subset = percent.mt < 25)
Stroma <- NormalizeData(Stroma)


Idents(Stroma) <- Stroma$RNA_snn_res.0.1
Stroma <- RenameIdents(Stroma, `0` = 'Endothelial cells 1',
`1` = 'Stromal cells 1.1',
`2` = 'Stromal cells 1.2',
`3` = 'Lymphatic endothelial cells',
`4` = 'Oligodendrocytes',
`5` = 'Endothelial cells 2',
`6` = 'Smooth muscle cells',
`7` = 'Stromal cells 2',
`8` = 'Proliferating cells',
`9` = 'Smooth muscle cells',
`10` = 'Epithelial cells',
`11` = 'Stromal cells 3',
`12` = 'Epithelial cells',
`13` = 'Macrophages',
`14` = 'Epithelial cells',
`15` = 'B cells',
`16` = 'B cells')
Stroma@meta.data$celltypes <- Idents(Stroma)
Stroma <- subset(Stroma, idents = c('B cells', 'Macrophages', 'Epithelial cells'), invert = TRUE)
Stroma <- FindNeighbors(Stroma, dims = 1:30, verbose = F)
Stroma <- FindClusters(Stroma, resolution = 1.2, verbose = F)
Stroma <- RunUMAP(Stroma, dims = 1:30, verbose = F)
Idents(Stroma) <- Stroma@meta.data$celltypes

saveRDS(Stroma, 'data/sc_final/Stroma_final.rds')




