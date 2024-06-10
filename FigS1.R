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

EpCam <- readRDS('~/Desktop/Hegazy lab/single_cell_2023/Seurat objects/EpCam_harmony.rds')
LPMC <- readRDS('~/Desktop/Hegazy lab/single_cell_2023/Seurat objects/LPMC_final.Rds')
Stroma <- readRDS('~/Desktop/Hegazy lab/single_cell_2023/Seurat objects/Stroma_final.rds')


#FigS1A

Idents(LPMC) <- LPMC$Hashtag
Idents(EpCam) <- EpCam$Hashtag
Idents(Stroma) <- Stroma$Hashtag

Stroma_f <- Stroma[,!Stroma$Hashtag %in% 'inconclusive']
LPMC_f <- LPMC[,!LPMC$Hashtag %in% 'inconclusive']
EpCam_f <- EpCam[,!EpCam$Hashtag %in% 'inconclusive']


stromal_umap_hash <- DimPlot(Stroma_f,
                             reduction = "umap",
                             label = FALSE) 

epcam_umap_hash <- DimPlot(EpCam_f,
                           reduction = "umap",
                           label = FALSE) 


lpmc_umap_hash <- DimPlot(LPMC_f,
                          reduction = "umap",
                          label = FALSE) 

lpmc_umap_hash + epcam_umap_hash + stromal_umap_hash 

#FigS1B

Idents(LPMC) <- LPMC$experiment
Idents(EpCam) <- EpCam$experiment
Idents(Stroma) <- Stroma$experiment

stromal_umap_exp <- DimPlot(Stroma,
                            reduction = "umap",
                            label = FALSE) 

epcam_umap_exp <- DimPlot(EpCam,
                          reduction = "umap",
                          label = FALSE) 


lpmc_umap_exp <- DimPlot(LPMC,
                         reduction = "umap",
                         label = FALSE) 

lpmc_umap_exp + epcam_umap_exp + stromal_umap_exp

#FigS1C

LPMC_df_hashtag <- LPMC@meta.data %>% group_by(celltypes, experiment, Hashtag) %>% 
  summarise(cells = n( ))

LPMC_df_hashtag <- LPMC_df_hashtag %>% filter(Hashtag != 'inconclusive')

LPMC_df_hashtag$celltypes_exp <- paste(LPMC_df_hashtag$celltypes, LPMC_df_hashtag$experiment, sep='_')

LPMC_df_hashtag[nrow(LPMC_df_hashtag) + 1,] = list('Plasma cells', 'LPMCs_G2', 'Hashtag', 0, 'Plasma cells_LPMCs_G2')
LPMC_df_hashtag[nrow(LPMC_df_hashtag) + 1,] = list('Neutrophils', 'LPMCs_G3', 'Hashtag', 0, 'Neutrophils_LPMCs_G3')
LPMC_df_hashtag[nrow(LPMC_df_hashtag) + 1,] = list("Ilc3","LPMCs_G2", 'Hashtag', 0, 'Ilc3_LPMCs_G2')

lpmc_stats <- ggbarplot(LPMC_df_hashtag, x = "celltypes", y = "cells", add = "mean_se",
                        color = "experiment",
                        palette = c('red', 'black'),
                        position = position_dodge(-0.8))+
  scale_x_discrete(limits = c("Inflammatory Monocytes", "Neutrophils", "Bcells", "DCs", "Cd4",
                              "NK cells", "Ilc1", "Ilc2", "Ilc3", "Tregs", "Cd8", "Cd8 effector",
                              "Plasma cells", "Proliferating", "gd-tcells")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  stat_compare_means(aes(group = celltypes_exp), label = "p.format", label.y = 1500)


Epcam_df_hashtag <- EpCam@meta.data %>% group_by(celltypes, experiment, Hashtag) %>% 
  summarise(cells = n( ))

Epcam_df_hashtag <- Epcam_df_hashtag %>% filter(Hashtag != 'inconclusive')

Epcam_df_hashtag$celltypes_exp <- paste(Epcam_df_hashtag$celltypes, Epcam_df_hashtag$experiment, sep='_')

epcam_stats <- ggbarplot(Epcam_df_hashtag, x = "celltypes", y = "cells", add = "mean_se",
          color = "experiment",
          palette = c('red', 'black'),
          position = position_dodge(-0.8))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  scale_x_discrete(limits = c('Enterocytes_1_0',
                              'Enterocytes_1_1',
                              'Enterocytes_2',
                              'Enterocytes_3',
                              'Enterocytes_4',
                              'Goblet',
                              'TA',
                              'Enteroendocrine',
                              'Proliferating epithelial cells',
                              'Tuft')) +
  stat_compare_means(aes(group = celltypes_exp), label = "p.format", label.y = 2000)

Stroma_df_hashtag <- Stroma@meta.data %>% group_by(celltypes, experiment, Hashtag) %>% 
  summarise(cells = n( ))

Stroma_df_hashtag <- Stroma_df_hashtag %>% filter(Hashtag != 'inconclusive')

Stroma_df_hashtag$celltypes_exp <- paste(Stroma_df_hashtag$celltypes, Stroma_df_hashtag$experiment, sep='_')

Stroma_df_hashtag[nrow(Stroma_df_hashtag) + 1,] = list('Lymphatic endothelial cells', "Stroma_G3", "Hashtag", 0, 'Lymphatic endothelial cells_Stroma_G3')
Stroma_df_hashtag[nrow(Stroma_df_hashtag) + 1,] = list('Oligodendrocytes', "Stroma_G3", "Hashtag", 0, 'Oligodendrocytes_Stroma_G3')
Stroma_df_hashtag[nrow(Stroma_df_hashtag) + 1,] = list('Proliferating cells', "Stroma_G3", "Hashtag", 0, 'Proliferating cells_Stroma_G3')
Stroma_df_hashtag[nrow(Stroma_df_hashtag) + 1,] = list( 'Stromal cells 1.2', "Stroma_G3", "Hashtag", 0, 'Stromal cells 1.2_Stroma_G3')
Stroma_df_hashtag[nrow(Stroma_df_hashtag) + 1,] = list('Stromal cells 2', "Stroma_G3", "Hashtag", 0, 'Stromal cells 2_Stroma_G3')


stroma_stats <- ggbarplot(Stroma_df_hashtag, x = "celltypes", y = "cells", add = "mean_se",
                          color = "experiment",
                          palette = c('red', 'black'),
                          position = position_dodge(-0.8))+
  scale_x_discrete(limits = c('Endothelial cells 1',
                              'Stromal cells 1.1',
                              'Stromal cells 1.2',
                              'Lymphatic endothelial cells',
                              'Oligodendrocytes',
                              'Endothelial cells 2',
                              'Smooth muscle cells',
                              'Stromal cells 2',
                              'Proliferating cells',
                              'Mesothelial cells')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  stat_compare_means(aes(group = celltypes_exp), label = "p.format", label.y = 2500)

lpmc_stats + epcam_stats + stroma_stats

#FigS1D

Stroma@meta.data$celltypes_exp <- paste(Stroma@meta.data$celltypes, Stroma@meta.data$experiment, sep='_')

osm_stroma <- FetchData(Stroma, vars=c('Osm', 'Hashtag', 'experiment', 'celltypes', 'celltypes_exp'))

data.df.stroma <- osm_stroma %>% group_by(Hashtag, experiment, celltypes, celltypes_exp) %>% 
  summarise(value=mean(Osm))
colnames(data.df.stroma)[5] <- "Osm"

data.df.stroma <- data.df.stroma %>% filter(Hashtag != c('inconclusive'))

###add lacking rows
data.df.stroma[nrow(data.df.stroma) + 1,] = list("Hashtag","Stroma_G3", 'Lymphatic endothelial cells', 'Lymphatic endothelial cells_Stroma_G3', 0)
data.df.stroma[nrow(data.df.stroma) + 1,] = list("Hashtag","Stroma_G3", 'Oligodendrocytes', 'Oligodendrocytes_Stroma_G3', 0)
data.df.stroma[nrow(data.df.stroma) + 1,] = list("Hashtag","Stroma_G3", 'Proliferating cells', 'Proliferating cells_Stroma_G3', 0)
data.df.stroma[nrow(data.df.stroma) + 1,] = list("Hashtag","Stroma_G3", 'Stromal cells 1.2', 'Stromal cells 1.2_Stroma_G3', 0)
data.df.stroma[nrow(data.df.stroma) + 1,] = list("Hashtag","Stroma_G3", 'Stromal cells 2', 'Stromal cells 2_Stroma_G3', 0)


stroma_osm_expr <- ggbarplot(data.df.stroma, x = "celltypes", y = "Osm", add = "mean_se",
          color = "experiment", palette = c('red', 'black'), 
          position = position_dodge(-0.8))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  stat_compare_means(aes(group = celltypes_exp), label = "p.format", label.y = 2)

#FigS1E

osm_epcam <- FetchData(EpCam, vars=c('Osm', 'Hashtag', 'experiment', 'celltypes', 'celltypes_exp'))

data.df.epcam <- osm_epcam %>% group_by(Hashtag, experiment, celltypes, celltypes_exp) %>% 
  summarise(value=mean(Osm))
colnames(data.df.epcam)[5] <- "Osm"

data.df.epcam <- data.df.epcam %>% filter(Hashtag != c('inconclusive'))

epcam_osm_expr <- ggbarplot(data.df.epcam, x = "celltypes", y = "Osm", add = "mean_se",
          color = "experiment", palette = c('red', 'black'), 
          position = position_dodge(-0.8))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  scale_x_discrete(limits = c('Enterocytes_1_0',
                              'Enterocytes_1_1',
                              'Enterocytes_2',
                              'Enterocytes_3',
                              'Enterocytes_4',
                              'Goblet',
                              'TA',
                              'Enteroendocrine',
                              'Proliferating epithelial cells',
                              'Tuft')) +
  stat_compare_means(aes(group = celltypes_exp), label = "p.format", label.y = 0.2)


#FigS1F

osmr_stroma <- FetchData(Stroma, vars=c('Osmr', 'Hashtag', 'experiment', 'celltypes', 'celltypes_exp'))

data.df.stroma <- osmr_stroma %>% group_by(Hashtag, experiment, celltypes, celltypes_exp) %>% 
  summarise(value=mean(Osmr))
colnames(data.df.stroma)[5] <- "Osmr"

data.df.stroma <- data.df.stroma %>% filter(Hashtag != c('inconclusive'))

data.df.stroma[nrow(data.df.stroma) + 1,] = list("Hashtag","Stroma_G3", 'Lymphatic endothelial cells', 'Lymphatic endothelial cells_Stroma_G3', 0)
data.df.stroma[nrow(data.df.stroma) + 1,] = list("Hashtag","Stroma_G3", 'Oligodendrocytes', 'Oligodendrocytes_Stroma_G3', 0)
data.df.stroma[nrow(data.df.stroma) + 1,] = list("Hashtag","Stroma_G3", 'Proliferating cells', 'Proliferating cells_Stroma_G3', 0)
data.df.stroma[nrow(data.df.stroma) + 1,] = list("Hashtag","Stroma_G3", 'Stromal cells 1.2', 'Stromal cells 1.2_Stroma_G3', 0)
data.df.stroma[nrow(data.df.stroma) + 1,] = list("Hashtag","Stroma_G3", 'Stromal cells 2', 'Stromal cells 2_Stroma_G3', 0)


stroma_osmr_expr <- ggbarplot(data.df.stroma, x = "celltypes", y = "Osmr", add = "mean_se",
          color = "experiment", palette = c('red', 'black'), 
          position = position_dodge(-0.8))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  scale_x_discrete(limits = c('Endothelial cells 1',
                              'Stromal cells 1.1',
                              'Stromal cells 1.2',
                              'Lymphatic endothelial cells',
                              'Oligodendrocytes',
                              'Endothelial cells 2',
                              'Smooth muscle cells',
                              'Stromal cells 2',
                              'Proliferating cells',
                              'Mesothelial cells')) +
  stat_compare_means(aes(group = celltypes_exp), label = "p.format", label.y = 2)

Idents(Stroma) <- Stroma@meta.data$celltypes_exp
dotplot <- DotPlot(object = Stroma, features = 'Osmr', cols = c("blue", "red", 'grey', 'yellow', 'green'), split.by = 'Hashtag')

df <- dotplot$data

df <- separate_wider_delim(df, cols = id, delim = "_", names = c('celltype', "group", "experiment", 'Hashtag'))

df$celltypes_exp <- paste(df$celltype, df$group, df$experiment, sep='_')
df <- df %>% filter(Hashtag != c('inconclusive'))
df[nrow(df) + 1,] = list(0.00000000, 0.0000000, 'Osmr', 'Lymphatic endothelial cells', 'Stroma', 'G3', 'Hashtag', 1, '#BEBEBE', 'Lymphatic endothelial cells_Stroma_G3')
df[nrow(df) + 1,] = list(0.00000000, 0.0000000, 'Osmr', 'Oligodendrocytes', 'Stroma', 'G3', 'Hashtag', 1, '#BEBEBE', 'Oligodendrocytes_Stroma_G3')
df[nrow(df) + 1,] = list(0.00000000, 0.0000000, 'Osmr', 'Proliferating cells', 'Stroma', 'G3', 'Hashtag', 1, '#BEBEBE', 'Proliferating cells_Stroma_G3')
df[nrow(df) + 1,] = list(0.00000000, 0.0000000, 'Osmr', 'Stromal cells 1.2', 'Stroma', 'G3', 'Hashtag', 1, '#BEBEBE', 'Stromal cells 1.2_Stroma_G3')
df[nrow(df) + 1,] = list(0.00000000, 0.0000000, 'Osmr', 'Stromal cells 2', 'Stroma', 'G3', 'Hashtag', 1, '#BEBEBE', 'Stromal cells 2_Stroma_G3')

stroma_osmr_cells <- ggbarplot(df, x = "celltype", y = "pct.exp", add = "mean_se",
          color = "experiment", palette = c('red', 'black'), 
          position = position_dodge(-0.8))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  scale_x_discrete(limits = c('Endothelial cells 1',
                              'Stromal cells 1.1',
                              'Stromal cells 1.2',
                              'Lymphatic endothelial cells',
                              'Oligodendrocytes',
                              'Endothelial cells 2',
                              'Smooth muscle cells',
                              'Stromal cells 2',
                              'Proliferating cells',
                              'Mesothelial cells')) +
  stat_compare_means(aes(group = celltypes_exp), label = "p.format", label.y = 100)

#FigS1G

osmr_LPMC <- FetchData(LPMC, vars=c('Osmr', 'Hashtag', 'experiment', 'celltypes', 'celltypes_exp'))

data.df.2 <- osmr_LPMC %>% group_by(Hashtag, experiment, celltypes, celltypes_exp) %>% 
  summarise(value=mean(Osmr))

data.df.2_f <- data.df.2 %>% filter(Hashtag != c('inconclusive'))

data.df_groupes <- osmr_LPMC %>% group_by(Hashtag, experiment, celltypes, celltypes_exp)

colnames(data.df.2_f)[5] <- "Osmr"

data.df_groupes <- osmr_LPMC %>% group_by(Hashtag, experiment, celltypes, celltypes_exp)

data.df.2_f[nrow(data.df.2_f) + 1,] = list("Hashtag","LPMCs_G3", 'Neutrophils', 'Neutrophils_LPMCs_G3', 0)
data.df.2_f[nrow(data.df.2_f) + 1,] = list("Hashtag","LPMCs_G2", 'Plasma cells', 'Plasma cells_LPMCs_G2', 0)
data.df.2_f[nrow(data.df.2_f) + 1,] = list("Hashtag","LPMCs_G2", 'Ilc3', 'Ilc3_LPMCs_G2', 0)


lpmc_osmr_expr <- ggbarplot(data.df.2_f, x = "celltypes", y = "Osmr", add = "mean_se",
          color = "experiment",
          palette = c('red', 'black'),
          position = position_dodge(-0.8))+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))+
  scale_x_discrete(limits = c("Inflammatory Monocytes", "Neutrophils", "Bcells", "DCs", "Cd4",
                              "NK cells", "Ilc1", "Ilc2", "Ilc3", "Tregs", "Cd8", "Cd8 effector",
                              "Plasma cells", "Proliferating", "gd-tcells")) +
  stat_compare_means(aes(group = celltypes_exp), label = "p.format", label.y = 2)
