library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggalluvial)
library(patchwork)

all_data <- readRDS('~/Desktop/Hegazy lab/single_cell_2023/Seurat objects/all_data_final.rds')
EpCam <- readRDS('~/Desktop/Hegazy lab/single_cell_2023/Seurat objects/EpCam_harmony.rds')
LPMC <- readRDS('~/Desktop/Hegazy lab/single_cell_2023/Seurat objects/LPMC_final.Rds')
Stroma <- readRDS('~/Desktop/Hegazy lab/single_cell_2023/Seurat objects/Stroma_final.rds')

#Fig1E

Idents(all_data) <- all_data$cells

umap_all_data <- DimPlot(all_data,
        reduction = "umap") 

#Fig1F

Idents(Stroma) <- Stroma$celltypes
Idents(LPMC) <- LPMC$celltypes
Idents(EpCam) <- EpCam$celltypes

stromal_umap <- DimPlot(Stroma,
                        reduction = "umap",
                        label = FALSE) 

epcam_umap <- DimPlot(EpCam,
                      reduction = "umap",
                      label = FALSE) 


lpmc_umap <- DimPlot(LPMC,
                     reduction = "umap",
                     label = FALSE) 

lpmc_umap + epcam_umap + stromal_umap 

#Fig1G

#Fig1G

LPMC_df <- LPMC@meta.data %>% group_by(experiment, celltypes) %>% 
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

Stroma_df <- Stroma@meta.data %>% group_by(experiment, celltypes) %>% 
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

EpCam_df <- EpCam@meta.data %>% group_by(experiment, celltypes) %>% 
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)


coul <- brewer.pal(10, "PRGn") 

coul <- colorRampPalette(coul)(14)

lpmc_plot <- ggplot(LPMC_df, aes(x = experiment, y = percent, fill = celltypes))+
  scale_x_discrete(limits = c('LPMCs_G3', 'LPMCs_G2')) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c('green4', 'hotpink3', 'pink', 'cornflowerblue', 'blue', 'lightblue',
                             'turquoise1', 'deepskyblue4', 'slategray3', 'slategray1', 'cornsilk', 
                             'gray', 'coral', 'darkolivegreen4', 'violet', 'goldenrod')) +
  geom_flow(aes(alluvium = celltypes), alpha= .9, 
            lty = 2, fill = "white", color = "black",
            curve_type = "linear", 
            width = .5) +
  geom_text(aes(label=paste0(sprintf("%1.1f", percent),"%")),
            position=position_stack(vjust=0.5)) 

stroma_colors <- colorRamp(brewer.pal(10,'Set3'))(10)
stroma_plot <- ggplot(Stroma_df, aes(x = experiment, y = percent, fill = celltypes))+
  scale_x_discrete(limits = c('Stroma_G3', 'Stroma_G2')) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "RdGy") +
  geom_flow(aes(alluvium = celltypes), alpha= .9, 
            lty = 2, fill = "white", color = "black",
            curve_type = "linear", 
            width = .5) +
  geom_text(aes(label=paste0(sprintf("%1.1f", percent),"%")),
            position=position_stack(vjust=0.5)) 


epcam_plot <- ggplot(EpCam_df, aes(x = experiment, y = percent, fill = celltypes))+
  scale_x_discrete(limits = c('EpCam_G3', 'EpCam_G2')) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "BrBG") +
  geom_flow(aes(alluvium = celltypes), alpha= .9, 
            lty = 2, fill = "white", color = "black",
            curve_type = "linear", 
            width = .5) +
  geom_text(aes(label=paste0(sprintf("%1.1f", percent),"%")),
            position=position_stack(vjust=0.5)) 

lpmc_plot +  epcam_plot + stroma_plot

#Fig1H

LPMC@meta.data$celltypes_exp <- paste(LPMC@meta.data$celltypes, LPMC@meta.data$experiment, sep='_')
osm_LPMC <- FetchData(LPMC, vars=c('Osm', 'Hashtag', 'experiment', 'celltypes', 'celltypes_exp'))

data.df.2 <- osm_LPMC %>% group_by(Hashtag, experiment, celltypes, celltypes_exp) %>% 
  summarise(value=mean(Osm))

data.df.2_f <- data.df.2 %>% filter(Hashtag != c('inconclusive'))

data.df_groupes <- osm_LPMC %>% group_by(Hashtag, experiment, celltypes, celltypes_exp)

colnames(data.df.2_f)[5] <- "Osm"

data.df_groupes <- osm_LPMC %>% group_by(Hashtag, experiment, celltypes, celltypes_exp)


###add lacking rows
data.df.2_f[nrow(data.df.2_f) + 1,] = list("Hashtag","LPMCs_G3", 'Neutrophils', 'Neutrophils_LPMCs_G3', 0)
data.df.2_f[nrow(data.df.2_f) + 1,] = list("Hashtag","LPMCs_G2", 'Plasma cells', 'Plasma cells_LPMCs_G2', 0)
data.df.2_f[nrow(data.df.2_f) + 1,] = list("Hashtag","LPMCs_G2", 'Ilc3', 'Ilc3_LPMCs_G2', 0)


lpmc_osm_expr <- ggbarplot(data.df.2_f, x = "celltypes", y = "Osm", add = "mean_se",
          color = "experiment",
          palette = c('red', 'black'),
          position = position_dodge(-0.8))+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))+
  scale_x_discrete(limits = c("Inflammatory Monocytes", "Neutrophils", "Bcells", "DCs", "Cd4",
                              "NK cells", "Ilc1", "Ilc2", "Ilc3", "Tregs", "Cd8", "Cd8 effector",
                              "Plasma cells", "Proliferating", "gd-tcells")) +
  stat_compare_means(aes(group = celltypes_exp), label = "p.format", label.y = 2)


Idents(LPMC) <- LPMC@meta.data$celltypes_exp
dotplot <- DotPlot(object = LPMC, features = 'Osm', cols = c("blue", "red", 'grey', 'yellow', 'green'), split.by = 'Hashtag')

df <- dotplot$data

df <- df[c(-4), ] 

df <- separate_wider_delim(df, cols = id, delim = "_", names = c('celltype', "group", "experiment", 'Hashtag'))

df$celltypes_exp <- paste(df$celltype, df$group, df$experiment, sep='_')
df <- df %>% filter(Hashtag != c('inconclusive'))
df[nrow(df) + 1,] = list(0.00000000, 0.0000000, 'Osm', 'Neutrophils', 'LPMCs', 'G3', 'Hashtag', 1, '#BEBEBE', 'Neutrophils_LPMCs_G3')
df[nrow(df) + 1,] = list(0.00000000, 0.0000000, 'Osm', 'Plasma cells', 'LPMCs', 'G2', 'Hashtag', 1, '#BEBEBE', 'Plasma cells_LPMCs_G2')
df[nrow(df) + 1,] = list(0.00000000, 0.0000000, 'Osm', 'Ilc3', 'LPMCs', 'G2', 'Hashtag', 1, '#BEBEBE', 'Ilc3_LPMCs_G2')


lpmc_osm_cells <- ggbarplot(df, x = "celltype", y = "pct.exp", add = "mean_se",
          color = "experiment",
          palette = c('red', 'black'),
          position = position_dodge(-0.8))+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))+
  scale_x_discrete(limits = c("Inflammatory Monocytes", "Neutrophils", "Bcells", "DCs", "Cd4",
                              "NK cells", "Ilc1", "Ilc2", "Ilc3", "Tregs", "Cd8", "Cd8 effector",
                              "Plasma cells", "Proliferating", "gd-tcells")) +
  stat_compare_means(aes(group = celltypes_exp), label = "p.format", label.y = 60)


#Fig1I

EpCam@meta.data$celltypes_exp <- paste(EpCam@meta.data$celltypes, EpCam@meta.data$experiment, sep='_')

osmr_epcam <- FetchData(EpCam, vars=c('Osmr', 'Hashtag', 'experiment', 'celltypes', 'celltypes_exp'))

data.df.epcam <- osmr_epcam %>% group_by(Hashtag, experiment, celltypes, celltypes_exp) %>% 
  summarise(value=mean(Osmr))
colnames(data.df.epcam)[5] <- "Osmr"

data.df.epcam <- data.df.epcam %>% filter(Hashtag != c('inconclusive'))

epcam_osmr_expr <- ggbarplot(data.df.epcam, x = "celltypes", y = "Osmr", add = "mean_se",
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
  stat_compare_means(aes(group = celltypes_exp), label = "p.format", label.y = 0.3)

Idents(EpCam) <- EpCam@meta.data$celltypes_exp
dotplot <- DotPlot(object = EpCam, features = 'Osmr', cols = c("blue", "red", 'grey', 'yellow', 'green'), split.by = 'Hashtag')

df <- dotplot$data

df$id <- gsub('Enterocytes_', 'Enterocytes.', df$id)
df$id <- gsub('Enterocytes.1_', 'Enterocytes.1.', df$id)



df <- separate_wider_delim(df, cols = id, delim = "_", names = c('celltype', "group", "experiment", 'Hashtag'))

df$celltypes_exp <- paste(df$celltype, df$group, df$experiment, sep='_')
df <- df %>% filter(Hashtag != c('inconclusive'))
df$celltype <- gsub('Enterocytes.', 'Enterocytes_', df$celltype)
df$celltypes_exp <- gsub('Enterocytes.', 'Enterocytes_', df$celltypes_exp)
df$celltype <- gsub('1.0', '1_0', df$celltype)
df$celltypes_exp <- gsub('1.0', '1_0', df$celltypes_exp)
df$celltype <- gsub('1.1', '1_1', df$celltype)
df$celltypes_exp <- gsub('1.1', '1_1', df$celltypes_exp)


epcam_osmr_cells <- ggbarplot(df, x = "celltype", y = "pct.exp", add = "mean_se",
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
  stat_compare_means(aes(group = celltypes_exp), label = "p.format", label.y = 25)