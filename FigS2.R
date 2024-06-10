library(Seurat)

EpCam <- readRDS('~/Desktop/Hegazy lab/single_cell_2023/Seurat objects/EpCam_harmony.rds')
LPMC <- readRDS('~/Desktop/Hegazy lab/single_cell_2023/Seurat objects/LPMC_final.Rds')
Stroma <- readRDS('~/Desktop/Hegazy lab/single_cell_2023/Seurat objects/Stroma_final.rds')


#FigS2A

Idents(LPMC) <- LPMC$celltypes

LPMC <- ScaleData(object = LPMC_ed, vars.to.regress = 'percent.mt', features = rownames(LPMC_ed), 
                     block.size = 2000)

lpmc.markers <- FindAllMarkers(object = LPMC, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)
top10_LPMC <- lpmc.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)

lpmc_heatmap <- DoHeatmap(object = LPMC, features = top10_LPMC$gene, label = F)

#FigS2B

Idents(EpCam) <- EpCam$celltypes

EpCam <- ScaleData(object = EpCam, vars.to.regress = 'percent.mt', features = rownames(EpCam), 
                   block.size = 2000)
epcam.markers <- FindAllMarkers(object = EpCam, only.pos = TRUE, min.pct = 0.25, 
                                thresh.use = 0.25)
top10_epcam <- epcam.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)

epcam_heatmap <- DoHeatmap(object = EpCam, features = top10_epcam$gene, label = F)


#FigS2C
Idents(Stroma) <- Stroma$celltypes

Stroma <- ScaleData(object = Stroma, vars.to.regress = 'percent.mt', features = rownames(Stroma), 
                    block.size = 2000)

stroma.markers <- FindAllMarkers(object = Stroma, only.pos = TRUE, min.pct = 0.25, 
                                 thresh.use = 0.25)
top10_stroma <- stroma.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)

stroma_heatmap <- DoHeatmap(object = Stroma, features = top10_stroma$gene, label = FALSE)