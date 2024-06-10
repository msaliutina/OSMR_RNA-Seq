library(Seurat)

LPMC <- readRDS('/Users/mariiasaliutina/Desktop/Hegazy lab/single_cell_2023/Seurat objects/LPMC.rds')

LPMC[["percent.mt"]] <- PercentageFeatureSet(LPMC, pattern = "^mt-")
LPMC <- subset(LPMC, subset = percent.mt < 25)
LPMC <- NormalizeData(LPMC)

#The cluster with CD45+ cells will be excluded due to unclear cell markers of this cluster
LPMC_ed <- subset(LPMC, idents = 'CD45+ cells', invert = T)

NK_cells <- FindSubCluster(LPMC_ed, cluster = 'NK cells', resolution = 0.1, graph.name = 'RNA_snn')

NK_cells <- SetIdent(NK_cells, value = 'sub.cluster')

NK_0_sub <- FindSubCluster(NK_cells, cluster = 'NK cells_0', resolution = 0.3, graph.name = 'RNA_snn')

NK_cells_0 <- WhichCells(NK_0_sub, idents = 'NK cells_0')
LPMC_ed <- SetIdent(LPMC_ed, cells = NK_cells_0, value = 'Ilc1')
NK_cells_1 <- WhichCells(NK_0_sub, idents = 'NK cells_1')
LPMC_ed <- SetIdent(LPMC_ed, cells = NK_cells_1, value = 'NK cells')
NK_cells_2 <- WhichCells(NK_0_sub, idents = 'NK cells_2')
LPMC_ed <- SetIdent(LPMC_ed, cells = NK_cells_2, value = 'Ilc3')
NK_cells_3 <- WhichCells(NK_0_sub, idents = 'NK cells_3')
LPMC_ed <- SetIdent(LPMC_ed, cells = NK_cells_3, value = 'Ilc2')
NK_cells_4 <- WhichCells(NK_0_sub, idents = 'NK cells_4')
LPMC_ed <- SetIdent(LPMC_ed, cells = NK_cells_4, value = 'Cd8 effector')


LPMC_ed$celltypes<- Idents(LPMC_ed)

saveRDS(LPMC_ed, 'Seurat objects/LPMC_final.Rds')

