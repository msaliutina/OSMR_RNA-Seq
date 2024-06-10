library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(ggplot2)
library(patchwork)
library(progeny)
library(pheatmap)
library(ComplexHeatmap)
library(colorRamp2)
library(dendextend)
library(EnhancedVolcano)
library(clusterProfiler)
library(enrichplot)
require(DOSE)
library(ggVennDiagram)
organism = 'org.Mm.eg.db'
library(organism, character.only = TRUE)


EpCam <- readRDS('~/Desktop/Hegazy lab/single_cell_2023/Seurat objects/EpCam_harmony_f.rds')

#Fig6A

Idents(EpCam) <- EpCam@meta.data$celltypes

enterocytes <- subset(EpCam, idents = c('Enterocytes_1_0', 'Enterocytes_1_1',
                                        'Enterocytes_2', 'Enterocytes_3', 'Enterocytes_4'))

Idents(object = enterocytes) <- enterocytes@meta.data$experiment

degs_extracted_all <- FindMarkers(enterocytes, ident.1 = "EpCam_G2", 
                                  ident.2= "EpCam_G3", min.cells.group = 1, min.cells.feature = 1,
                                  min.pct = 0,logfc.threshold = 0, only.pos = FALSE)


genes <- c("Cyp2c55", "Ces2c", "Ahnak2", "Ggh", "Irf8",
           "Scd2", "Mfsd2a", "Tppp3", "Mcl1", "Ms4a12",
           "Ripk3", "Cyp2d12", "Rpl35", "Ifit2", "Sgk1",
           "Pck1",
           "Tnni1",
           "Manf", "Myl7",
           "Slc30a10", "Vegfa",
           "Riiad1", "Lhfpl2", "Sidt1",
           "Ppp1r14a", "Filip1l",
           "Ctss", "Parp3",
           "Cd74",
           "Ccl8",
           "Ero1l",
           "Apol6", "Steap4", "Gbp4",
           "Cd47", "Tbc1d9",
           "C3",
           "Fads2", "Tac1", "Nlrc5", "Ndufc2",
           "Cbr2", "Pdzd3",
           "Adcy2",
           "Timp3", "Efna1", "Duoxa1", "Ifi44", "Tnfaip8l3",
           "Cyp2c69Iyd", "Noct", "Spag17", "Sap30l",
           "Cyp3a25",
           "Cyp3a44", "Angptl4", "Pdk2",
           "Fads1", "Slc41a3",
           "Tnfrsf8", "Gm45145",
           "Wtip", "Cap2", "Taf4b", "Olfr56",
           "Sun3", "Zbtb16", "Dnajc5", "Pkib", "Adgrf3", "S100a8", "Osmr")

volcano_sc <- EnhancedVolcano(degs_extracted_all,
                lab = rownames(degs_extracted_all),
                selectLab = genes,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                legendPosition = 'top')

#Fig6B

ranks <- degs_extracted_all$avg_log2FC
names(ranks) <- rownames(degs_extracted_all)
gene_list<-na.omit(ranks)
gene_list = sort(gene_list, decreasing = TRUE)
gse_sc <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             minGSSize = 3, 
             maxGSSize = 300, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "BH")

pathways_to_show <- c('response to type II interferon',
                      'tumor necrosis factor production',
                      'cellular response to lipopolysaccharide',
                      'neutrophil chemotaxis',
                      'granulocyte chemotaxis',
                      'leukocyte chemotaxis',
                      'cell chemotaxis',
                      'regulation of cytokine-mediated signaling pathway',
                      'regulation of leukocyte mediated immunity',
                      'response to lipopolysaccharide',
                      'cytokine production involved in immune response',
                      'antigen processing and presentation',
                      'type II interferon production',
                      'defense response to bacterium')

gsea_dotplot <- dotplot(gse_sc, showCategory=pathways_to_show, split=".sign") + facet_grid(.~.sign)

#Fig6C

entero_merge_df <- data.frame(Cell = names(Idents(enterocytes)), 
                              CellType = as.character(Idents(enterocytes)),
                              stringsAsFactors = FALSE)

entero_merge_progeny <- progeny(enterocytes, scale=FALSE, organism="Mouse", top=500, perm=1, 
                                return_assay = TRUE)

entero_merge_progeny <- Seurat::ScaleData(entero_merge_progeny, assay = "progeny") 

progeny_scores_epith <- 
  as.data.frame(t(GetAssayData(entero_merge_progeny, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 


progeny_scores_epith <- inner_join(progeny_scores_epith, entero_merge_df)

summarized_progeny_scores_epith <- progeny_scores_epith %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

summarized_progeny_scores_epith <- summarized_progeny_scores_epith %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)
progenyBreaks = c(seq(min(summarized_progeny_scores_epith), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_epith)/paletteLength, 
                      max(summarized_progeny_scores_epith), 
                      length.out=floor(paletteLength/2)))

dend15 <- t(summarized_progeny_scores_epith) %>% dist %>% hclust(method = "average") %>% as.dendrogram

dend15 <- dend15 %>%
  rotate(c('JAK-STAT', 'PI3K', 'MAPK', 'WNT', 'TNFa', 'VEGF', 'EGFR', 'NFkB', 'TGFb', 'Hypoxia',
           'p53', 'Estrogen', 'Trail', 'Androgen'))

pheat_metadata <- data.frame(c('Infamed mice', 'Steady state mice'))
row.names(pheat_metadata) <- colnames(t(summarized_progeny_scores_epith))
colnames(pheat_metadata) <- 'group'
f1 = colorRamp2((progenyBreaks), colorRampPalette(c("Darkblue", "white","red"))(101))

ha = HeatmapAnnotation(df = pheat_metadata)

progeny_sc <- Heatmap(t(summarized_progeny_scores_epith), name = "mat", cluster_rows = dend15,
        col = f1, rect_gp = gpar(col = "black", lwd = 2))

#Fig6D

epith_counts <- read.table("~/Desktop/Hegazy lab/bulk_rnaseq_2023/gene_counts_osmr_epith.txt", sep = ',', header = T,
                           row.names = 'X') 
epith_counts <- epith_counts %>% filter(gene_biotype == 'protein_coding')
epith_counts[c(1, 2)] <- NULL

group <- c(rep('OSMR', 4), rep('WT', 5))
sample <- colnames(epith_counts)
metadata <- data.frame(sample, group)

ddsMat <- DESeqDataSetFromMatrix(countData = epith_counts,
                                 colData = metadata,
                                 design = ~ group)
dds <- estimateSizeFactors(ddsMat)

dds <- DESeq(dds)
rld <- rlog(dds)

res_unshrunken <- results(dds, contrast= c('group', 'OSMR', 'WT'), alpha = 0.05)
res <- lfcShrink(dds, contrast=c('group', 'OSMR', 'WT'), res=res_unshrunken, type='ashr')
write.csv(res, '~/Desktop/Hegazy lab/bulk_rnaseq_2023/epith_degs.csv', row.names = T)

volcano_bulk_epith <- EnhancedVolcano(res,
                                      lab = rownames(res),
                                      x = 'log2FoldChange',
                                      y = 'pvalue',
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE,
                                      legendPosition = 'top')

#Fig6E

original_gene_list <- res$log2FoldChange

names(original_gene_list) <- rownames(res)

gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
gse_epith <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             minGSSize = 3, 
             maxGSSize = 300, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "BH")

pathways_to_show_bulk <- c('response to type II interferon',
                      'neutrophil chemotaxis',
                      'granulocyte chemotaxis',
                      'cell chemotaxis',
                      'type II interferon production',
                      'response to tumor necrosis factor',
                      'defense response to bacterium',
                      'cellular response to lipopolysaccharide',
                      'tumor necrosis factor superfamily cytokine production',
                      'regulation of leukocyte mediated immunity',
                      'acute inflammatory response',
                      'inflammatory response to wounding',
                      'acute-phase response',
                      'myeloid leukocyte migration',
                      'tumor necrosis factor superfamily cytokine production')

gse2 <- pairwise_termsim(gse_epith)
treeplot_epith <- treeplot(gse2, showCategory = pathways_to_show_bulk)

#Fig6F
normalized_epith_counts <- counts(dds, normalized=TRUE)
epith_matrix <- as.matrix(normalized_epith_counts)
calc_epith <- progeny(
  epith_matrix,
  scale = TRUE,
  organism = "Mouse",
  top = 100,
  perm = 1,
  verbose = FALSE,
  z_scores = FALSE,
  get_nulldist = FALSE,
  assay_name = "RNA",
  return_assay = FALSE)

calc_df_epith <- data.frame(calc_epith)
calc_df_epith$sum <- c(rep('VilCre', 4), rep('WT-inflammed', 5))
VilCre <- calc_epith[c(1:4), ]
VilCre <- colMeans(VilCre)
WT <- calc_epith[c(5:9), ]
WT <- colMeans(WT)

calc_means_epith <- data.frame(VilCre, WT)

row_dend = as.dendrogram(hclust(dist(calc_means_epith)))

dend15 <- calc_means_epith %>% dist %>% hclust(method = "average") %>% as.dendrogram

dend15 <- dend15 %>%
  rotate(c('JAK-STAT', 'WNT', 'PI3K', 'NFkB', 'TNFa', 'Trail', 'VEGF', 'Hypoxia', 'EGFR', 'Estrogen',
           'p53', 'TGFb', 'Androgen', 'MAPK'))


f1 = colorRamp2((progenyBreaks), colorRampPalette(c("Darkblue", "white","red"))(101))

ha = HeatmapAnnotation(df = pheat_metadata_vilcre,
                       col = list('group' = c('VillinCre/ERT2 x Osmrfl/fl' = 'black', 'VillinCre/ERT2 x Osmrfl/WT' = 'gold')))

progeny_epith <- Heatmap(calc_means_epith, name = "mat", cluster_rows = dend15,
        top_annotation = ha,
        col = f1,
        rect_gp = gpar(col = "black", lwd = 2))

#Fig6G

res_bulk <- data.frame(res) %>% dplyr::filter(padj < 0.05)
res_single_cell <- degs_extracted_all %>% dplyr::filter(p_val_adj < 0.05)


genes <- list('bulk data' = rownames(res_bulk),
              'single cell data' = rownames(res_single_cell))

venn <- ggVennDiagram(genes, category.names = names(genes))
venn + scale_x_continuous(expand = expansion(mult = .2))

intersection <- Reduce(intersect, genes)

#Fig6H

ego <- enrichGO(gene          = intersection,
                OrgDb         = organism,
                keyType = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

pathways_intersection <-c('regulation of innate immune response',
                          'response to type II interferon',
                          'defense response to virus',
                          'canonical NF-kappaB signal transduction',
                          'type II interferon production',
                          'leukocyte migration',
                          'innate immune response activating cell surface receptor signaling pathway',
                          'antigen processing and presentation of peptide antigen',
                          'regulation of interleukin-1 production',
                          'myeloid leukocyte activation',
                          'regulation of cytokine production involved in inflammatory response',
                          'non-canonical NF-kappaB signal transduction',
                          'leukocyte chemotaxis',
                          'negative regulation of CD4-positive, alpha-beta T cell activation',
                          'antigen processing and presentation')

barplot_GO <- barplot(ego, 
        drop = TRUE, 
        showCategory = pathways_intersection, 
        title = "GO Biological Pathways",
        font.size = 8)

#Fig6I

pathways <- c("granulocyte chemotaxis", "leukocyte chemotaxis", "myeloid leukocyte migration", "neutrophil chemotaxis")

gseaplot2(gse_epith, geneSetID = c(52, 62, 69, 100), pvalue_table = TRUE)

