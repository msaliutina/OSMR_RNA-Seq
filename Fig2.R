library(DESeq2)
library(enrichplot)
library(ggplot2)
library(progeny)
library(tidyverse)
library(pheatmap)
library(DOSE)
library(EnhancedVolcano)
library(reshape2)
library(ComplexHeatmap)
library(colorRamp2)
library(dendextend)
library(biomaRt)
library(AnnotationDbi)
library(clusterProfiler)
organism = 'org.Mm.eg.db'
library(organism, character.only = TRUE)


tissue_seq <- read.table("~/Desktop/Hegazy lab/bulk_rnaseq_2023/gene_counts_osmr.txt")

tissue_seq <- tissue_seq %>% filter(str_detect(gene_biotype, 'protein_coding'))

tissue_seq[c(1:3)] <- NULL 

group <- c(rep('VilCre', 5), rep('WT-inflammed', 4), rep('Cad5Cre', 4), 
           rep('Col1a2Cre', 5), rep('WT', 4))
sample <- colnames(tissue_seq)
metadata <- data.frame(sample, group)

ddsMat <- DESeqDataSetFromMatrix(countData = tissue_seq,
                                 colData = metadata,
                                 design = ~ group)
dds <- estimateSizeFactors(ddsMat)

dds <- DESeq(dds)

total_genes <- c('Osm', 'Il11', 'Osmr', 'Il33', 'Pdpn', 'Icam1', 'Col1a1', 'Il24', 'Il6',
                 'Ccl4', 'Ccl3', 'Ccl2', 'Csf3', 'Cxcl11', 'Cxcl10', 'Cxcl9',
                 'Ccl8', 'Ifng', 'Il22', 'Cxcl5', 'Cxcl2', 'Cxcl3', 'Cxcl1', 'Il1a', 'Il1b', 'Il17a',
                 'Tlr2', 'Tlr4', 'Tlr5', 'Nfkb1', 'Rela', 'Nod2', 'Ifng',
                 'Ifngr1', 'Ifngr2',
                 'Stat1', 'Tbx21', 'Il4', 'Il12a', 'Il12b', 'Il12rb2', 'Il12rb1', 
                 'Stat4', 'Il18', 'Il18r1', 'Il18rap', 'Jun', 'Tnf', 'Il6',
                 'Il1a', 'Il1b', 'Tgfb1', 'Tgfb2', 'Tgfb3', 
                 'Smad2', 'Smad3', 'Stat3', 'Il21', 'Il21r', 'Il23a',
                 'Foxp3', 'Il17a', 'Il4ra', 'Il2rg', 'Stat6', 'Gata3',
                 'Il10', 'Nfatc1')

res_unshrunken_vil <- results(dds, contrast= c('group', 'VilCre', 'WT'), alpha = 0.05)
res_vil <- lfcShrink(dds, contrast=c('group', 'VilCre', 'WT'), res=res_unshrunken_vil, type='ashr')

res_unshrunken_wt <- results(dds, contrast= c('group', 'WT-inflammed', 'WT'), alpha = 0.05)
res_wt <- lfcShrink(dds, contrast=c('group', 'WT-inflammed', 'WT'), res=res_unshrunken_wt, type='ashr')

res_unshrunken_vil_infl <- results(dds, contrast= c('group', 'VilCre', 'WT-inflammed'), alpha = 0.05)
res_vil_infl <- lfcShrink(dds, contrast=c('group', 'VilCre', 'WT-inflammed'), res=res_unshrunken_vil_infl, type='ashr')

res_unshrunken_cad <- results(dds, contrast= c('group', 'Cad5Cre', 'WT'), alpha = 0.05)
res_cad <- lfcShrink(dds, contrast=c('group', 'Cad5Cre', 'WT'), res=res_unshrunken_cad, type='ashr')

res_unshrunken_cad_infl <- results(dds, contrast= c('group', 'Cad5Cre', 'WT-inflammed'), alpha = 0.05)
res_cad_infl <- lfcShrink(dds, contrast=c('group', 'Cad5Cre', 'WT-inflammed'), res=res_unshrunken_cad_infl, type='ashr')

res_unshrunken_col <- results(dds, contrast= c('group', 'Col1a2Cre', 'WT'), alpha = 0.05)
res_col <- lfcShrink(dds, contrast=c('group', 'Col1a2Cre', 'WT'), res=res_unshrunken_col, type='ashr')

res_unshrunken_col_infl <- results(dds, contrast= c('group', 'Col1a2Cre', 'WT-inflammed'), alpha = 0.05)
res_col_infl <- lfcShrink(dds, contrast=c('group', 'Col1a2Cre', 'WT-inflammed'), res=res_unshrunken_col_infl, type='ashr')

#Fig2G
comparison <- c('IECΔOsmr vs WT inflamed',
                'IECΔOsmr vs steady state',
                'EndoΔOsmr vs WT inflamed',
                'EndoΔOsmr vs steady state',
                'StromaΔOsmr x Osmrfl/fl vs WT inflamed',
                'StromaΔOsmr x Osmrfl/fl vs steady state',
                'WT inflamed vs steady state')
upregulated <- c(3066, 647, 4, 3360, 328, 4366, 4435)
downregulated <- c(3753, 871, 0, 3706, 463, 4187, 3940)
non_significant <- c(10781, 16082, 17596, 10534, 16809, 9047, 9225)

de_stat <- data.frame(comparison, upregulated, downregulated, non_significant)

de_stat <- melt(de_stat[,c('comparison', 'upregulated','downregulated','non_significant')],id.vars = 1)

ggplot(de_stat,aes(x = comparison,y = value, fill=variable)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_text(aes(label = value),
            position = position_dodge(0.9),
            vjust = -1)

#Fig2H

ibd_genes <- c("Tlr5", "Il12rb2", "Foxp3", "Il4", "Il12a", "Stat6", "Jun", "Smad2", 
                        "Nfkb1", "Smad3", "Nfatc1", "Tgfb2", "Tgfb3", "Col1a1", "Cxcl2", "Cxcl3", 
                        "Ifng", "Tnf", "Cxcl1", "Il1b", "Cxcl10", "Cxcl9", "Cxcl5", "Osmr", "Ifngr1", 
                        "Il23a", "Stat3", "Il18", "Ifngr2", "Rela", "Nod2", "Tlr4", "Tlr2", "Il10", 
                        "Icam1", "Il18r1", "Il33", "Tgfb1", "Pdpn", "Stat1", "Gata3", "Ccl4", 
                        "Il18rap", "Il4ra", "Stat4", "Il21", "Tbx21", "Ccl3", "Il2rg", "Il22", 
                        "Osm", "Il17a", "Il12b", "Il12rb1", "Il21r", "Il11", "Il1a",
                        "Csf3", "Il6", "Ccl2", "Ccl8")


rld <- rlog(dds, blind=FALSE)
mat_nat <- assay(rld)[ibd_genes, ]
mat_nat <- mat_nat - rowMeans(mat_nat)

group <- c(rep('VilCre', 5), rep('WT-inflammed', 4), rep('Cad5Cre', 4), 
           rep('Col1a2Cre', 5), rep('WT', 4))
VilCre <- mat_nat[, c(1:5)]
VilCre <- rowMeans(VilCre)
WT_inflammed <- mat_nat[, c(6:9)]
WT_inflammed <- rowMeans(WT_inflammed)
Cad5Cre <- mat_nat[, c(10:13)]
Cad5Cre <- rowMeans(Cad5Cre)
Col1a2Cre <- mat_nat[, c(14:18)]
Col1a2Cre <- rowMeans(Col1a2Cre)
WT <- mat_nat[, c(19:22)]
WT <- rowMeans(WT)

calc_means <- data.frame(VilCre, WT_inflammed, Cad5Cre, Col1a2Cre,  WT)
pheat_metadata <- data.frame(c('VillinCre/ERT2 x Osmrfl/fl', 'WT_inflammed', 'Cdh5Cre/ERT2 x Osmrfl/fl',
                               'Col1a2Cre/ERT2 x Osmrfl/fl', 'WT'))
row.names(pheat_metadata) <- colnames(calc_means)
colnames(pheat_metadata) <- 'group'

paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

ha = HeatmapAnnotation(df = pheat_metadata)

heatmap <- Heatmap(calc_means, name = "mat",
        #cluster_rows = dend15,
        #cluster_columns = dend15,
        top_annotation = ha,
        col = myColor)

#Fig2I

original_gene_list <- res_vil_infl$log2FoldChange

names(original_gene_list) <- rownames(res_vil_infl)

gene_list<-na.omit(original_gene_list)
gene <- names(gene_list)[abs(gene_list) > 1]


gene_entrez <- AnnotationDbi::select(org.Mm.eg.db, keys=gene, columns='ENTREZID', keytype='SYMBOL')
gene_entrez <- na.omit(gene_entrez)

kk <- enrichKEGG(gene         = gene_entrez$ENTREZID,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05)

pathways_kegg <- c('Cytokine-cytokine receptor interaction - Mus musculus (house mouse)',
                   'IL-17 signaling pathway - Mus musculus (house mouse)',
                   'Chemokine signaling pathway - Mus musculus (house mouse)',
                   'TNF signaling pathway - Mus musculus (house mouse)',
                   'Cell adhesion molecules - Mus musculus (house mouse)',
                   'Inflammatory bowel disease - Mus musculus (house mouse)',
                   'NF-kappa B signaling pathway - Mus musculus (house mouse)',
                   'Th17 cell differentiation - Mus musculus (house mouse)',
                   'JAK-STAT signaling pathway - Mus musculus (house mouse)')

dotplot(kk, showCategory=pathways_kegg)
