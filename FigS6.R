library(Seurat)
library(DESeq2)
library(EnhancedVolcano)
library(ComplexHeatmap)

all_data <- readRDS('/Users/mariiasaliutina/Desktop/Hegazy lab/single_cell_2023/Seurat objects/all_data_final.rds')

#FigS6A

DotPlot(object = all_data, features = c('Vil1', 'Col1a2', 'Cdh5'), cols="RdBu") 

#FigS6D

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

rlog_counts <- rlog(dds)
rlogMatrix_counts <- assay(rlog_counts)

sampleCor <- cor(rlogMatrix_counts)
colnames(sampleCor) <- gsub('_count', '', colnames(sampleCor))
rownames(sampleCor) <- gsub('_count', '', rownames(sampleCor))

dend15 <- sampleCor %>% dist %>% hclust(method = "average") %>% as.dendrogram

dend15 <- dend15 %>%
  rotate(c('P13E1Col_1_6', 'P14E1Vil_2_4', 'P13E1Col_1_4', 'P14E2Vil_2_5', 'P13E1Col_1_5',
           'P13E2Cad_1_1', 'P13E2Cad_1_5', 'P13E2Cad_2_3', 'P13E2Cad_1_3', 'P13E2Cad_1_4',
           'P13E2Cad_2_4', 'P13E1Col_1_7', 'P14E2Vil_1_2', 'P13E1Col_1_3', 'P14E2Vil_1_3',
           'P14E1Vil_1_3', 'P14E2Vil_1_1', 'P14E1Vil_1_1', 'a9', 'a10', 'a20', 'a19'))

pheat_metadata <- data.frame(c(rep('IECΔOsmr', 5), rep('WT inflamed', 4),
                               rep('EndoΔOsmr', 4),
                               rep('StromaΔOsmr', 5),
                               rep('Steady state', 4)))
row.names(pheat_metadata) <- metadata$sample
row.names(pheat_metadata) <- gsub('_count', '', rownames(pheat_metadata))
colnames(pheat_metadata) <- 'group'
paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

ha = HeatmapAnnotation(df = pheat_metadata)

Heatmap(sampleCor, name = "mat",
        #cluster_rows = dend15,
        #cluster_columns = dend15,
        top_annotation = ha,
        col = myColor)


#FigS6E

res_unshrunken_vil <- results(dds, contrast= c('group', 'VilCre', 'WT'), alpha = 0.05)
res_vil <- lfcShrink(dds, contrast=c('group', 'VilCre', 'WT'), res=res_unshrunken_vil, type='ashr')

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

col_volcano <- EnhancedVolcano(res_col,
                               lab = rownames(res_col),
                               x = 'log2FoldChange',
                               y = 'pvalue',
                               title = NULL,
                               subtitle = "StromaΔOsmr vs steady state",
                               subtitleLabSize = 15,
                               legendPosition = 'right')

col_volcano_infl <- EnhancedVolcano(res_col_infl,
                                    lab = rownames(res_col_infl),
                                    x = 'log2FoldChange',
                                    y = 'pvalue',
                                    title = NULL,
                                    subtitle = "StromaΔOsmr vs WT inflamed",
                                    subtitleLabSize = 15,
                                    legendPosition = 'right')

#FigS6F

cad_volcano <- EnhancedVolcano(res_cad,
                               lab = rownames(res_cad),
                               x = 'log2FoldChange',
                               y = 'pvalue',
                               title = NULL,
                               subtitle = "EndoΔOsmr vs steady state",
                               subtitleLabSize = 15,
                               legendPosition = 'right')

cad_volcano_infl <- EnhancedVolcano(res_cad_infl,
                                    lab = rownames(res_cad_infl),
                                    x = 'log2FoldChange',
                                    y = 'pvalue',
                                    title = NULL,
                                    subtitle = "Cad5Cre/ERT2 x Osmrfl/fl vs WT inflamed",
                                    subtitleLabSize = 15,
                                    legendPosition = 'right')

#FigS6G

vil_volcano <- EnhancedVolcano(res_wt,
                           lab = rownames(res_wt),
                           x = 'log2FoldChange',
                           y = 'pvalue',
                           title = NULL,
                           subtitle = 'VIECΔOsmr vs steady state',
                           subtitleLabSize = 15,
                           legendPosition = 'none')

vil_volcano_infl <- EnhancedVolcano(res_vil_infl,
                                 lab = rownames(res_vil_infl),
                                 x = 'log2FoldChange',
                                 y = 'pvalue',
                                 title = NULL,
                                 selectLab = total_genes,
                                 subtitle = 'IECΔOsmr vs WT inflamed',
                                 subtitleLabSize = 15,
                                 max.overlaps = 200,
                                 legendPosition = 'none')
