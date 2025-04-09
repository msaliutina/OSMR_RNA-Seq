library(GEOquery)
library(limma)
library(tidyverse)
library(GSVA)
library(GSEABase)
library(readxl)


unifi_ust <- getGEO('GSE206285', GSEMatrix = T)

metadata <- pData(phenoData(unifi_ust[[1]]))

unifi_ust_gse <- unifi_ust[[1]]
expr_matrix_unifi <- exprs(unifi_ust_gse)
feature.data <- unifi_ust$GSE206285_series_matrix.txt.gz@featureData@data

expr_matrix_unifi <- data.frame(expr_matrix_unifi)
expr_matrix_unifi <- expr_matrix_unifi %>%
  rownames_to_column(var = "ID")
expr_matrix_unifi <- inner_join(expr_matrix_unifi, feature.data, by = "ID")

expr_matrix_unifi_mean <- expr_matrix_unifi %>%
  group_by(`Gene Symbol`) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

rownames(expr_matrix_unifi_mean) <- expr_matrix_unifi_mean$`Gene Symbol`

expr_matrix_unifi_mean$`Gene Symbol` <- NULL

geneSetList <- list(Top50 = c( "ACOD1", "CCL7", "CLEC4E", "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL9",
                               "CXCR2", "IDO1", "IL36G", "KRT6B", "LRG1", "LTF", "LY6D", "MFSD2A",
                               "MMP7", "MYOT", "NANOS2", "NLRP12", "OLR1", "PLET1", "PPEF1", "PRSS22",
                               "PRSS27", "REG3G", "RETNLB", "S100A8", "S100A9", "SERPINB11", "SLC30A10",
                               "SPINK6", "TACSTD2", "TMEM92", "TNFRSF8", "TRIM10",
                               "CLCA4C-PS", "CYP2C55", "CYP2C69", "GM4841", "GM5483", "KRT90", "LY6I",
                               "PIRA11", "PLET1OS", "REG3B", "SPRR2H", "SPRR2K", "STFA2", "TGTP1"))

#keep only UC samples
expr_matrix_unifi_uc <- expr_matrix_unifi_mean[, c(1:533)]

expr_matrix_unifi_uc_m <- as.matrix(expr_matrix_unifi_mean)

rownames(expr_matrix_unifi_uc_m) <- rownames(expr_matrix_unifi_mean)
expr_matrix_unifi_uc_m <- expr_matrix_unifi_uc_m[, 2:533]
expr_matrix_unifi_uc_n <- matrix(as.numeric(expr_matrix_unifi_uc_m), nrow = nrow(expr_matrix_unifi_uc_m), ncol = ncol(expr_matrix_unifi_uc_m))
rownames(expr_matrix_unifi_uc_n) <- rownames(expr_matrix_unifi_uc_m)
colnames(expr_matrix_unifi_uc_n) <- colnames(expr_matrix_unifi_uc_m)

geneSetList_fig <- list(Top50 = c(
  "DMBT1", "SERPINA3", "REG1A", "REG1B", "HTRA1", "NOS2", "SOCS3", "CFI",
  "DSG3", "PLA2G2A", "TNIP3", "TIFA", "HTR3A", "LYPD1", "IFITM3", "TGM2",
  "SLC9B2", "IFITM1", "PDZK1IP1", "GRIN3A", "SERPINA1", "ALPPL2", "STMN3",
  "WNT5A", "SLC5A8", "OLFM4", "HEG1", "FSTL1", "FLNC",
  "ENKUR", "C2CD4A", "CARD14", "GPR37", "MMP10", "IL1R1", "CITED4", "PRUNE2",
  "ECEL1P2", "GBP4", "KCNN2", "CASP5", "TRPV6", "KRT6A", "NAV3", "TMEM173",
   "STEAP4",
  "IL13RA2", "PLEKHF1", "LYPD5"
))

enrichmentScores <- gsva(expr_matrix_unifi_uc_n, geneSetList, method='ssgsea')

enrichment_vector <- enrichmentScores["Top50", ]

gene_expression_vector <- expr_matrix_unifi_uc_n["OSMR", ]

min_val <- min(enrichment_vector)
max_val <- max(enrichment_vector)
norm_scores <- (enrichment_vector - min_val) / (max_val - min_val)
enrichment_vector_norm <- norm_scores * 2 - 1

cor_result <- cor.test(enrichment_vector_norm, gene_expression_vector, method = "pearson")
print(cor_result)


plot_data <- data.frame(
  Enrichment = enrichment_vector_norm,
  Expression = gene_expression_vector
)

r_value <- round(cor_result$estimate, 2)
p_value <- round(cor_result$p.value, 3)

p <- ggplot(plot_data, aes(x = Enrichment, y = Expression)) +
  geom_point(color = "darkred", size = 2) +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  labs(
    title = "Correlation between IL-22 Enrichment Score and OSMR v2 Supp.Fig.2",
    subtitle = paste("Pearson r =", r_value, ", p =", p_value),
    x = "Enrichment Score (Top50 gene set)",
    y = "GeneX Expression"
  ) +
  theme_minimal()


