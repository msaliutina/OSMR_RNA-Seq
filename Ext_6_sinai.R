library(GEOquery)
library(limma)
library(tidyverse)
library(GSVA)
library(GSEABase)
library(DESeq2)


sinai_uc <- load('data/corr_analysis/sinai/rawreads_colon_uc.Rdata')

sinai_uc <- get(sinai_uc)

rownames_uc <- rownames(sinai_uc)

sinai_uc <- as.data.frame(lapply(sinai_uc, function(x) 
  if(is.numeric(x)) round(x, 0) else x))

rownames(sinai_uc) <- rownames_uc

sinai_cd <- load('data/corr_analysis/sinai/rawreads_ileum_cd.Rdata')

sinai_cd <- get(sinai_cd)

rownames_cd <- rownames(sinai_cd)
sinai_cd <- as.data.frame(lapply(sinai_cd, function(x) 
  if(is.numeric(x)) round(x, 0) else x))

rownames(sinai_cd) <- rownames_cd

sinai_uc_cd <- merge(sinai_uc, sinai_cd, by="row.names")

diagnosis <- c(rep('UC', 293), rep('CD', 162))

rownames(sinai_uc_cd) <- sinai_uc_cd$Row.names
sinai_uc_cd <- sinai_uc_cd[, -1]

design <- data.frame(colnames(sinai_uc_cd), diagnosis)

dds <- DESeqDataSetFromMatrix(countData = sinai_uc_cd,
                              colData = design,
                              design = ~ diagnosis)

dds <- vst(dds)

vst_df_sinai <- assay(dds) 

#split the cohort to UC and CD subsets
vst_sinai_uc <- vst_df_sinai[, c(1:293)]
vst_sinai_cd <- vst_df_sinai[, c(294:455)]


geneSetList <- list(Top50 = c(
  "DMBT1", "SERPINA3", "REG1A", "REG1B", "HTRA1", "NOS2", "SOCS3", "CFI",
  "DSG3", "PLA2G2A", "TNIP3", "TIFA", "HTR3A", "LYPD1", "IFITM3", "TGM2",
  "SLC9B2", "IFITM1", "PDZK1IP1", "GRIN3A", "SERPINA1", "ALPPL2", "STMN3",
  "WNT5A", "SLC5A8", "OLFM4", "HEG1", "FSTL1", "FLNC",
  "ENKUR", "C2CD4A", "CARD14", "GPR37", "MMP10", "IL1R1", "CITED4", "PRUNE2",
  "ECEL1P2", "GBP4", "KCNN2", "CASP5", "TRPV6", "KRT6A", "NAV3", "TMEM173",
  "STEAP4",
  "IL13RA2", "PLEKHF1", "LYPD5"
))


enrichmentScores <- gsva(vst_sinai_uc, geneSetList, method='ssgsea')

enrichment_vector <- enrichmentScores["Top50", ]

gene_expression_vector <- vst_sinai_uc["OSMR", ]

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
  geom_point(color = "darkgrey", size = 2) +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  labs(
    title = "Sinai cohort, UC",
    subtitle = paste("Pearson r =", r_value, ", p =", p_value),
    x = "Enrichment Score (Top50 gene set)",
    y = "GeneX Expression"
  ) +
  theme_minimal()

#same for CD samples

enrichmentScores <- gsva(vst_sinai_cd, geneSetList, method='ssgsea')

enrichment_vector <- enrichmentScores["Top50", ]

gene_expression_vector <- vst_sinai_cd["OSMR", ]

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
  geom_point(color = "darkgrey", size = 2) +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  labs(
    title = "Sinai cohort, CD",
    subtitle = paste("Pearson r =", r_value, ", p =", p_value),
    x = "Enrichment Score (Top50 gene set)",
    y = "GeneX Expression"
  ) +
  theme_minimal()




