library(limma)
library(tidyverse)
library(GSVA)
library(GSEABase)
library(readxl)

#CD data
cd_counts <- read.table('data/corr_analysis/IBDome/TB_ileum_cd.txt')

rownames(cd_counts) <- cd_counts$V1

colnames(cd_counts) <- cd_counts[1, ]

cd_counts <- cd_counts[-1, -1]

cd_counts_m <- as.matrix(cd_counts)

cd_counts_m <- apply(cd_counts_m, 2, as.numeric)
rownames(cd_counts_m) <- rownames(cd_counts)

#IL22 gene signature from Pavlidis et al.
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
enrichmentScores <- gsva(cd_counts_m, geneSetList, method='ssgsea')

# Extract the enrichment score vector for your gene set
enrichment_vector <- enrichmentScores["Top50", ]

# Extract the expression vector for your gene of interest, e.g., "GeneX"
gene_expression_vector <- cd_counts_m["OSMR", ]


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
  geom_point(color = "darkblue", size = 2) +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  labs(
    title = "Correlation between IL-22 Enrichment Score and OSMR v2 Supp.Fig.2",
    subtitle = paste("Pearson r =", r_value, ", p =", p_value),
    x = "Normalized Enrichment Score (Top50 gene set)",
    y = "GeneX Expression"
  ) +
  theme_minimal()

#UC data
uc_counts <- read.table('data/corr_analysis/IBDome/TB_colon_uc.txt')

rownames(uc_counts) <- uc_counts$V1

colnames(uc_counts) <- uc_counts[1, ]

uc_counts <- uc_counts[-1, -1]

uc_counts_m <- as.matrix(uc_counts)

uc_counts_m <- apply(uc_counts_m, 2, as.numeric)
rownames(uc_counts_m) <- rownames(uc_counts)

enrichmentScores <- gsva(uc_counts_m, geneSetList, method='ssgsea')

enrichment_vector <- enrichmentScores["Top50", ]

gene_expression_vector <- uc_counts_m["OSMR", ]


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
  geom_point(color = "darkblue", size = 2) +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  labs(
    title = "Correlation between IL-22 Enrichment Score and OSMR v2 Supp.Fig.2",
    subtitle = paste("Pearson r =", r_value, ", p =", p_value),
    x = "Normalized Enrichment Score (Top50 gene set)",
    y = "GeneX Expression"
  ) +
  theme_minimal()

