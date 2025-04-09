###load packages 
library(Seurat)
library(biomaRt)
library(ktplots)
library(SingleCellExperiment)
library(AnnotationDbi)
library(Orthology.eg.db)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(tidyverse)

##downsampling

all_data <- readRDS('data/sc_final/all_data_final.rds')

Idents(all_data) <- all_data$celltypes

all_data_downsamples <- all_data[, sample(colnames(all_data), size =10000, replace=F)]

all_data_norm_counts <- all_data_downsamples@assays[['RNA']]@data

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

mapfun <- function(mousegenes){
  gns <- mapIds(org.Mm.eg.db, mousegenes, "ENTREZID", "SYMBOL")
  mapped <- AnnotationDbi::select(Orthology.eg.db, keys = gns, columns = c("Homo_sapiens", "Mus_musculus"), keytype = "Mus_musculus")
  naind <- is.na(mapped$Homo_sapiens)
  hsymb <- mapIds(org.Hs.eg.db, as.character(mapped$Homo_sapiens[!naind]), "SYMBOL", "ENTREZID")
  out <- data.frame(Mouse_symbol = mousegenes, mapped, Human_symbol = NA)
  out$Human_symbol[!naind] <- hsymb
  out
}

genesV2 = mapfun(rownames(all_data_norm_counts))
all_data_norm_counts_ed <- merge(all_data_norm_counts, genesV2, by.x = 0, by.y = "Mouse_symbol")
all_data_norm_counts_ed <- na.omit(all_data_norm_counts_ed)
rownames(all_data_norm_counts_ed) <- all_data_norm_counts_ed$Human_symbol
all_data_norm_counts_ed[c(1, 10002, 10003, 10004)] <- NULL
write.table(all_data_norm_counts_ed, 'cellphonedb_count_ed_lpmc.txt', sep='\t', quote=F, row.names = T)

meta_data <- all_data_downsamples@meta.data

meta_data <- cbind(rownames(all_data_downsamples@meta.data), all_data_downsamples@meta.data[,'celltypes', drop=F])

names(meta_data)[names(meta_data) == 'celltypes'] <- 'cluster'

##Meta.txt file containing information of cell annotation  
write.table(meta_data, 'cellphonedb_meta_wtDSS_ed_lpmc.txt', sep='\t', quote=F, row.names=F)



###Run CellPhoneDB in terminal 
#activate conda environment 
#source /.../cell_phone_db/bin/activate
#generate new directory for CellPhoneDB output 
#mkdir /.../ouputs_wtDSS 
#run CellPhoneDB 
#cellphonedb method statistical_analysis /.../cellphonedb_meta_wtDSS_ed_lpmc.txt /.../cellphonedb_count_ed_lpmc.txt --output-path=/.../ouputs_wtDSS --project-name="Mouse_OSMR" --counts-data=symbol 

#Ext_4F

means <- read.delim(file = "data/cellphonedb_output/means.txt", check.names = FALSE)
pvals <- read.delim(file = "data/cellphonedb_output/pvalues.txt", check.names = FALSE)
decon <- read.delim(file = 'data/cellphonedb_output/deconvoluted.txt', check.names = FALSE)

plot_cpdb_heatmap(pvals=pvals, cell_types=c("Inflammatory Monocytes", "Neutrophils", "Bcells", "DCs",
                                            "Cd4", "NK cells", "Ilc1", "Tregs",
                                            "Cd8", "Cd8 effector", "Ilc3", "Plasma cells",
                                            "Proliferating", "gd-tcells", "Ilc2", 'Enterocytes_4',
                                            'Enterocytes_3', 'Enterocytes_2', 'Enterocytes_1_1', 'Proliferating epithelial cells',
                                            'Enterocytes_1_0', 'TA', 'Goblet', 'Enteroendocrine', 'Tuft'), 
                  cellheight = 10, cellwidth = 10)

plot_cpdb3(
  scdata = SCE,
  cell_type1 = "Enterocytes_1_0|Enterocytes_1_1|Enterocytes_2|Enterocytes_3|Enterocytes_4",
  cell_type2 = 'Inflammatory Monocytes|Neutrophils|DCs',
  celltype_key = "celltypes", 
  means = means,
  pvals = pvals,
  deconvoluted = decon 
)