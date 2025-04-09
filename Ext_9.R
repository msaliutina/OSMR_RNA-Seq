library(readxl)
library(EnhancedVolcano)
library(clusterProfiler)
library(tidyverse)

hca7 <- read_excel('data/cell_lines/results_Voom_DEG_OSM_CRC.xlsx', 2)
hca7 <- hca7 %>% filter(biotype == "protein_coding")
hca7_24h <- hca7 %>% filter(adj.P.Val < 0.05 & logFC > 0.5)

#Ext_9C
EnhancedVolcano(hca7,
                lab = hca7$symbol,
                x = 'logFC',
                y = 'adj.P.Val',
                title = NULL,
                subtitle = "",
                subtitleLabSize = 15,
                legendPosition = 'right')

hca_genes_ed <- c("Abcb1b", "Ablim2", "Acp3", "Actbl2", "Actrt3", "Adgrf1", "Adgrg6", "Adm", "Akap6", "Ampd3", "Antxr2", 
                  "Anxa10", "Apol6", "Apol8", "Arid5b", "Asph", "Atad5", "Atp2b4", "Atp7b", "Batf3", "Bccip", "Bcl2l15", 
                  "Bcl3", "Bcl6", "Bhlhe40", "Bmal2", "Bmpr1b", "Btbd16", "C1qtnf6", "C1s1", "Calcrl", "Caprin2", "Casp4",
                  "Ccdc68", "Ccdc70", "Ccdc71l", "Cd2ap", "Cd44", "Cd47", "Cd59a", "Cdc42ep4", "Cdcp1", "Cdh3", "Cebpd", "Cfb", 
                  "Cfh", "Cfi", "Chst15", "Ciita", "Cldn1", "Cldn16", "Clip2", "Cox6b2", "Cpm", "Cpne8", "Crispld2", "Crybg3", 
                  "Crygc", "Csf1", "Ctr9", "Ctsl", "Cxcl1", "Cxcl2", "Cxcl3", "Cyp24a1", "Deptor", "Dgke", "Dixdc1", "Dmbt1", 
                  "Dnajc12", "Dpyd", "Dst", "Dtx3l", "Duox1", "Duoxa2", "Dusp23", "Dusp4", "Edrf1", "Eif5a2", "Epas1", "Erp27", 
                  "Esf1", "Ets1", "Etv6", "Eva1c", "Fam89a", "Fcmr", "Fgb", "Fgfr1", "Fgg", "Fyco1", "Gab1", "Galc", "Galnt12", 
                  "Galnt2", "Gbp3", "Gclm", "Gcnt3", "Gna15", "Gng4", "Hdac5", "Heg1", "Helb", "Hif1a", "Hipk3", "Hivep2", "Hk2",
                  "Hmbox1", "Hrh1", "Htr3a", "Htra1", "Icam1", "Ier5l", "Iffo2", "Ifi204", "Ifitm3", "Igfbp1", "Igflr1", "Igsf23",
                  "Il13ra1", "Il13ra2", "Il15ra", "Il18bp", "Il1rap", "Il24", "Il7", "Insc", "Irag1", "Irf1", "Irx2", "Itga1", "Itga3",
                  "Jakmip1", "Jcad", "Kantr", "Kcnk15", "Kcnmb4", "Kcnq3", "Kdsr", "Klhl29", "Krt6a", "Ktn1", "Lactb2", "Lancl3", "Lcn2", 
                  "Letm2", "Lipg", "Lpar1", "Lrg1", "Lurap1l", "Maff", "Map3k5", "Map3k6", "Mbd4", "Mc1r", "Mcc", "Mecom", "Mmp13",
                  "Mtrf1l", "Mvp", "Mx2", "Myo1c", "Naaladl2", "Nampt", "Nceh1", "Ncf4", "Ncoa6", "Nectin2", "Nid1", "Nnmt", "Nr1d1", 
                  "Nr5a2", "Nrg1", "Nrp1", "Osmr", "Oxsr1", "P2rx5", "P2ry6", "P3h2", "Paep", "Parp9", "Pcdh7", "Pcolce2", "Pcsk5", 
                  "Pdp1", "Pdzk1ip1", "Peli2", "Pim3", "Pitpnc1", "Piwil4", "Pkia", "Pkib", "Plat", "Plcxd3", "Plekhs1", "Plin3", "Plpp3",
                  "Ppbp", "Prdm8", "Prkch", "Pros1", "Prss2", "Prss23", "Ptger2", "Ptp4a1", "Ptp4a3", "Ptpn22", "Ptprb", "Ptpre", 
                  "Rab27a", "Rab3a", "Rab3b", "Rab7b", "Ralgapa1", "Raph1", "Rasl11a", "Rassf8", "Rbm7", "Reg3a", "Rexo2", "Rgs5", 
                  "Rlbp1", "Rnf125", "Rnf141", "Rnf213", "Rps6kb1", "Rtel1", "S100a3", "S100a9", "Samd11", "Sbno2", "Scg5", "Sectm1b",
                  "Serpina1c", "Serpinb2", "Serpinb3c", "Serpinb7", "Serpinb8", "Serpinb9", "Sh2d1b1", "Sh3glb1", "Sh3tc1", "Shoc1", 
                  "Slc16a14", "Slc16a7", "Slc22a4", "Slc2a13", "Slc9b2", "Slco2a1", "Snapc1", "Socs2", "Socs3", "Sord", "Sp110", 
                  "Spink1", "Spry4", "St3gal1", "Stat3", "Sting1", "Stom", "Strip2", "Tacc1", "Tars3", "Tbc1d2", "Tbx20", "Tcaf2",
                  "Tead4", "Tgfbi", "Tgm2", "Tlk2", "Tll2", "Tmc5", "Tmed6", "Tmem158", "Tmem30b", "Tmem45b", "Tnfaip8l3", "Tnfrsf8", 
                  "Tnik", "Tnnc1", "Tpk1", "Trabd2b", "Trib2", "Trim15", "Trim40", "Trim58", "Trim69", "Tsc22d1", "Tubb3", 
                  "Tubd1", "Ubash3b", "Vwa5a", "Wfdc3", "Wnt9a", "Xkr9", "Ypel2", "Zfp84")


#Ext_9D
ids_hca <- c(rep('Osm stimulation', 297))
gene_set_hca <- data.frame(ids_hca, hca_genes_ed) 

res_epith <- read.csv('data/epith_degs.csv')

original_gene_list <- res_epith$log2FoldChange

names(original_gene_list) <- res_epith$X

gene_list<-na.omit(original_gene_list)

gene_list = sort(gene_list, decreasing = TRUE)

em2 <- GSEA(gene_list, TERM2GENE = gene_set_hca)

gseaplot2(em2, geneSetID = 1, pvalue_table = TRUE)
