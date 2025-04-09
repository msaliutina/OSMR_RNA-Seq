library(Seurat)

all_data <- readRDS('data/sc_final/all_data_final.rds')
EpCam <- readRDS('data/sc_final/EpCam_harmony.rds')
LPMC <- readRDS('data/sc_final/LPMC_final.Rds')
Stroma <- readRDS('data/sc_final/Stroma_final.rds')

#FigS2A

genes_to_dotplot <- c('Il6', 'Il1b', 'Il11', 'Cxcl2', 'Csf3',
                      'Il1a', 'Cxcl3', 'Cxcl1', 'Cxcl5', 'Il17a', 'Ifng', 'Cxcl9', 'Tnf',
                      'Ccl2', 'Cxcl10', 'Csf2', 'Ccl4', 'Ccl3',
                      'Il23a', 'Il12b', 'Il12a', 'Il22', 'Reg3b', 'Reg3g', 'S100a9', 'Cd74', 'Socs1', 
                      'Socs3', 'Mki67', 'Lyz2', 'S100a8', 'H2-Aa', 'H2-Ab1', 'Duoxa2', 'Aqp8', 'Stat1', 'Stat3', 'Isg15', 'Epcam',
                      'Ptprc')

DefaultAssay(all_data) <- 'RNA'

mylevels <- c('Ilc2', 'Ilc1', 'Ilc3', 'gd-tcells', 'Proliferating', 'Plasma cells',
              'Cd8 effector', 'Cd8', 'Tregs', 'NK cells', 'Cd4', 'DCs',
              'Bcells', 'Neutrophils', 'Inflammatory Monocytes',
              'Mesothelial cells', 'Stromal cells 1.1', 'Endothelial cells 2',
              'Stromal cells 2', 'Lymphatic endothelial cells', 'Proliferating cells',
              'Smooth muscle cells', 'Stromal cells 1.2', 'Endothelial cells 1',
              'Oligodendrocytes', 'Enterocytes_3', 'Tuft', 'Enteroendocrine',
              'Enterocytes_0', 'TA', 'Paneth', 'Goblet', 'Enterocytes_1',
              'Enterocytes_2')

dotplot <- DotPlot(object = all_data, features = genes_to_dotplot, split.by = 'experiment', cols="RdBu")
dotplot + theme(axis.text.x = element_text(angle = 45, hjust=1))

#FigS2B

cytokine_receptors <- c('Ackr2', 'Ccrl2', 'Cmklr1', 'Cntfr', 'Crlf1', 'Csf2ra', 'Csf2rb', 'Csf3r', 'Epor',
                        'Flt3', 'Ifnar1', 'Ifnar2', 'Ifngr1', 'Ifngr2', 'Ifnlr1', 'Il10ra', 'Il10rb', 'Il11ra1',
                        'Il12rb1', 'Il12rb2', 'Il13ra1', 'Il13ra2', 'Il15ra', 'Il17ra', 'Il17rb', 'Il17rc', 'Il17rd',
                        'Il17re', 'Il18r1', 'Il18rap', 'Il1r1', 'Il1r2', 'Il1rap', 'Il1rapl2', 'Il1rl1', 'Il1rl2', 'Il20ra',
                        'Il20rb', 'Il21r', 'Il22ra1', 'Il22ra2', 'Il23r', 'Il27ra', 'Il2ra', 'Il2rb', 'Il2rg', 'Il31ra',
                        'Il3ra', 'Il4ra', 'Il5ra', 'Il6ra', 'Il6st', 'Il7r', 'Il9r', 'Leprot', 'Lifr', 'Osmr')

dotplot_epcam <- DotPlot(object = EpCam, features = rev(cytokine_receptors), cols="RdBu", split.by = 'experiment') 
dotplot_epcam + theme(axis.text.x = element_text(angle = 45, hjust=1))


#review


genes_to_dotplot <- c('Il22', 'Il6', 'Il1b', 'Il11', 'Cxcl2', 'Csf3',
                      'Il1a', 'Cxcl3', 'Cxcl1', 'Cxcl5', 'Il17a', 'Ifng', 'Cxcl9', 'Tnf',
                      'Ccl2', 'Cxcl10', 'Csf2', 'Ccl4', 'Ccl3',
                      'Il23a', 'Il12b', 'Il12a', 'Il22', 'Reg3b', 'Reg3g', 'S100a9', 'Cd74', 'Socs1', 
                      'Socs3', 'Mki67', 'Lyz2', 'S100a8', 'H2-Aa', 'H2-Ab1', 'Duoxa2', 'Aqp8', 'Stat1', 'Stat3', 'Isg15', 'Epcam',
                      'Ptprc')

DefaultAssay(all_data) <- 'RNA'

mylevels <- c('Ilc2', 'Ilc1', 'Ilc3', 'gd-tcells', 'Proliferating', 'Plasma cells',
              'Cd8 effector', 'Cd8', 'Tregs', 'NK cells', 'Cd4', 'DCs',
              'Bcells', 'Neutrophils', 'Inflammatory Monocytes',
              'Mesothelial cells', 'Stromal cells 1.1', 'Endothelial cells 2',
              'Stromal cells 2', 'Lymphatic endothelial cells', 'Proliferating cells',
              'Smooth muscle cells', 'Stromal cells 1.2', 'Endothelial cells 1',
              'Oligodendrocytes', 'Enterocytes_3', 'Tuft', 'Enteroendocrine',
              'Enterocytes_0', 'TA', 'Paneth', 'Goblet', 'Enterocytes_1',
              'Enterocytes_2')

dotplot <- DotPlot(object = all_data, features = genes_to_dotplot, split.by = 'experiment', cols="RdBu")
dotplot + theme(axis.text.x = element_text(angle = 45, hjust=1))

