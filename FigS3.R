#For the access to the preliminary code as well as the dataset please check the CD (Kong et al.) and UC (Smilie et al. datasets)

#FigS5A

DotPlot(epi.seur, features = 'OSMR', cols=c('red', 'blue', 'white'), split.by = 'Health')
DotPlot(imm.seur, features = c('OSM', 'OSMR'), cols=c('red', 'blue', 'white'), split.by = 'Health')
#FigS5B

col_epi_dotplot <- DotPlot(object = colon_epi, features = c('OSMR'), cols="RdBu", split.by = 'Type', scale = FALSE)

col_imm_dotplot <- DotPlot(object = colon_imm, features = c('OSM', 'OSMR'), cols="RdBu", split.by = 'Type', scale = FALSE)

col_epi_dotplot + col_imm_dotplot

ileum_epi_dotplot <- DotPlot(object = ileum_epi, features = c('OSMR'), cols="RdBu", split.by = 'Type', scale = FALSE)

ileum_imm_dotplot <- DotPlot(object = ileum_imm, features = c('OSM', 'OSMR'), cols="RdBu", split.by = 'Type', scale = FALSE)



