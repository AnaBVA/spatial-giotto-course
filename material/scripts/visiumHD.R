##%######################################################%##
#                                                          #
####              Descarga el set de datos              ####
#                                                          #
##%######################################################%##

# Usaremos el set de datos de cerebro de ratón de la tecnología VisiumHD,
# disponible en: https://www.10xgenomics.com/datasets/visium-hd-cytassist-6p5mm-mouse-brain

# Para descargar, ejecuta estos comandos en la terminal de R:
dir.create("visiumHD_data/")

download.file(url = "https://cf.10xgenomics.com/samples/spatial-exp/4.1.0/Visium_HD_6p5mm_Mouse_Brain/Visium_HD_6p5mm_Mouse_Brain_segmented_outputs.tar.gz",
              destfile = "visiumHD_data/Visium_HD_6p5mm_Mouse_Brain_segmented_outputs.tar.gz")

download.file(url = "https://cf.10xgenomics.com/samples/spatial-exp/4.1.0/Visium_HD_6p5mm_Mouse_Brain/Visium_HD_6p5mm_Mouse_Brain_binned_outputs.tar.gz",
              destfile = "visiumHD_data/Visium_HD_6p5mm_Mouse_Brain_binned_outputs.tar.gz")

# Descomprime los archivos
untar(tarfile = "visiumHD_data/Visium_HD_6p5mm_Mouse_Brain_segmented_outputs.tar.gz", 
      exdir = "visiumHD_data/")

untar(tarfile = "visiumHD_data/Visium_HD_6p5mm_Mouse_Brain_binned_outputs.tar.gz", 
      exdir = "visiumHD_data/")

# Alternativamente, puedes descargar los archivos con este comando en la terminal de linux
# wget https://cf.10xgenomics.com/samples/spatial-exp/4.1.0/Visium_HD_6p5mm_Mouse_Brain/Visium_HD_6p5mm_Mouse_Brain_segmented_outputs.tar.gz
# wget https://cf.10xgenomics.com/samples/spatial-exp/4.1.0/Visium_HD_6p5mm_Mouse_Brain/Visium_HD_6p5mm_Mouse_Brain_binned_outputs.tar.gz

##%######################################################%##
#                                                          #
####       Configura las instrucciones de Giotto        ####
#                                                          #
##%######################################################%##

library(Giotto)

instructions <- createGiottoInstructions(save_dir = "visiumhd_plots",
                                         save_plot = TRUE,
                                         show_plot = TRUE,
                                         return_plot = FALSE,
                                         python_path = NULL # alternativamente, "/local/python/path/python" 
)

##%######################################################%##
#                                                          #
####                Lee los datos usando                ####
####       la resolución de segmentación celular        ####
#                                                          #
##%######################################################%##

visium_hd <- createGiottoVisiumHDObjectCell(
    segmented_outputs_dir = "visiumHD_data/segmented_outputs/",
    load_polygons = "cell",
    load_expression = TRUE,
    load_transcripts = FALSE,
    instructions = instructions)

##%######################################################%##
#                                                          #
####                 Explora el objeto                  ####
#                                                          #
##%######################################################%##

# Cuántas células y genes hay en el objeto? 
visium_hd # 19278, 17542 (features, cells)

spatPlot2D(visium_hd,
           point_size = 1)

pDataDT(visium_hd)

##%######################################################%##
#                                                          #
####                 Control de calidad                 ####
#                                                          #
##%######################################################%##

visiumhd_statistics <- addStatistics(gobject = visium_hd,
                                     expression_values = "raw")

spatPlot2D(gobject = visiumhd_statistics, 
           cell_color = "nr_feats", 
           point_size = 1,
           color_as_factor = FALSE)

# Distribución de genes por célula
filterDistributions(gobject = visiumhd_statistics, 
                    detection = "cells")

##%######################################################%##
#                                                          #
####                      Filtrado                      ####
#                                                          #
##%######################################################%##

visium_hd <- filterGiotto(gobject = visium_hd,
                          expression_threshold = 1,
                          feat_det_in_min_cells = 100,
                          min_det_feats_per_cell = 500,
                          expression_values = "raw",
                          verbose = TRUE)

##%######################################################%##
#                                                          #
####                   Normalización                    ####
#                                                          #
##%######################################################%##

visium_hd <- normalizeGiotto(gobject = visium_hd)

##%######################################################%##
#                                                          #
####                    Estadísticas                    ####
#                                                          #
##%######################################################%##

visium_hd <- addStatistics(gobject = visium_hd)

spatPlot2D(gobject = visium_hd, 
           show_image = FALSE, 
           point_size = 1,
           point_alpha = 0.7,
           cell_color = "nr_feats", 
           color_as_factor = FALSE)

##%######################################################%##
#                                                          #
####                        HVF                         ####
#                                                          #
##%######################################################%##

visium_hd <- calculateHVF(visium_hd,
                          save_plot = TRUE)

##%######################################################%##
#                                                          #
####              Reducción de dimensiones              ####
#                                                          #
##%######################################################%##

visium_hd <- runPCA(visium_hd)

screePlot(gobject = visium_hd, 
          ncp = 30)

plotPCA(gobject = visium_hd)

##%######################################################%##
#                                                          #
####                        UMAP                        ####
#                                                          #
##%######################################################%##

visium_hd <- runUMAP(visium_hd, 
                     dimensions_to_use = 1:10)

plotUMAP(gobject = visium_hd)

##%######################################################%##
#                                                          #
####                        tSNE                        ####
#                                                          #
##%######################################################%##

visium_hd <- runtSNE(visium_hd, 
                     dimensions_to_use = 1:10)

plotTSNE(gobject = visium_hd)

##%######################################################%##
#                                                          #
####                  Vecinos cercanos                  ####
#                                                          #
##%######################################################%##

visium_hd <- createNearestNetwork(gobject = visium_hd, 
                                  dimensions_to_use = 1:10)

##%######################################################%##
#                                                          #
####                 Leiden clustering                  ####
#                                                          #
##%######################################################%##

visium_hd <- doLeidenCluster(gobject = visium_hd, 
                             resolution = 0.4)

spatPlot2D(visium_hd, 
           cell_color = "leiden_clus",
           point_size = 1)

##%######################################################%##
#                                                          #
####         Enriquecimiento de tipos celulares         ####
#                                                          #
##%######################################################%##

gran_markers <- c("Nr3c2", "Gabra5", "Tubgcp2", "Ahcyl2",
                  "Islr2", "Rasl10a", "Tmem114", "Bhlhe22", 
                  "Ntf3", "C1ql2")

oligo_markers <- c("Efhd1", "H2-Ab1", "Enpp6", "Ninj2",
                   "Bmp4", "Tnr", "Hapln2", "Neu4",
                   "Wfdc18", "Ccp110")        

di_mesench_markers <- c("Cartpt", "Scn1a", "Lypd6b",  "Drd5",
                        "Gpr88", "Plcxd2", "Cpne7", "Pou4f1",
                        "Ctxn2", "Wnt4")

sign_matrix <- makeSignMatrixPAGE(sign_names = c("Granule_neurons",
                                                 "Oligo_dendrocytes",
                                                 "di_mesenchephalon"),
                                  sign_list = list(gran_markers,
                                                   oligo_markers,
                                                   di_mesench_markers))

visium_hd <- runPAGEEnrich(gobject = visium_hd, 
                           sign_matrix = sign_matrix)

spatCellPlot2D(gobject = visium_hd,
               spat_enr_names = "PAGE",
               cell_annotation_values = c("Granule_neurons",
                                          "Oligo_dendrocytes",
                                          "di_mesenchephalon"),
               cow_n_col = 2,
               coord_fix_ratio = 1, 
               point_size = 0.5, 
               show_legend = TRUE) 

##%######################################################%##
#                                                          #
####                   Usa MapMyCells                   ####
#                                                          #
##%######################################################%##

# Prepara la matriz de expresión con el formato sugerido en
# https://brain-map.org/bkp/analyze/mapmycells/files

raw_matrix <- getExpression(visium_hd,
                            values = "raw",
                            output = "matrix")

write.csv(t(as.matrix(raw_matrix)), file = "visiumhd_matrix.csv",
          quote = FALSE, row.names = TRUE)

# Usa la herramienta en linea para la anotación celular del cerebro
# https://knowledge.brain-map.org/mapmycells/process/?_gl=1*1yz7j6s*_gcl_au*MTYxNTc1MTE4My4xNzkxNzQ2NzM5Li0uLS4xNzkxNzQ2NzM5LjE4MTAyNDE3NTkuMTc5MTc0NjczOS4xNzkxNzQ3MzU0*_ga*MjgzMjY0NTM1LjE3OTE3NDY3Mzk.*_ga_MFZ8RZNMY4*czE3OTE3NDY3MzgkbzEkZzEkdDE3OTE3NDczODgkajYwJGwwJGgxOTI1MDgwMjU3
# O puedes ejecutar el código de forma local usando este código
# https://github.com/AllenInstitute/cell_type_mapper
# Usa la taxonomía de referencia 10X Whole Mouse Brain y el algoritmo de 
# mapeo jerárquico. Ejecuta la anotación y descarga los resultados.

mapmycells_table <- read.csv("mapmycells/visiumhd_matrixcsv_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1791747544978.csv",
                             skip = 4)
colnames(mapmycells_table)[1] <- "cell_ID"

visium_hd <- addCellMetadata(visium_hd,
                             new_metadata = mapmycells_table[, c("cell_ID", "class_name", "class_bootstrapping_probability")],
                             by_column = TRUE)

spatPlot2D(visium_hd,
           cell_color = "class_name",
           show_image = FALSE,
           point_size = 1,
           cell_color_code = c("01 IT-ET Glut" = "red",
                               "02 NP-CT-L6b Glut" = "red",
                               "03 OB-CR Glut" = "red",
                               "04 DG-IMN Glut" = "red",
                               "13 CNU-HYa Glut" = "red",
                               "14 HY Glut" = "red",
                               "16 HY MM Glut" = "red",
                               "17 MH-LH Glut" = "red",
                               "18 TH Glut" = "red",
                               "19 MB Glut" = "red",
                               "23 P Glut" = "red",
                               "24 MY Glut" = "red",
                               "29 CB Glut" = "red",
                               "05 OB-IMN GABA" = "blue", 
                               "06 CTX-CGE GABA" = "blue",
                               "07 CTX-MGE GABA" = "blue",
                               "08 CNU-MGE GABA" = "blue",
                               "09 CNU-LGE GABA" = "blue",
                               "10 LSX GABA" = "blue",
                               "11 CNU-HYa GABA" = "blue",
                               "12 HY GABA" = "blue",
                               "20 MB GABA" = "blue", 
                               "26 P GABA" = "blue",
                               "27 MY GABA" = "blue",
                               "28 CB GABA" = "blue",
                               "21 MB Dopa" = "pink",
                               "30 Astro-Epen" = "green",
                               "31 OPC-Oligo" = "orange",
                               "32 OEC" = "yellow",
                               "33 Vascular" = "grey",
                               "34 Immune" = "purple"
                               ))

##%######################################################%##
#                                                          #
####    Lectura de datos con diferentes resoluciones    ####
#                                                          #
##%######################################################%##

# NOTA: reinicia tu sesión antes de ejecutar la siguiente función

instructions <- createGiottoInstructions(save_dir = "visiumbins_plots",
                                         save_plot = TRUE,
                                         show_plot = TRUE,
                                         return_plot = FALSE,
                                         python_path = NULL # alternativamente, "/local/python/path/python" 
)

# Resolución bin 16um
visium_hd_16 <- createGiottoVisiumHDObjectBin(
    binned_outputs_dir = "visiumHD_data/binned_outputs/",
    bin = 16,
    load_expression = TRUE,
    load_transcripts = FALSE,
    tissue_only = TRUE,
    instructions = instructions
)

# Explora el objeto, cuántos genes y células fueron detectados en esta
# resolución? 19346, 100979 (features, cells)
visium_hd_16

pDataDT(visium_hd_16)

spatPlot2D(visium_hd_16,
           show_image = TRUE,
           point_size = 0.1)

# Resolución bin 8um
visium_hd_8 <- createGiottoVisiumHDObjectBin(
    binned_outputs_dir = "visiumHD_data/binned_outputs/",
    bin = 8,
    load_expression = TRUE,
    load_transcripts = FALSE,
    tissue_only = TRUE
)

# Explora el objeto, cuántos genes y células fueron detectados en esta
# resolución? 19346, 399612 (features, cells)
visium_hd_8

pDataDT(visium_hd_8)

spatPlot2D(visium_hd_8,
           show_image = TRUE,
           point_size = 0.05)

# Filtrado
visium_hd_16 <- filterGiotto(visium_hd_16,
                             feat_det_in_min_cells = 100,
                             min_det_feats_per_cell = 100)

# Normalizacion
visium_hd_16 <- normalizeGiotto(visium_hd_16,
                                scale_cells = FALSE,
                                scale_feats = FALSE)

# Estadísticas
visium_hd_16 <- addStatistics(visium_hd_16)

spatPlot2D(gobject = visium_hd_16, 
           cell_color = "nr_feats", 
           point_size = 0.5,
           color_as_factor = FALSE)

