##%######################################################%##
#                                                          #
####              Descarga el set de datos              ####
#                                                          #
##%######################################################%##

# Usaremos el set de datos de cancer de mama de la tecnología Xenium,
# disponible en: https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast

# Para descargar, ejecuta estos comandos en la terminal de R:
dir.create("xenium_data/")

download.file(url = "https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip",
              destfile = "xenium_data/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip")

# Descomprime los archivos
untar(tarfile = "xenium_data/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip", 
      exdir = "xenium_data/")

##%######################################################%##
#                                                          #
####       Configura las instrucciones de Giotto        ####
#                                                          #
##%######################################################%##

library(Giotto)

instructions <- createGiottoInstructions(
    save_dir = "xenium_plots",
    save_plot = TRUE,
    show_plot = TRUE,
    return_plot = FALSE,
    python_path = NULL
)

##%######################################################%##
#                                                          #
####               Crea el objeto Giotto                ####
#                                                          #
##%######################################################%##

data_path <- "xenium_data/outs/"

# Si tienes limitaciones de memoria RAM, usa esta versión sin transcritos
xenium_gobj <- createGiottoXeniumObject(
    xenium_dir = data_path,
    qv_threshold = 20, # qv of 20 is the default and also what 10x uses
    feat_type = "rna",
    load_expression = TRUE,
    load_transcripts = FALSE,
    instructions = instructions
)

# Si tienes mayor memoria RAM, usa esta versión para leer los transcritos
xenium_gobj <- createGiottoXeniumObject(
    xenium_dir = data_path,
    qv_threshold = 20, # qv of 20 is the default and also what 10x uses
    feat_type = "rna",
    load_expression = FALSE,
    load_transcripts = TRUE,
    instructions = instructions
)

##%######################################################%##
#                                                          #
####                 Explora el objeto                  ####
#                                                          #
##%######################################################%##

xenium_gobj

spatPlot2D(xenium_gobj,
           plot_method = "scattermore", # scattermore acelera la creación de la gráfica
           point_shape = "no_border",
           point_size = 1
)

##%######################################################%##
#                                                          #
####             Agrega los transcritos que             ####
####             corresponden a cada célula             ####
#                                                          #
##%######################################################%##

# Calcula el sobrelape de los transcritos (rna) con el área de las células
xenium_gobj <- calculateOverlap(xenium_gobj,
                                spat_info = "cell",
                                feat_info = "rna")

# Crea la matriz de expresión usando los transcritos que sobrelapan con
# las células
xenium_gobj <- overlapToMatrix(xenium_gobj,
                               poly_info = "cell",
                               feat_info = "rna",
                               name = "raw")

# Explora de nuevo el objeto, notas alguna diferencia?
xenium_gobj

##%######################################################%##
#                                                          #
####                      Filtrado                      ####
#                                                          #
##%######################################################%##

xenium_gobj <- filterGiotto(xenium_gobj,
                            spat_unit = "cell",
                            expression_threshold = 1,
                            feat_det_in_min_cells = 3,
                            min_det_feats_per_cell = 5)

##%######################################################%##
#                                                          #
####                   Normalización                    ####
#                                                          #
##%######################################################%##

xenium_gobj <- normalizeGiotto(xenium_gobj)

##%######################################################%##
#                                                          #
####                    Estadísticas                    ####
#                                                          #
##%######################################################%##

xenium_gobj <- addStatistics(xenium_gobj)

##%######################################################%##
#                                                          #
####              Reducción de dimensiones              ####
#                                                          #
##%######################################################%##

# NOTA: en este caso NO calculamos HVFs porque tenemos muy pocos genes.

xenium_gobj <- runPCA(xenium_gobj, 
                      feats_to_use = NULL) 

screePlot(xenium_gobj, 
          ncp = 20)

##%######################################################%##
#                                                          #
####                        UMAP                        ####
#                                                          #
##%######################################################%##

xenium_gobj <- runUMAP(xenium_gobj, 
                       dimensions_to_use = 1:10)

##%######################################################%##
#                                                          #
####                        tSNE                        ####
#                                                          #
##%######################################################%##

xenium_gobj <- runtSNE(xenium_gobj, 
                       dimensions_to_use = 1:10)

##%######################################################%##
#                                                          #
####                  Vecinos cercanos                  ####
#                                                          #
##%######################################################%##

xenium_gobj <- createNearestNetwork(xenium_gobj, 
                                    dimensions_to_use = 1:10, 
                                    k = 10)

##%######################################################%##
#                                                          #
####                 Leiden clustering                  ####
#                                                          #
##%######################################################%##

xenium_gobj <- doLeidenCluster(xenium_gobj,
                               resolution = 0.25, 
                               n_iterations = 100)

spatPlot2D(xenium_gobj,
           plot_method = "scattermore",
           cell_color = "leiden_clus",
           point_size = 0.5,
           point_shape = "no_border",
           background_color = "black"
)

plotPCA(xenium_gobj, 
        cell_color = "leiden_clus",
        point_size = 0.5)

plotUMAP(xenium_gobj, 
         cell_color = "leiden_clus",
         point_size = 0.5)

plotTSNE(xenium_gobj, 
         cell_color = "leiden_clus",
         point_size = 0.5)

##%######################################################%##
#                                                          #
####             Visualización subcelular               ####
#                                                          #
##%######################################################%##

# En esta versión, graficamos el contorno de las células y rellenamos
# usando el Leiden cluster id
spatInSituPlotPoints(xenium_gobj,
                     polygon_feat_type = "cell",
                     polygon_alpha = 1,
                     polygon_line_size = 0.01,
                     polygon_color = "black",
                     polygon_fill = "leiden_clus",
                     polygon_fill_as_factor = TRUE
)

# Visualiza algunos transcritos
# NOTA: este paso solamente funciona si leíste los transcritos al inicio
spatInSituPlotPoints(xenium_gobj,
                     feats = list("rna" = c("LUM", "CXCR4", "ITGAX")),
                     feats_color_code = c(
                         "LUM" = "cyan",
                         "CXCR4" = "magenta",
                         "ITGAX" = "yellow"),
                     point_size = 0.1,
                     xlim = c(1000, 2000),
                     ylim = c(-3000, -2000),
                     plot_last = "polygons",
                     polygon_feat_type = "cell",
                     polygon_alpha = 0.3,
                     polygon_line_size = 0.01,
                     polygon_color = "black",
                     polygon_fill = "leiden_clus",
                     polygon_fill_as_factor = TRUE,
                     show_image = FALSE,
                     image_name = "dapi"
)

##%######################################################%##
#                                                          #
####      Anotación celular asistida por scRNAseq       ####
#                                                          #
##%######################################################%##

# Los datos de Xenium detectan transcritos con una alta sensibilidad,
# sin embargo, la segmentación juega un papel muy importante en la 
# limpieza de las lecturas. Como resultado, algunos (o muchos) transcritos
# pueden llegar a ser anotados en la célula incorrecta.

# En esta paso, usaremos un scRNAseq como referencia para la anotación celular.
# Descarga el dataset de referencia usando el siguiente comando:
download.file(url = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7782nnn/GSM7782696/suppl/GSM7782696%5F5p%5Fcount%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5",
              destfile = "xenium_data/GSM7782696_5p_count_filtered_feature_bc_matrix.h5")

## Pre-procesa el scRNAseq
expr_data_path <- "xenium_data/GSM7782696_5p_count_filtered_feature_bc_matrix.h5"

sc <- createGiottoObject(expression = get10Xmatrix_h5(expr_data_path)$`Gene Expression`)

fids <- featIDs(sc)
mt_genes <- fids[grep("^MT-.*", fids)]

sc <- addFeatsPerc(sc, 
                   expression_values = "raw", 
                   feats = mt_genes, 
                   vector_name = "mito")
ggplot(spatValues(sc, feats = "mito"), aes(x = "cells", y = mito)) + 
    geom_violin()

meta <- pDataDT(sc)
meta[, sum(mito > 50)] # 5006 # filter these out
meta[, sum(mito < 50)] # 17963
sc <- subset(sc, subset = mito < 50)
ggplot(spatValues(sc, feats = "mito"), aes(x = "cells", y = mito)) +
    geom_violin()

## Procesa el objeto
sc <- sc |>
  filterGiotto(
    expression_threshold = 1,
    feat_det_in_min_cells = 20,
    min_det_feats_per_cell = 200
  ) |>
  normalizeGiotto() |>
  calculateHVF() |>
  runPCA()

screePlot(sc, ncp = 50)

sc <- runUMAP(sc, 
              dimensions_to_use = seq(30))
plotUMAP(sc)

sc <- sc |>
  createNearestNetwork(dimensions_to_use = seq(30)) |>
  doLeidenCluster(res = 1.5) # overclustering

plotUMAP(sc, 
          cell_color = "leiden_clus", 
          point_size = 0.6, 
          point_border_stroke = 0
)

# Análisis de comunicación celular


