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
# Descarga el dataset de referencia y léelo usando el siguiente comando:
sc <- loadGiotto("sc")

plotUMAP(sc, cell_color = "cell_types")


# Crea el objeto conjunto de Xenium + scRNAseq

## Identifica los genes en común
ufids <- intersect(featIDs(xenium_gobj), featIDs(sc)) 

## Remueve los elementos no comunes entre los objetos
sc_join <- giotto() |> 
  setGiotto(sc[[c("expression", "spatial_locs"), "raw"]]) 
xen_join <- giotto() |>
  setGiotto(xenium_gobj[[c("expression", "spatial_locs"), "raw", spat_unit = "cell"]])

sc_join <- addCellMetadata(sc_join,
                           new_metadata = pDataDT(sc)[,c("cell_ID", "cell_types")],
                           by_column = TRUE)

j <- joinGiottoObjects(
  list(sc_join[ufids], xen_join[ufids]), 
  gobject_names = c("sc", "xen"))

instructions(j) <- instructions

## Procesa el objeto
j <- filterGiotto(j,
                  expression_threshold = 1, 
                  feat_det_in_min_cells = 1, 
                  min_det_feats_per_cell = 10)

j <- normalizeGiotto(j) # pearson and quantile norms do not work well for this

j <- runPCA(j, 
            feats_to_use = NULL)

plotPCA(j,
        cell_color = "list_ID",
        point_border_stroke = 0,
        point_size = 0.3)

screePlot(j, 
          ncp = 50)

j <- runUMAP(j, 
             dimensions_to_use = seq(10))

plotUMAP(j, 
         cell_color = "list_ID",
         point_size = 0.3, 
         point_border_stroke = 0, 
         cell_color_code = c("red", "blue"))

## Procesa el objeto removiendo el efecto de batch
j <- runGiottoHarmony(j,
                      vars_use = "list_ID",
                      dim_reduction_name = "pca",
                      dimensions_to_use = seq(10),
                      name = "harmony")

j <- runUMAP(j, 
             dim_reduction_to_use = "harmony",
             dim_reduction_name = "harmony",
             name = "umap_harmony",
             dimensions_to_use = seq(10))

plotUMAP(j, 
         cell_color = "cell_types",
         dim_reduction_name = "umap_harmony", 
         cell_color_code = c(
           "turquoise", "azure4", "pink", "tomato", 
           "chartreuse3", "lightslateblue", "firebrick4", 
           "blue", "magenta", "cornflowerblue"),
         point_size = 0.5, 
         select_cells = spatIDs(j, subset = list_ID == "sc"),
         point_border_stroke = 0, 
         other_point_size = 0.3)

## Transferencia de etiquetas
j <- labelTransfer(j, 
                   source_cell_ids = spatIDs(j, subset = list_ID == "sc"),
                   k = 10,
                   labels = "cell_types",
                   reduction_method = "harmony",
                   reduction_name = "harmony",
                   dimensions_to_use = seq(10))

plotUMAP(j, 
         cell_color = "trnsfr_cell_types", 
         dim_reduction_name = "umap_harmony", 
         cell_color_code = c(
           "turquoise", "azure4", "pink", "tomato", 
           "chartreuse3", "lightslateblue", "firebrick4", 
           "blue", "magenta", "cornflowerblue"),
         point_size = 0.3,
         point_border_stroke = 0)

## Agrega la anotación en el Xenium original
ann_meta <- pDataDT(j)
ann_meta <- ann_meta[list_ID == "xen",]
ann_meta[, cell_ID := gsub("^xen-", "", cell_ID)]
ann_meta <- ann_meta[, .(cell_ID, trnsfr_cell_types, trnsfr_cell_types_prob)]

xenium_gobj <- addCellMetadata(xenium_gobj, 
                               new_metadata = ann_meta, 
                               by_column = TRUE, 
                               column_cell_ID = "cell_ID")

spatInSituPlotPoints(xenium_gobj, 
                     polygon_fill = "trnsfr_cell_types", 
                     polygon_fill_as_factor = TRUE,
                     polygon_line_size = 0,
                     polygon_fill_code = c(
                       "turquoise", "azure4", "pink", "tomato", 
                       "chartreuse3", "lightslateblue", "firebrick4", 
                       "blue", "magenta", "cornflowerblue"))

