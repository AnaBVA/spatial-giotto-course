##%######################################################%##
#                                                          #
####              Descarga el set de datos              ####
#                                                          #
##%######################################################%##

# Usaremos el set de datos de cerebro de ratón de la tecnología Visium,
# disponible en: https://www.10xgenomics.com/datasets/mouse-brain-section-coronal-1-standard-1-1-0

# Para descargar, ejecuta estos comandos en la terminal de R:
dir.create("visium_data/")

download.file(url = "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain_raw_feature_bc_matrix.tar.gz",
              destfile = "visium_data/V1_Adult_Mouse_Brain_raw_feature_bc_matrix.tar.gz")

download.file(url = "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain_spatial.tar.gz",
              destfile = "visium_data/V1_Adult_Mouse_Brain_spatial.tar.gz")

# Descomprime los archivos

untar(tarfile = "visium_data/V1_Adult_Mouse_Brain_raw_feature_bc_matrix.tar.gz", 
      exdir = "visium_data/")

untar(tarfile = "visium_data/V1_Adult_Mouse_Brain_spatial.tar.gz", 
      exdir = "visium_data/")

##%######################################################%##
#                                                          #
####       Configura las instrucciones de Giotto        ####
#                                                          #
##%######################################################%##

library(Giotto)

instructions <- createGiottoInstructions(save_dir = "visium_plots",
                                         save_plot = TRUE,
                                         show_plot = TRUE,
                                         return_plot = FALSE,
                                         python_path = NULL # alternativamente, "/local/python/path/python" 
                                         )

##%######################################################%##
#                                                          #
####               Crea el objeto Giotto                ####
#                                                          #
##%######################################################%##

data_path <- "visium_data/"

visium_brain <- createGiottoVisiumObject(visium_dir = data_path,
                                         expr_data = "raw",
                                         png_name = "tissue_lowres_image.png",
                                         gene_column_index = 2,
                                         instructions = instructions)

##%######################################################%##
#                                                          #
####                 Explora el objeto                  ####
#                                                          #
##%######################################################%##

# Visión general del objeto
visium_brain

# Expresión
showGiottoExpression(visium_brain)

# Metadatos de célula
showGiottoCellMetadata(visium_brain)
pDataDT(visium_brain)

# Metadatos de los genes/proteínas
showGiottoFeatMetadata(visium_brain)

# Imágenes
showGiottoImageNames(visium_brain)

# Reducción de dimensiones
showGiottoDimRed(visium_brain)

# Instrucciones
instructions(visium_brain)

##%######################################################%##
#                                                          #
####          Filtra los spots sobre el tejido          ####
#                                                          #
##%######################################################%##

spatPlot2D(gobject = visium_brain, 
           cell_color = "in_tissue", 
           point_size = 1,
           cell_color_code = c("0" = "lightgrey", "1" = "blue"), 
           show_image = TRUE, 
           image_name = "image")

cell_metadata <- pDataDT(visium_brain)
in_tissue_barcodes <- cell_metadata[in_tissue == 1]$cell_ID

visium_brain <- subsetGiotto(visium_brain, 
                             cell_ids = in_tissue_barcodes)

spatPlot2D(gobject = visium_brain, 
           cell_color = "in_tissue", 
           point_size = 1,
           cell_color_code = c("0" = "lightgrey", "1" = "blue"), 
           show_image = TRUE, 
           image_name = "image")

##%######################################################%##
#                                                          #
####                 Control de calidad                 ####
#                                                          #
##%######################################################%##

# Crea un objecto con el número de genes y transcritos por célula/spot
visium_brain_statistics <- addStatistics(gobject = visium_brain,
                                         expression_values = "raw")

spatPlot2D(gobject = visium_brain_statistics, 
           cell_color = "nr_feats", 
           point_size = 1,
           color_as_factor = FALSE)

# Distribución de genes por spot/célula
filterDistributions(gobject = visium_brain_statistics, 
                    detection = "cells")

# Distribución del número de células con cierto número de genes
filterDistributions(gobject = visium_brain_statistics, 
                    detection = "feats")

# Evalúa diferentes puntos de corte al filtrar
filterCombinations(gobject = visium_brain_statistics,
                   expression_thresholds = c(1, 2, 3),
                   feat_det_in_min_cells = c(50, 100, 200),
                   min_det_feats_per_cell = c(500, 1000, 1500))

##%######################################################%##
#                                                          #
####                      Filtrado                      ####
#                                                          #
##%######################################################%##

visium_brain <- filterGiotto(gobject = visium_brain,
                             expression_threshold = 1,
                             feat_det_in_min_cells = 50,
                             min_det_feats_per_cell = 1000,
                             expression_values = "raw",
                             verbose = TRUE)

##%######################################################%##
#                                                          #
####                   Normalización                    ####
#                                                          #
##%######################################################%##

visium_brain <- normalizeGiotto(gobject = visium_brain, 
                                scalefactor = 6000, 
                                verbose = TRUE)

##%######################################################%##
#                                                          #
####                    Estadísticas                    ####
#                                                          #
##%######################################################%##

visium_brain <- addStatistics(gobject = visium_brain)

spatPlot2D(gobject = visium_brain, 
           show_image = TRUE, 
           point_alpha = 0.7,
           cell_color = "nr_feats", 
           color_as_factor = FALSE)

##%######################################################%##
#                                                          #
####                        HVF                         ####
#                                                          #
##%######################################################%##

# Regresión de Loess es usada cuando la relación entre la expresión promedio
# y la varianza son no lineales o pueden ser descritas por un modelo no paramétrico

visium_brain <- calculateHVF(gobject = visium_brain, 
                             method = "cov_loess",
                             save_plot = TRUE)

# Residuales Pearson son usados para la estabilización de la varianza 
# (para lidiar con el ruido técnico) y resaltar los genes con alta dispersión

visium_brain <- calculateHVF(gobject = visium_brain, 
                             method = "var_p_resid",
                             save_plot = TRUE)

# Grupos de covarianza (default), son usagos cuando la variabilidad de la expresión 
# difiere a través de los niveles de expresión o las regiones espaciales, sin
# asumir una relación específica entre la expresión promedio y la varianza.

visium_brain <- calculateHVF(gobject = visium_brain, 
                             method = "cov_groups",
                             save_plot = TRUE)

##%######################################################%##
#                                                          #
####              Reducción de dimensiones              ####
#                                                          #
##%######################################################%##

visium_brain <- runPCA(gobject = visium_brain)

screePlot(gobject = visium_brain, 
          ncp = 30)

plotPCA(gobject = visium_brain)

##%######################################################%##
#                                                          #
####                        UMAP                        ####
#                                                          #
##%######################################################%##

visium_brain <- runUMAP(visium_brain, 
                        dimensions_to_use = 1:10)

plotUMAP(gobject = visium_brain)

##%######################################################%##
#                                                          #
####                        tSNE                        ####
#                                                          #
##%######################################################%##

visium_brain <- runtSNE(visium_brain, 
                        dimensions_to_use = 1:10)

plotTSNE(gobject = visium_brain)

##%######################################################%##
#                                                          #
####                  Vecinos cercanos                  ####
#                                                          #
##%######################################################%##

# sNN network (default)
visium_brain <- createNearestNetwork(gobject = visium_brain, 
                                     dimensions_to_use = 1:10, 
                                     k = 15)

# K-Nearest Neighbors (kNN)
visium_brain <- createNearestNetwork(gobject = visium_brain, 
                                     dimensions_to_use = 1:10, 
                                     k = 15,
                                     type = "kNN")

##%######################################################%##
#                                                          #
####                 Leiden clustering                  ####
#                                                          #
##%######################################################%##

visium_brain <- doLeidenCluster(gobject = visium_brain, 
                                resolution = 0.4, 
                                n_iterations = 1000)

plotPCA(gobject = visium_brain,
        cell_color = "leiden_clus")

plotUMAP(gobject = visium_brain,
         cell_color = "leiden_clus", 
         show_NN_network = TRUE, 
         point_size = 2.5)

plotTSNE(gobject = visium_brain,
         cell_color = "leiden_clus", 
         show_NN_network = TRUE, 
         point_size = 2.5)

spatPlot2D(visium_brain, 
           cell_color = "leiden_clus",
           point_size = 1,
           coord_fix_ratio = 1)

# Podemos resaltar un grupo de interés
spatPlot2D(visium_brain, 
           cell_color = "leiden_clus",
           select_cell_groups = "8", 
           point_size = 1,
           coord_fix_ratio = 1, 
           show_other_cells = TRUE,
           cell_color_code = c("8" = "red"), 
           other_cell_color = "grey")

##%######################################################%##
#                                                          #
####              Integración con scRNAseq              ####
#                                                          #
##%######################################################%##

## Descarga el dataset de referencia
GiottoData::getSpatialDataset(dataset = "scRNA_mouse_brain", 
                              directory = data_path)

## Lee el scRNAseq como objeto Giotto
sc_expression <- file.path(data_path, "brain_sc_expression_matrix.txt.gz")
sc_metadata <- file.path(data_path, "brain_sc_metadata.csv")

g_sc <- createGiottoObject(expression = sc_expression,
                           instructions = instructions)

## Agrega los metadatos
cell_metadata <- data.table::fread(sc_metadata)

g_sc <- addCellMetadata(g_sc, 
                        new_metadata = cell_metadata[, c("Class", "Subclass")])

## Normalización
g_sc <- normalizeGiotto(g_sc)

##%######################################################%##
#                                                          #
####         Enriquecimiento de tipos celulares         ####
#                                                          #
##%######################################################%##

## Encuentra los genes que caracterizan el comportamiento de cada tipo celular 
## en el scRNAseq de referencia

markers_scran <- findMarkers_one_vs_all(gobject = g_sc, 
                                        method = "scran",
                                        expression_values = "normalized",
                                        cluster_column = "Class", 
                                        min_feats = 3)

topgenes_scran <- markers_scran[, head(.SD, 10), by = "cluster"]

celltypes <- levels(factor(markers_scran$cluster)) 

sign_list <- list()
for (i in 1:length(celltypes)){
    sign_list[[i]] = topgenes_scran[which(topgenes_scran$cluster == celltypes[i]),]$feats
}

sign_matrix <- makeSignMatrixPAGE(sign_names = celltypes,
                                  sign_list = sign_list)

visium_brain <- runPAGEEnrich(gobject = visium_brain, 
                              sign_matrix = sign_matrix)

spatCellPlot2D(gobject = visium_brain,
               spat_enr_names = "PAGE",
               cell_annotation_values = celltypes,
               cow_n_col = 3,
               coord_fix_ratio = 1, 
               point_size = 1, 
               show_legend = TRUE)

##%######################################################%##
#                                                          #
####                   Deconvolución                    ####
#                                                          #
##%######################################################%##

library(spacexr)
library(SummarizedExperiment)
library(SpatialExperiment)

# Crea el objeto referencia desde el scRNAseq 
sc_counts <- getExpression(g_sc,
                           values = "raw",
                           output = "matrix")

cell_metadata <- pDataDT(g_sc)
sc_cell_types <- factor(cell_metadata$Class)
names(sc_cell_types) <- cell_metadata$cell_ID

reference_se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = sc_counts),
  colData = data.frame(cell_type = sc_cell_types)
)

# Crea el objeto espacial desde el spatial transcriptomics
spatial_locs <- getSpatialLocations(visium_brain,
                                    output = "data.table")

st_coords <- data.frame(x = spatial_locs$sdimx,
                        y = spatial_locs$sdimy)
rownames(st_coords) <- spatial_locs$cell_ID

st_counts <- getExpression(visium_brain,
                           values = "raw",
                           output = "matrix")

spatial_spe <- SpatialExperiment::SpatialExperiment(
  assay = st_counts,
  spatialCoords = as.matrix(st_coords)
)

# Crea y ejecuta RCTD
myRCTD <- createRctd(spatial_experiment = spatial_spe, 
                     reference_experiment = reference_se
)

results_spe <- runRctd(myRCTD, rctd_mode = "full", max_cores = 1)

# Da formato a los resultados
results_full <- results_spe@assays@data$weights

x <- as.data.frame(t(as.matrix(results_full)))
norm_weights <- sweep(x, 1, rowSums(x), '/') 
norm_weights <- cbind(data.frame(cell_ID = rownames(results_spe@colData)),
                      norm_weights)


# Visualiza el resultado
library(ggplot2)
library(scatterpie)

rctd_proportions <- merge(spatial_locs, norm_weights)

ggplot() +
  geom_scatterpie(aes(sdimx, sdimy),
                  data = rctd_proportions,
                  cols = colnames(x),
                  pie_scale = 0.5,
                  linewidth = 0) +
  labs(fill = "") +
  theme_minimal()
ggsave("visium_results/deconvolution_plot.png")

##%######################################################%##
#                                                          #
####           Genes con patrones espaciales            ####
#                                                          #
##%######################################################%##

# createSpatialNetwork() crea una red espacial que conecta a las células/spots
# basados en su distancia física. Existen dos métodos: Delaunay y kNN.
# En el método Delaunay, los vecinos son elegidos por la triangulación Delaunay
# y un criterio de máxima distancia.
# En el método kNN, el número de vecinos puede ser determinado por el parámetro k,
# o por la máxima distancia entre células.

visium_brain <- createSpatialNetwork(gobject = visium_brain,
                                     method = "kNN", 
                                     k = 6,
                                     maximum_distance_knn = 400,
                                     name = "spatial_network")

spatPlot2D(gobject = visium_brain,  
           show_network= TRUE,
           network_color = "blue", 
           spatial_network_name = "spatial_network")

# Usa the el método de ranking binario para ordenar los genes que muestran 
# un patrón espacial. Este paso puede tomar algunos minutos.

ranktest <- binSpect(visium_brain, 
                     bin_method = "rank",
                     calc_hub = TRUE, 
                     hub_min_int = 5,
                     spatial_network_name = "spatial_network")

spatFeatPlot2D(visium_brain, 
               expression_values = "scaled",
               feats = ranktest$feats[1:4], 
               cow_n_col = 2, 
               point_size = 1)

##%######################################################%##
#                                                          #
####              Módulos de co-expresión               ####
#                                                          #
##%######################################################%##

# Encuentra los genes con patrones espaciales similares y crea módulos de co-expresión

# Toma los primeros 500 genes con patrón espacial
ext_spatial_genes <- ranktest[1:500,]$feats

# Calcula las distacias entre los genes
spat_cor_netw_DT <- detectSpatialCorFeats(
  visium_brain,
  method = "network",
  spatial_network_name = "spatial_network",
  subset_feats = ext_spatial_genes)

# Identifica los genes con mayor correlación a un cierto gene
top10_genes <- showSpatialCorFeats(spat_cor_netw_DT, 
                                   feats = "Mbp", 
                                   show_top_feats = 10)

spatFeatPlot2D(visium_brain, 
               expression_values = "scaled",
               feats = as.character(top10_genes$variable[1:4]), 
               point_size = 1.5)

# Agrupa los genes con mayor correlación, define el número de clusters usando 
# el parámetro k 
spat_cor_netw_DT <- clusterSpatialCorFeats(spat_cor_netw_DT, 
                                           name = "spat_netw_clus", 
                                           k = 20)

# Visualiza la correlación de los genes
heatmSpatialCorFeats(visium_brain,
                     spatCorObject = spat_cor_netw_DT,
                     use_clus_name = "spat_netw_clus",
                     heatmap_legend_param = list(title = NULL))

# Ordena los clusters
netw_ranks <- rankSpatialCorGroups(
  visium_brain,
  spatCorObject = spat_cor_netw_DT, 
  use_clus_name = "spat_netw_clus")

# Calcula el score de enriquecimiento de los metagenes por cluster
cluster_genes_DT <- showSpatialCorFeats(spat_cor_netw_DT, 
                                        use_clus_name = "spat_netw_clus", 
                                        show_top_feats = 1)

cluster_genes <- cluster_genes_DT$clus 
names(cluster_genes) <- cluster_genes_DT$feat_ID

visium_brain <- createMetafeats(visium_brain, 
                                feat_clusters = cluster_genes, 
                                name = "cluster_metagene")

spatCellPlot(visium_brain,
             spat_enr_names = "cluster_metagene",
             cell_annotation_values = netw_ranks$clusters,
             point_size = 1, 
             cow_n_col = 5)

##%######################################################%##
#                                                          #
####              Calcula Leiden clusters               ####
####           usando la información espacial           ####
#                                                          #
##%######################################################%##

# Toma los primeros genes por cada módulo de co-expresión
coexpr_dt <- data.table::data.table(
  genes = names(spat_cor_netw_DT$cor_clusters$spat_netw_clus),
  cluster = spat_cor_netw_DT$cor_clusters$spat_netw_clus)

data.table::setorder(coexpr_dt, cluster)

top30_coexpr_dt <- coexpr_dt[, head(.SD, 30) , by = cluster]

spatial_genes <- top30_coexpr_dt$genes

##%######################################################%##
#                                                          #
####      Re-calcula la dimensión de reducciones y      ####
####      los clusters usando los genes espaciales      ####
#                                                          #
##%######################################################%##

visium_brain <- runPCA(gobject = visium_brain,
                       feats_to_use = spatial_genes,
                       name = "custom_pca")

visium_brain <- runUMAP(visium_brain, 
                        dim_reduction_name = "custom_pca", 
                        dimensions_to_use = 1:20,
                        name = "custom_umap")

visium_brain <- createNearestNetwork(gobject = visium_brain,
                                     dim_reduction_name = "custom_pca",
                                     dimensions_to_use = 1:20, 
                                     k = 5,
                                     name = "custom_NN")

visium_brain <- doLeidenCluster(gobject = visium_brain, 
                                network_name = "custom_NN",
                                resolution = 0.15, 
                                n_iterations = 1000,
                                name = "custom_leiden")

spatPlot2D(visium_brain, 
           cell_color = "custom_leiden",
           point_size = 3)
