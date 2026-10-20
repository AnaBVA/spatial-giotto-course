library(Giotto)

# Carga un objeto visium y un vizgen pequeño desde el paquete GiottoData
mini_visium <- GiottoData::loadGiottoMini("visium")
mini_vizgen <- GiottoData::loadGiottoMini("vizgen")

# Imprime el objeto
mini_visium
mini_vizgen

# Lista la información disponible
list_cell_metadata(mini_visium)
list_cell_metadata(mini_vizgen)

list_expression(mini_visium)
list_expression(mini_vizgen)

list_expression_names(vizmini, 
                      spat_unit = "z1", 
                      feat_type = "rna")

list_dim_reductions(mini_visium)
list_dim_reductions(mini_vizgen)

# Da un vistazo a los componentes del objeto
showGiottoCellMetadata(mini_visium)
showGiottoCellMetadata(mini_vizgen)

showGiottoExpression(mini_visium)

showGiottoImageNames(mini_visium)

# Extrae la información del objeto
getExpression(mini_visium)
getExpression(mini_visium,
              values = "normalized")

getExpression(mini_vizgen,
              spat_unit = "aggregate")
getExpression(mini_vizgen,
              spat_unit = "aggregate",
              output = "matrix")

getDimReduction(mini_visium,
                reduction_method = "pca")
getDimReduction(mini_visium,
                reduction_method = "umap",
                output = "data.table")








