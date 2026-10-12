# Instala Giotto
remotes::install_github("giotto-suite/GiottoUtils")
remotes::install_github("giotto-suite/GiottoClass")
remotes::install_github("giotto-suite/GiottoVisuals")
remotes::install_github("giotto-suite/Giotto")

# Instala el ambiente de Giotto
Giotto::installGiottoEnvironment()

# Revisa que el ambiente se ha instalado, debe dar TRUE
Giotto::checkGiottoEnvironment()

# Instala otros paquetes para el pipeline
BiocManager::install("spacexr")
BiocManager::install("SummarizedExperiment")
BiocManager::install("SpatialExperiment")
install.packages("ggplot2")
install.packages("scatterpie")
