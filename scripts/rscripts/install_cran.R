#####################################################
#       Packages that need to be installed          #
#####################################################
# From CRAN
packages_cran <- c("tidyverse",
                   "Seurat",
                   "SeuratObject",
                   "patchwork",
                   "ggplot2",
                   "circlize",
                   "reactable",
                   "sctransform",
                   "shiny",
                   "shinyWidgets",
                   "shinyFeedback",
                   "rclipboard",
                   "future",
                   "ggthemes",
                   "metap",
                   "shinycssloaders",
                   "DT",
                   "dplyr",
                   "hdf5r")

for (i in packages_cran){
    install.packages(i, dep = TRUE)
}
