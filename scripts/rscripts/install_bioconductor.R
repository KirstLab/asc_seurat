#####################################################
#       Packages that need to be installed          #
#####################################################
# From Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

packages_bioc <- c("ComplexHeatmap",
                   "tradeSeq",
                   "SingleCellExperiment",
                   "slingshot",
                   "multtest",
                   "biomaRt",
                   "topGO")

for (i in packages_bioc){
    BiocManager::install(i)
}
