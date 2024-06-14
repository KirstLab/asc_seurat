require(scCustomize)
require(cowplot)

pbmc.data <- Read10X("data/example_PBMC/")

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

sing_cell_data <- pbmc

sc_meta <- as.data.frame(sing_cell_data[[]])
sc_meta$cellcluster <- base::rownames(sc_meta)
sc_meta <- sc_meta[, c("cellcluster", "seurat_clusters")]

table_n <- as.data.frame(base::table(sc_meta$seurat_clusters))
table_n <- table_n %>%
    dplyr::rename(Cluster = "Var1",
                  `N. of cells` = "Freq")


table_n <- table_n %>%
    dplyr::rename(Cluster = "Var1",
                  `N. of cells` = "Freq")
# Assuming your table is stored as a data frame named "table_n"
ggplot(table_n, aes(x = Cluster, y = `N. of cells`)) +
    geom_bar(stat = "identity", fill = "#008ac2ff") +
    labs(x = "Cluster", y = "Number of Cells", title = "Number of Cells per Cluster") +
    theme_minimal() +
    geom_text(aes(label = `N. of cells`), vjust = -0.5, size = 7) +
    theme(
        axis.title = element_text(size = 20),   # Axis titles
        axis.text = element_text(size = 16),    # Axis tick labels
        plot.title = element_blank(),   # Plot title
        legend.text = element_text(size = 12),  # Legend text (if you have a legend)
        legend.title = element_text(size = 14)  # Legend title (if you have a legend)
    )

