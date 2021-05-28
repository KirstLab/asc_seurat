lognorm_function <- function(name,
                             sc_data = sc_data,
                             to_filter = to_filter,
                             min_count = min_count,
                             max_count = max_count,
                             max_mito_perc = max_mito_perc,
                             scale_factor = scale_factor,
                             filter_clusters = filter_clusters,
                             filter_clusters_opt = filter_clusters_opt,
                             most_var_method = most_var_method,
                             n_of_var_genes = n_of_var_genes) {

    showNotification("Normalizing the data (LogNormalize)",
                     duration = NULL,
                     id = "m3")

    sc_data <- base::subset(sc_data,
                      subset = nFeature_RNA > min_count &
                          nFeature_RNA < max_count &
                          percent.mt < max_mito_perc)

    if ( filter_clusters == 1 ) {

        if (filter_clusters_opt == 0) { # select

            ret_data <- base::subset(sc_data,
                               cells = to_filter)

        } else if (filter_clusters_opt == 1) { # exclude

            ret_data <- base::subset(sc_data,
                               cells = to_filter)

        }

    } else {
        ret_data <- sc_data
    }

    ret_data <- NormalizeData(ret_data,
                              normalization.method = "LogNormalize",
                              scale.factor = scale_factor)

    on.exit(removeNotification(id = "m3"), add = TRUE)

    showNotification("Scalling the data",
                     duration = NULL,
                     id = "m4")

    ret_data <- FindVariableFeatures(ret_data,
                                     selection.method = most_var_method,
                                     nfeatures = n_of_var_genes)

    allgenes <- base::rownames(ret_data)
    ret_data <- ScaleData(ret_data, features = allgenes)

    on.exit(removeNotification(id = "m4"), add = TRUE)

    showNotification("Running PCA",
                     duration = NULL,
                     id = "m5")

    ret_data <-  RunPCA(ret_data,
                        features = VariableFeatures(object = ret_data), verbose = F)

    on.exit(removeNotification(id = "m5"), add = TRUE)

    return(ret_data)
}
