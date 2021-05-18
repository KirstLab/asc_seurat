SCTransform_function <- function(name,
                                 sc_data = sc_data,
                                 to_filter = to_filter,
                                 min_count = min_count,
                                 max_count = max_count,
                                 max_mito_perc = max_mito_perc,
                                 filter_clusters = filter_clusters,
                                 filter_clusters_opt = filter_clusters_opt) {

    showNotification("Normalizing the data (SCTransform)",
                     duration = NULL,
                     id = "m3")

    sc_data <- subset(sc_data,
                      subset = nFeature_RNA > min_count &
                          nFeature_RNA < max_count &
                          percent.mt < max_mito_perc)

    if ( filter_clusters == 1 ) {

        if (filter_clusters_opt == 0) { # select

            ret_data <- subset(sc_data,
                               cells = to_filter)

        } else if (filter_clusters_opt == 1) { # exclude

            ret_data <- subset(sc_data,
                               cells = to_filter)

        }

    }     else {
        ret_data <- sc_data
    }

    ret_data <- SCTransform(ret_data, vars.to.regress = "percent.mt", verbose = FALSE)

    on.exit(removeNotification(id = "m3"), add = TRUE)

    showNotification("Running PCA",
                     duration = NULL,
                     id = "m5")

    ret_data <-  RunPCA(ret_data,
                        verbose = F)
    
    on.exit(removeNotification(id = "m5"), add = TRUE)

    return(ret_data)

}
