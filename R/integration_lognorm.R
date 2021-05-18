integration_lognorm <- function(config_csv = config_csv,
                                project_name = project_name,
                                int_regex_mito = int_regex_mito,
                                scale_factor_tab2 = scale_factor_tab2,
                                most_var_method_integration = most_var_method_integration,
                                n_of_var_genes_integration = n_of_var_genes_integration,
                                n_of_PCs_integration = n_of_PCs_integration) {

    sing_cell_list <- list()
    for( i in 1:nrow(config_csv) ) {

        data_10x_raw <- Seurat::Read10X( data.dir = paste0("./data/", config_csv[ i, 1 ] ) )

        data_10x <- Seurat::CreateSeuratObject(counts = data_10x_raw,
                                       project = project_name,
                                       min.cells = as.numeric(config_csv[ i, 3 ]),
                                       min.features = as.numeric(config_csv[ i, 4 ]))

        data_10x <- AddMetaData(data_10x,
                                as.character(config_csv[ i, 2 ]), # name of sample
                                col.name = 'treat')

        data_10x <- Seurat::PercentageFeatureSet(data_10x,
                                         pattern = int_regex_mito,
                                         col.name = "percent.mt")

        data_10x <- base::subset(data_10x,
                           subset = percent.mt < as.numeric(config_csv[ i, 6 ]) &
                               nFeature_RNA < as.numeric(config_csv[ i, 5 ]))

        data_10x <- NormalizeData(data_10x,
                                  verbose = T,
                                  normalization.method = "LogNormalize",
                                  scale.factor = scale_factor_tab2)

        data_10x <- FindVariableFeatures(data_10x,
                                         selection.method = most_var_method_integration,
                                         nfeatures = n_of_var_genes_integration,
                                         verbose = T)

        sing_cell_list[[i]] <- data_10x

    }

    features <- SelectIntegrationFeatures(object.list = sing_cell_list)

    sc_data.anchors <- FindIntegrationAnchors( object.list = sing_cell_list,
                                               anchor.features = features,
                                               )

    sc_data.combined <- IntegrateData( anchorset = sc_data.anchors )

    sc_data.combined

}