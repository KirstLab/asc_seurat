integration_sctransform <- function(config_csv = config_csv,
                                    project_name = project_name,
                                    int_regex_mito = int_regex_mito) {

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

        sing_cell_list[[i]] <- data_10x

    }

    sing_cell_list <- lapply(X = sing_cell_list,
                             FUN = SCTransform,
                             method = "glmGamPoi",
                             vars.to.regress = "percent.mt")
    
    features <- SelectIntegrationFeatures(object.list = sing_cell_list,
                                          nfeatures = 3000)

    sing_cell_list <- PrepSCTIntegration(object.list = sing_cell_list,
                                         anchor.features = features)

    sc_anchors <- FindIntegrationAnchors(object.list = sing_cell_list,
                                         normalization.method = "SCT",
                                         anchor.features = features)

    sc_combined.sct <- IntegrateData(anchorset = sc_anchors,
                                     normalization.method = "SCT")

    return(sc_combined.sct)

}