# Asc-Seurat
# Version 2.1
set.seed(1407)
options(shiny.sanitize.errors = FALSE)

# CRAN
suppressMessages( require(tidyverse) )
suppressMessages( require(Seurat) )
suppressMessages( require(SeuratObject) )
suppressMessages( require(patchwork) )
suppressMessages( require(ggplot2) )
suppressMessages( require(circlize) )
suppressMessages( require(reactable) )
suppressMessages( require(sctransform) )
suppressMessages( require(shiny) )
suppressMessages( require(shinyWidgets) )
suppressMessages( require(shinyFeedback) )
suppressMessages( require(rclipboard) )
suppressMessages( require(future) )
suppressMessages( require(ggthemes) )
suppressMessages( require(metap) )
suppressMessages( require(shinycssloaders) )
suppressMessages( require(DT) )
suppressMessages( require(dplyr) )
suppressMessages( require(hdf5r) )
suppressMessages( require(scales) )
suppressMessages( require(utils) )
suppressMessages( require(vroom) )

# Bioconductor
suppressMessages( require(ComplexHeatmap) )
suppressMessages( require(SingleCellExperiment) )
suppressMessages( require(slingshot) )
suppressMessages( require(multtest) )
suppressMessages( require(biomaRt) )
suppressMessages( require(topGO) )

# dynverse packages
suppressMessages( require(dynplot) )
suppressMessages( require(dynwrap) )
suppressMessages( require(dynfeature) )

## New packages in version 2.1
suppressMessages( require(glmGamPoi) ) # Bioconductor
##

## Allows parallelization of some of Seurat's functions
#plan("multicore")

if (dir.exists('/app/user_work')) {
    
    source("/app/R/feature_plot.R")
    source("/app/R/adapted_slingshot.R")
    source("/app/R/biomart_section.R")
    source("/app/R/lognormalization_function.R")
    source("/app/R/SCTransform_function.R")
    source("/app/R/server_functions.R")
    source("/app/R/functions_load_markers_list.R")
    source("/app/R/finding_markers.R")
    source("/app/R/improved_dot_and_violin_plots.R")
    source("/app/R/integration_lognorm.R")
    source("/app/R/integration_SCTransform.R")
    
} else {
    source("R/feature_plot.R")
    source("R/adapted_slingshot.R")
    source("R/biomart_section.R")
    source("R/lognormalization_function.R")
    source("R/SCTransform_function.R")
    source("R/server_functions.R")
    source("R/functions_load_markers_list.R")
    source("R/finding_markers.R")
    source("R/improved_dot_and_violin_plots.R")
    source("R/integration_lognorm.R")
    source("R/integration_SCTransform.R") 
}

###########################
### Load phytozome mart ###
###########################
phytozome_mart <- new("Mart",
                      biomart = "phytozome_mart",
                      vschema = "zome_mart",
                      host    = "https://phytozome.jgi.doe.gov:443/biomart/martservice")

function(input, output, session) {
    
    if (dir.exists('/app/user_work')) {
        setwd('/app/user_work')
    }
    
    #####################################
    ######   Tab 1 - Clustering    ######
    #####################################
    
    observeEvent(input$min_features, {
        
        updateNumericInput(inputId = "min_count", value = input$min_features)
        
    })
    
    observeEvent(input$min_count, {
        
        updateNumericInput(inputId = "max_count", value = (input$min_count + 2000))
        
    })
    
    output$select_sample_tab1 = renderUI({
        
        dir_list <- list.dirs('./data', recursive=FALSE)
        
        div(class = "option-group",
            shinyWidgets::pickerInput(
                inputId = "sample_folder_tab1",
                label = "Select the sample to use",
                choices = sort(dir_list),
                multiple = FALSE,
                options = list(`actions-box` = TRUE)
            ))
    })
    
    single_cell_data_reac <- eventReactive(input$load_10X, {
        
        shinyFeedback::feedbackWarning("min_cells", is.na(input$min_cells), "Required value")
        shinyFeedback::feedbackWarning("min_features", is.na(input$min_features), "Required value")
        shinyFeedback::feedbackWarning("mito_regex", !shiny::isTruthy(input$mito_regex), "Required value")
        
        req(input$min_cells)
        req(input$min_features)
        req(input$mito_regex)
        
        showNotification("Loading the data",
                         duration = NULL,
                         id = "p1")
        
        sing_cell_data.data <- Seurat::Read10X(data.dir = req(input$sample_folder_tab1))
        
        # Initialize the Seurat object with the raw (non-normalized data).
        sing_cell_data <- Seurat::CreateSeuratObject(counts = sing_cell_data.data,
                                                     project = input$proj_name,
                                                     min.cells = input$min_cells,
                                                     min.features = input$min_features)
        
        # Calculate the % of mithocondrial contamination
        sing_cell_data <- Seurat::PercentageFeatureSet(sing_cell_data,
                                                       pattern = input$mito_regex,
                                                       col.name = "percent.mt")
        
        return(sing_cell_data)
        
    })
    
    output$VlnPlot <- renderPlot({
        
        data_set <- req( single_cell_data_reac() )
        on.exit(removeNotification(id = "p1"), add = TRUE)
        
        Seurat::VlnPlot(data_set,
                        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                        ncol = 3,
                        split.plot = F)
        
    })
    
    observeEvent(input$run_vinplot, {
        
        data_sc <- req( single_cell_data_reac() )
        
        test_cond <- if( !is.na(input$max_count) && !is.na(input$min_count) ) {
            
            shinyFeedback::feedbackWarning("max_count",
                                           input$min_count > input$max_count,
                                           "No cells will be selected by appling this parameters!")
            
            shinyFeedback::feedbackWarning("min_count",
                                           input$min_count > input$max_count,
                                           "No cells will be selected by appling this parameters!")
            
            validate(need(input$min_count < input$max_count, "Error: No cells will be selected by appling this parameters!"))
            
        }
        
        if ( !is.na(input$min_count) ) {
            
            data_sc <- base::subset(data_sc,
                                    subset = nFeature_RNA > input$min_count)
            
        }
        
        if ( !is.na(input$max_count) ) {
            
            shinyFeedback::feedbackWarning("max_count",
                                           input$max_count <= 0,
                                           "No cells will be selected by appling this parameters!")
            
            validate(need(input$max_count > 0, "Error: No cells will be selected by appling this parameters!"))
            
            data_sc <- base::subset(data_sc,
                                    subset =  nFeature_RNA < input$max_count)
            
        }
        
        if ( !is.na(input$max_mito_perc) ) {
            
            shinyFeedback::feedbackWarning("max_mito_perc",
                                           input$max_mito_perc <= 0,
                                           "No cells will be selected by appling this parameters!")
            
            validate(need(input$max_mito_perc > 0, "Error: No cells will be selected by appling this parameters!"))
            
            data_sc <- base::subset(data_sc,
                                    subset =  percent.mt < input$max_mito_perc)
            
        }
        
        single_cell_data_filt <- data_sc
        
        output$VlnPlot_filt <- renderPlot({
            
            Seurat::VlnPlot(single_cell_data_filt,
                            features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                            ncol = 3,
                            split.plot = F)
            
        })
        
    })
    
    output$p1_down <- downloadHandler(
        
        filename = function() {
            paste("Violin_plot", ".", input$p1_format, sep = "")
        },
        content = function(file) {
            
            withProgress(message = "Please wait, preparing the data for download.",
                         value = 0.5, {
                             
                             data_sc <- req( single_cell_data_reac() )
                             
                             test_cond <- if( !is.na(input$max_count) && !is.na(input$min_count) ) {
                                 
                                 shinyFeedback::feedbackWarning("max_count",
                                                                input$min_count > input$max_count,
                                                                "No cells will be selected by appling this parameters!")
                                 
                                 shinyFeedback::feedbackWarning("min_count",
                                                                input$min_count > input$max_count,
                                                                "No cells will be selected by appling this parameters!")
                                 
                                 validate(need(input$min_count < input$max_count, "Error: No cells will be selected by appling this parameters!"))
                                 
                             }
                             
                             if ( !is.na(input$min_count) ) {
                                 
                                 data_sc <- base::subset(data_sc,
                                                         subset = nFeature_RNA > input$min_count)
                                 
                             }
                             
                             if ( !is.na(input$max_count) ) {
                                 
                                 shinyFeedback::feedbackWarning("max_count",
                                                                input$max_count <= 0,
                                                                "No cells will be selected by appling this parameters!")
                                 
                                 validate(need(input$max_count > 0, "Error: No cells will be selected by appling this parameters!"))
                                 
                                 data_sc <- base::subset(data_sc,
                                                         subset =  nFeature_RNA < input$max_count)
                                 
                                 
                                 
                             }
                             
                             if ( !is.na(input$max_mito_perc) ) {
                                 
                                 shinyFeedback::feedbackWarning("max_mito_perc",
                                                                input$max_mito_perc <= 0,
                                                                "No cells will be selected by appling this parameters!")
                                 
                                 validate(need(input$max_mito_perc > 0, "Error: No cells will be selected by appling this parameters!"))
                                 
                                 data_sc <- base::subset(data_sc,
                                                         subset =  percent.mt < input$max_mito_perc)
                                 
                             }
                             
                             height <- as.numeric( req( input$p1_height) )
                             width <- as.numeric( req( input$p1_width) )
                             res <- as.numeric( req( input$p1_res) )
                             
                             p <- Seurat::VlnPlot(data_sc,
                                                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                                                  ncol = 3,
                                                  split.plot = F)
                             
                             ggplot2::ggsave(file,
                                             p,
                                             height=height,
                                             width=width,
                                             units="cm",
                                             dpi=res)
                             
                         })
            
        }
    )
    
    single_cell_data_pca <- eventReactive( c(input$run_pca, input$rerun_after_filtering), {
        
        sc_data <- req( single_cell_data_reac() )
        
        if (input$filter_clusters == 1) { to_filter <- req( to_filter() )}
        
        req(input$min_count)
        req(input$max_count)
        req(input$max_mito_perc)
        req(input$filter_clusters)
        req(input$filter_clusters_opt)
        req(input$most_var_method)
        
        if (input$normaliz_method == 0) { #LogNormalize
            
            shinyFeedback::feedbackWarning("scale_factor", is.na(input$scale_factor), "Required value")
            shinyFeedback::feedbackWarning("n_of_var_genes", is.na(input$n_of_var_genes), "Required value")
            
            req(input$scale_factor)
            req(input$n_of_var_genes)
            
            sc_data2 <- lognorm_function("Test",
                                         to_filter = to_filter,
                                         sc_data = sc_data,
                                         min_count = input$min_count,
                                         max_count = input$max_count,
                                         max_mito_perc = input$max_mito_perc,
                                         scale_factor = input$scale_factor,
                                         filter_clusters = input$filter_clusters,
                                         filter_clusters_opt = input$filter_clusters_opt,
                                         most_var_method = input$most_var_method,
                                         n_of_var_genes = input$n_of_var_genes)
            
        } else if (input$normaliz_method == 1) { #"SCTransform"
            
            sc_data2 <- SCTransform_function("Test",
                                             to_filter = to_filter,
                                             sc_data = sc_data,
                                             min_count = input$min_count,
                                             max_count = input$max_count,
                                             max_mito_perc = input$max_mito_perc,
                                             filter_clusters = input$filter_clusters,
                                             filter_clusters_opt = input$filter_clusters_opt)
            
            # This won't be used until the DE analysis and visualization. Here is a good place to perform the normalization since it works even after the user filter the clusters of interest.
            
            sc_data2 <- NormalizeData(sc_data2,
                                      assay = "RNA",
                                      normalization.method = "LogNormalize",
                                      scale.factor = 10000)
            
            all_genes <- rownames(sc_data2@assays$RNA)
            sc_data2 <- Seurat::ScaleData(sc_data2,
                                          assay = "RNA",
                                          features = all_genes)
            
        }
        
        sc_data2
        
    })
    
    observeEvent(c(input$run_pca, input$rerun_after_filtering), {
        
        output$n_of_PCAs <- renderPlot({
            
            data_sc <- req( single_cell_data_pca() )
            
            Seurat::ElbowPlot(data_sc, ndims = 50, reduction = "pca")
            
        })
        
    })
    
    output$p2_down <- downloadHandler(
        
        filename = function() {
            paste("Elbow_Plot", ".", input$p2_format, sep = "")
        },
        content = function(file) {
            
            withProgress(message = "Please wait, preparing the data for download.",
                         value = 0.5, {
                             
                             shinyFeedback::feedbackWarning("p2_height", is.na(input$p2_height), "Required value")
                             shinyFeedback::feedbackWarning("p2_width", is.na(input$p2_width), "Required value")
                             
                             height <- as.numeric( req( input$p2_height) )
                             width <- as.numeric( req( input$p2_width) )
                             res <- as.numeric(input$p2_res)
                             
                             p <- Seurat::ElbowPlot(req( single_cell_data_pca()), ndims = 50, reduction = "pca")
                             
                             ggplot2::ggsave(file,
                                             p,
                                             height=height,
                                             width=width,
                                             units="cm",
                                             dpi=res)
                             
                         })
        }
        
    )
    
    ############################################################
    ## This section allows the filtering of clusters of cells ##
    ############################################################
    clusters_single_cell_data_reso_umap <- reactive({
        
        sc_data <- req( single_cell_data_reso_umap() )
        as.numeric( unique( as.character(sc_data@meta.data$seurat_clusters ) ) )
        
    })
    
    output$cluster_list_ui <- renderUI ({
        
        clusters <- req( clusters_single_cell_data_reso_umap() )
        
        shinyWidgets::pickerInput(
            inputId = "cluster_list",
            label = "Choose clusters to select or exclude",
            choices = sort(clusters),
            multiple = TRUE,
            options = list(`actions-box` = TRUE)
        )
        
    })
    
    to_filter <- reactive({
        
        shinyFeedback::feedbackWarning("cluster_list",
                                       is.null(input$cluster_list),
                                       "Required value")
        req(input$cluster_list)
        
        if ( input$filter_clusters_opt == 0 ) { # "select"
            
            to_filter <- base::subset(req( single_cell_data_reso_umap() ),
                                      idents = as.numeric(input$cluster_list))
            
            to_filter_ch <- to_filter@meta.data
            to_filter_ch <- base::rownames(to_filter_ch)
            
            to_filter_ch
            
        } else if ( input$filter_clusters_opt == 1 ) { # exclude
            
            to_filter <- base::subset(req( single_cell_data_reso_umap() ),
                                      idents = as.numeric(input$cluster_list),
                                      invert = TRUE)
            
            to_filter_ch <- to_filter@meta.data
            to_filter_ch <- base::rownames(to_filter_ch)
            
            to_filter_ch
        }
        
    })
    
    ############################################################
    
    single_cell_data_reso_umap <- eventReactive(input$run_clustering, {
        
        shinyFeedback::feedbackWarning("n_of_PCs", is.na(input$n_of_PCs), "Required value")
        shinyFeedback::feedbackWarning("resolution_clust", is.na(input$resolution_clust), "Required value")
        
        req(input$n_of_PCs)
        req(input$resolution_clust)
        
        showNotification("Running the clustering step",
                         duration = NULL,
                         id = "m6")
        
        data_sc <- req( single_cell_data_pca() )
        
        data_sc <- Seurat::FindNeighbors(data_sc, dims = 1:input$n_of_PCs)
        
        data_sc <- Seurat::FindClusters(data_sc, resolution = input$resolution_clust)
        
        sc_data <- Seurat::RunUMAP(data_sc, dims = 1:input$n_of_PCs)
        sc_data <- Seurat::RunTSNE(sc_data, dims = 1:input$n_of_PCs)
        
        sc_data
    })
    
    observeEvent(input$run_clustering, {
        
        output$tSNE <- renderPlot({
            
            Seurat::DimPlot(req( single_cell_data_reso_umap() ), reduction = "tsne", label = T, pt.size = .1)
            
        })
        
        output$umap <- renderPlot({
            
            Seurat::DimPlot(req( single_cell_data_reso_umap()), reduction = "umap", label = T, pt.size = .1)
            
        })
        
        output$cluster_size <-  renderPrint({
            
            sing_cell_data <- req( single_cell_data_reso_umap() )
            
            sc_meta <- as.data.frame(sing_cell_data[[]])
            sc_meta$cellcluster <- base::rownames(sc_meta)
            sc_meta <- sc_meta[, c("cellcluster", "seurat_clusters")]
            
            on.exit(removeNotification(id = "m6"), add = TRUE)
            
            base::table(sc_meta$seurat_clusters)
        })
        
    })
    
    output$p3_down <- downloadHandler(
        
        filename = function() {
            paste("clustering_plot_",  input$p3_down_opt, ".", input$p3_format, sep = "")
        },
        content = function(file) {
            
            withProgress(message = "Please wait, preparing the data for download.",
                         value = 0.5, {
                             
                             shinyFeedback::feedbackWarning("p3_height", is.na(input$p3_height), "Required value")
                             shinyFeedback::feedbackWarning("p3_width", is.na(input$p3_width), "Required value")
                             
                             height <- as.numeric( req( input$p3_height) )
                             width <- as.numeric( req( input$p3_width) )
                             res <- as.numeric(input$p3_res)
                             
                             if ( input$p3_down_opt == "UMAP" ) {
                                 
                                 p <- Seurat::DimPlot( req( single_cell_data_reso_umap() ), reduction = "umap", label = T, pt.size = .1)
                                 
                             } else if ( input$p3_down_opt == "t-SNE") {
                                 
                                 p <- Seurat::DimPlot(req( single_cell_data_reso_umap() ), reduction = "tsne", label = T, pt.size = .1)
                                 
                             }
                             
                             ggplot2::ggsave(file,
                                             p,
                                             height=height,
                                             width=width,
                                             units="cm",
                                             dpi=res)
                             
                         })
        }
        
    )
    
    output$downloadRDS <- downloadHandler(
        
        filename = function() {
            paste("clustered_dataset", ".rds", sep = "")
        },
        content = function(file) {
            
            withProgress(message = "Please wait, preparing the data for download.",
                         value = 0.5, {
                             
                             saveRDS( req( single_cell_data_reso_umap() ), file)
                             
                         })
        }
        
    )
    
    # We only change the assay here because it will allow the users to filter the clusters and still use the SCT assay for the PCA and clustering (when using SCTransform normalization. It has no affect if using lognormalization).
    single_cell_data_reso_umap_to_DE_vis <- reactive({
        
        sc_data <- single_cell_data_reso_umap()
        
        DefaultAssay(sc_data) <- "RNA"
        sc_data
        
    })
    
    output$find_markers_clust_id_tab1_ui <- renderUI ({
        
        clusters <- req( clusters_single_cell_data_reso_umap() )
        
        shinyWidgets::pickerInput(
            inputId = "find_markers_clust_id_tab1",
            label = "Select the cluster of interest",
            choices = sort(clusters),
            multiple = FALSE,
            options = list(`actions-box` = TRUE))
        
    })
    
    output$find_markers_clust_ID1_tab1_ui <- renderUI ({
        
        clusters <- req( clusters_single_cell_data_reso_umap() )
        
        shinyWidgets::pickerInput(
            inputId = "find_markers_clust_ID1_tab1",
            label = "Select the cluster of interest",
            choices = sort(clusters),
            multiple = F,
            options = list(`actions-box` = TRUE))
        
    })
    
    output$find_markers_clust_ID2_tab1_ui <- renderUI ({
        
        clusters <- req( clusters_single_cell_data_reso_umap() )
        
        # Exclude the cluster already selected in the option 1, since the two cluster must be different
        clusters <- clusters[ !clusters == input$find_markers_clust_ID1_tab1]
        
        shinyWidgets::pickerInput(
            inputId = "find_markers_clust_ID2_tab1",
            label = "Select the cluster(s) to compare",
            choices = sort(clusters),
            multiple = T,
            options = list(`actions-box` = TRUE,
                           "max-options-group" = 1),
            selected = sort(clusters)[1])
        
    })
    
    markers_tab1 <- eventReactive( input$run_ident_markers_tab1, {
        
        shinyFeedback::feedbackWarning("find_markers_tab1_return_thresh",
                                       is.na(input$find_markers_tab1_return_thresh),
                                       "Required value")
        shinyFeedback::feedbackWarning("find_markers_tab1_logfc_threshold",
                                       is.na(input$find_markers_tab1_logfc_threshold),
                                       "Required value")
        shinyFeedback::feedbackWarning("find_markers_tab1_min_pct",
                                       is.na(input$find_markers_tab1_min_pct),
                                       "Required value")
        
        req(input$find_markers_tab1_return_thresh)
        req(input$find_markers_tab1_logfc_threshold)
        req(input$find_markers_tab1_min_pct)
        
        sc_data <- req( single_cell_data_reso_umap() )
        
        if (input$find_markers_tab1_opt == 2) { #distinguishing a cluster from other
            
            ident_1 = input$find_markers_clust_ID1_tab1
            ident_2 = input$find_markers_clust_ID2_tab1
            
        } else if (input$find_markers_tab1_opt == 1 ) {
            
            ident_1 = input$find_markers_clust_id_tab1
            
        }
        
        showNotification("Identifing markers or D.E. genes",
                         duration = NULL,
                         id = "tab1_n1")
        
        finding_markers("finding_markers_tab1",
                        sc_data,
                        assay_choice = "RNA",
                        find_markers_tab1_opt = input$find_markers_tab1_opt,
                        find_markers_tab1_return.thresh = input$find_markers_tab1_return_thresh,
                        find_markers_tab1_logfc.threshold = input$find_markers_tab1_logfc_threshold,
                        find_markers_tab1_min.pct = input$find_markers_tab1_min_pct,
                        find_markers_tab1_test.use = input$find_markers_tab1_test.use,
                        ident.1 = ident_1,
                        ident.2 = ident_2,
                        find_markers_tab1_filt_pos = input$find_markers_tab1_filt_pos
        )
        
    })
    
    output$markers_tab1_react <- renderReactable({
        
        markers_tab1 <- req( markers_tab1() )
        
        if (input$find_markers_tab1_opt == 1) {
            
            markers_tab1$cluster <- isolate(input$find_markers_clust_id_tab1)
            markers_tab1 <- markers_tab1[, c( 1, ncol(markers_tab1), 2 : ( ncol(markers_tab1) - 1) ) ]
            
        } else if (input$find_markers_tab1_opt == 2) {
            
            markers_tab1$cluster <- paste0( isolate(input$find_markers_clust_ID1_tab1),
                                            "_vs_" ,
                                            isolate(input$find_markers_clust_ID2_tab1))
            markers_tab1 <- markers_tab1[, c( 1, ncol(markers_tab1), 2 : ( ncol(markers_tab1) - 1) ) ]
            
        } else if (input$find_markers_tab1_opt == 0) {
            
            markers_tab1 <- markers_tab1[, c(1, ncol(markers_tab1), 2 : ( ncol(markers_tab1) - 1) ) ]
            
        }
        
        on.exit(removeNotification(id = "tab1_n1"), add = TRUE)
        
        my_reactable(markers_tab1)
        
    })
    
    output$download_markers_tab1 <- downloadHandler(
        
        filename = function() {
            
            if ( input$find_markers_tab1_opt == 0 ) {
                
                paste("list_of_markers_all_clusters", ".csv", sep = "")
                
            } else if (input$find_markers_tab1_opt == 1) {
                
                paste("list_of_markers_", "cluster", input$find_markers_clust_id_tab1, ".csv", sep = "")
                
            } else if (input$find_markers_tab1_opt == 2) {
                
                paste("list_of_markers_", "cluster", input$find_markers_clust_ID1_tab1, "vs" , input$find_markers_clust_ID2_tab1, ".csv", sep = "")
                
            }
            
        },
        content = function(file) {
            
            withProgress(message = "Please wait, preparing the data for download.",
                         value = 0.5, {
                             
                             markers_tab1 <- req( markers_tab1() )
                             
                             if (input$find_markers_tab1_opt == 1) {
                                 
                                 markers_tab1$cluster <- input$find_markers_clust_id_tab1
                                 markers_tab1 <- markers_tab1[, c( 1, ncol(markers_tab1), 2 : ( ncol(markers_tab1) - 1) ) ]
                                 
                             } else if (input$find_markers_tab1_opt == 2) {
                                 
                                 markers_tab1$cluster <- paste0(input$find_markers_clust_ID1_tab1, "_vs_" , input$find_markers_clust_ID2_tab1)
                                 markers_tab1 <- markers_tab1[, c( 1, ncol(markers_tab1), 2 : ( ncol(markers_tab1) - 1) ) ]
                                 
                             } else if (input$find_markers_tab1_opt == 0) {
                                 
                                 markers_tab1 <- markers_tab1[, c(1, ncol(markers_tab1), 2 : ( ncol(markers_tab1) - 1) ) ]
                                 
                             }
                             
                             write.csv(markers_tab1, file, row.names = FALSE)
                             
                         })
        }
        
    )
    
    features <- reactive({
        
        req(input$markers_list)
        ext <- tools::file_ext(input$markers_list$name)
        markers_list_file <- input$markers_list$datapath
        
        if ( input$markers_list_header_opt == "" ) {
            
            shinyFeedback::feedbackWarning("markers_list_header_opt",
                                           TRUE,
                                           "Please, inform if the file has a header.")
            
        }
        
        req(input$markers_list_header_opt)
        
        read_file_markers("tab1_readfile",
                          markers_list_file = markers_list_file,
                          feed_ID ="markers_list",
                          ext = ext,
                          header_opt = input$markers_list_header_opt)
    })
    
    observeEvent(input$load_markers, {
        
        output$marker_group_selec = renderUI({
            
            pickerInput_markers_group("features_group", genes = req( features() ) )
            
        })
        
        output$marker_genes_selec = renderUI({
            
            features_group <- req(input$features_group)
            id_choice <- req(input$genes_ids)
            
            pickerInput_markers_genes("selected_genes",
                                      genes = req( features() ),
                                      features_group = features_group,
                                      id_choice = id_choice)
        })
        
    })
    
    filt_features <- eventReactive(input$run_heatmap, {
        
        features_f <- req( features() )
        features_group <- req(input$features_group)
        
        features_f <- dplyr::filter(features_f,
                                    Group %in% features_group)
        
        if (input$filter_genes_q == 0) {
            
            shinyFeedback::feedbackWarning("selected_genes",
                                           is.null(input$selected_genes),
                                           "Please, select one or more genes")
            req(input$selected_genes)
            
            if (input$genes_ids == "ID") {
                
                features_f <- features_f[features_f$GeneID %in% input$selected_genes, ]
                
            } else if (input$genes_ids == "name") {
                
                features_f <- features_f[features_f$Name %in% input$selected_genes, ]
                
            }
            
        }
        
        features_f
        
    })
    
    sc_data_av_react <- eventReactive(input$run_heatmap, {
        
        sc_data_av <- Seurat::AverageExpression( req( single_cell_data_reso_umap() ),
                                                 assays = "RNA",
                                                 slot = input$slot_selection_heatmap)
        
        sc_data_av <- as.matrix(sc_data_av[[1]])
        
        on.exit(removeNotification(id = "m7"), add = TRUE)
        sc_data_av
        
    })
    
    observeEvent(input$run_heatmap, {
        
        # This calculate the height of the plot based on the n of genes so the plot can be adjusted to fit all genes
        heatmap_n_genes <- reactive({
            
            req(filt_features())
            
            heatmap_genes(features = filt_features(), sc_data_av = req( sc_data_av_react()) )
            
        })
        heatmap_Height <- reactive( 150 + ( 20 * req( heatmap_n_genes() ) ) )
        
        heat_map_prep <- reactive({
            req(filt_features())
            
            features <- filt_features()
            sc_data_av <- req( sc_data_av_react() )
            
            features_selec <- as.data.frame(unique(features$GeneID))
            
            if (nrow(features) == 1) {
                
                sc_data_av_feat <- sc_data_av[base::rownames(sc_data_av) %in% features_selec[, 1], ]
                sc_data_av_feat <- as.matrix(t(sc_data_av_feat))
                base::rownames(sc_data_av_feat) <- features_selec[1, 1]
                
            } else {
                
                sc_data_av_feat <- sc_data_av[base::rownames(sc_data_av) %in% features_selec[, 1], ]
                
            }
            sc_data_av_feat
        })
        
        output$heat_map <- renderPlot({
            
            req(heat_map_prep())
            heat_map_prep <- heat_map_prep()
            
            ComplexHeatmap::Heatmap(heat_map_prep, border = TRUE,
                                    rect_gp = gpar(col = "white", lwd = 2),
                                    column_title = "Clusters",
                                    column_title_side = "bottom",
                                    name = "Expression",
                                    show_row_dend = T)
            
            
        }, height = heatmap_Height())
        
        output$heat_map_ui <- renderUI({
            
            plotOutput("heat_map", height = req( heatmap_Height() ) )
            
        })
        
        output$marker_to_feature_plot = renderUI({
            
            feature_plots_gene_selection("selected_genes_for_feature_plot",
                                         genes_names = req( filt_features() ),
                                         sc_data_av_react = req( sc_data_av_react() ),
                                         genes_ids = input$genes_ids)
            
        })
        
        output$p4_down <- downloadHandler(
            
            filename = function() {
                paste("Heatmap", ".", input$p4_format, sep = "")
            },
            content = function(file) {
                
                withProgress(message = "Please wait, preparing the data for download.",
                             value = 0.5, {
                                 
                                 shinyFeedback::feedbackWarning("p4_height", is.na(input$p4_height), "Required value")
                                 shinyFeedback::feedbackWarning("p4_width", is.na(input$p4_width), "Required value")
                                 
                                 height <- as.numeric( req( input$p4_height) )
                                 width <- as.numeric( req( input$p4_width) )
                                 res <- as.numeric(input$p4_res)
                                 
                                 p <- ComplexHeatmap::Heatmap(req( heat_map_prep() ), border = TRUE,
                                                              rect_gp = gpar(col = "white", lwd = 2),
                                                              column_title = "Clusters",
                                                              column_title_side = "bottom",
                                                              name = "Expression",
                                                              show_row_dend = T)
                                 
                                 # Complex heatmap does not work well with ggplot2::ggsave. So, we use grid.grap to make it compatible.
                                 gb = grid.grabExpr(draw(p))
                                 
                                 ggplot2::ggsave(file,
                                                 gb,
                                                 height=height,
                                                 width=width,
                                                 units="cm",
                                                 dpi=res)
                                 
                             })
            }
            
        )
        
    })
    
    features_selec <- eventReactive(input$run_feature_plot, {
        
        shinyFeedback::feedbackWarning("selected_genes_for_feature_plot",
                                       is.null(input$selected_genes_for_feature_plot),
                                       "Please, select one or more genes.")
        req(input$selected_genes_for_feature_plot)
        
        features <- req( filt_features() )
        
        if (input$genes_ids == "ID") {
            
            features_f <- dplyr::filter(features,
                                        GeneID %in% input$selected_genes_for_feature_plot)
            
        }
        else if (input$genes_ids == "name") {
            
            features_f <- dplyr::filter(features,
                                        Name %in% input$selected_genes_for_feature_plot)
            
        }
        
        as.character(unique(features_f$GeneID))
        
    })
    
    observeEvent(input$run_feature_plot, {
        
        output$feature_plot <- renderUI({
            
            feat_length <- length(req( features_selec() ) )
            
            plot_output_list <- lapply(1:feat_length, function(i) {
                plotname <- paste("plot1", i, sep="")
                plotOutput(plotname, height = 300)
            })
            
            # Convert the list to a tagList - this is necessary for the list of items
            # to display properly.
            do.call(tagList, plot_output_list)
        })
        
        output$feature_plot_dark <- renderUI({
            
            feat_length <- length(req( features_selec()))
            
            plot_output_list <- lapply(1:feat_length, function(i) {
                plotname <- paste("plot2", i, sep="")
                plotOutput(plotname, height = 300)
            })
            
            do.call(tagList, plot_output_list)
        })
        
        output$umap2 <- renderUI({
            
            feat_length <- length(req( features_selec()))
            
            plot_output_list <- lapply(1:feat_length, function(i) {
                plotname <- paste("plot3", i, sep="")
                plotOutput(plotname, height = 300)
            })
            
            do.call(tagList, plot_output_list)
        })
        
        output$run_vln_plot <- renderUI({
            
            feat_length <- length(req( features_selec()))
            
            plot_output_list <- lapply(1:feat_length, function(i) {
                plotname <- paste("plot4", i, sep="")
                plotOutput(plotname, height = 300)
            })
            
            do.call(tagList, plot_output_list)
        })
        
        output$run_dot_plot <- renderUI({
            
            feat_length <- length(req( features_selec()))
            
            plot_output_list <- lapply(1:feat_length, function(i) {
                plotname <- paste("plot5", i, sep="")
                plotOutput(plotname, height = 300)
            })
            
            do.call(tagList, plot_output_list)
        })
        
        sc_data <- req( single_cell_data_reso_umap() )
        for (i in 1:length( req( features_selec() ) ) ) {
            
            features <- req( features_selec() )
            local({
                
                my_i <- i
                
                plotname <- paste("plot1", my_i, sep="")
                output[[plotname]] <- renderPlot({
                    
                    assay_id <- "RNA"
                    
                    minimal <- min(sc_data[[assay_id]]@data[features[my_i], ])
                    maximal <- max(sc_data[[assay_id]]@data[features[my_i], ])
                    
                    suppressMessages( Seurat::FeaturePlot(sc_data,
                                                          cols = c("lightgrey", "red"),
                                                          features = features[my_i],
                                                          slot = input$slot_selection_feature_plot,
                                                          reduction = "umap") +
                                          scale_colour_gradient2(limits=c(minimal, maximal),
                                                                 midpoint = maximal / 2,
                                                                 low = "gray80",
                                                                 mid = "gold",
                                                                 high = "red") )
                    
                })
                
                plotname <- paste("plot2", my_i, sep="")
                output[[plotname]] <- renderPlot({
                    
                    assay_id <- "RNA"
                    
                    minimal <- min(sc_data[[assay_id]]@data[features[my_i], ])
                    maximal <- max(sc_data[[assay_id]]@data[features[my_i], ])
                    
                    suppressMessages( Seurat::FeaturePlot(sc_data,
                                                          cols = c("lightgrey", "red"),
                                                          features = features[my_i],
                                                          slot = input$slot_selection_feature_plot,
                                                          reduction = "umap") +
                                          scale_colour_gradient2(limits=c(minimal, maximal),
                                                                 midpoint = maximal / 2,
                                                                 low = "gray80",
                                                                 mid = "gold",
                                                                 high = "red") ) +
                        DarkTheme()
                    
                })
                
                plotname <- paste("plot3", my_i, sep="")
                output[[plotname]] <- renderPlot({
                    
                    Seurat::DimPlot(req( single_cell_data_reso_umap() ), reduction = "umap", label = T, pt.size = .1)
                    
                })
                
                plotname <- paste("plot4", my_i, sep="")
                output[[plotname]] <- renderPlot({
                    
                    Seurat::VlnPlot(sc_data,
                                    features = features[my_i],
                                    slot = input$slot_selection_feature_plot)
                    
                })
                
                plotname <- paste("plot5", my_i, sep="")
                output[[plotname]] <- renderPlot({
                    
                    Seurat::DotPlot(sc_data,
                                    features = features[my_i],
                                    cols = c("lightgrey", "red")) +
                        scale_x_discrete(position = "top") +
                        xlab("")+
                        ylab("") +
                        ggtitle("") +
                        theme(
                            axis.title.x=element_blank(),
                            axis.title.y=element_blank(),
                            axis.text.y = element_text( size = rel(1) ),
                            axis.text.x = element_text( angle = 90) ) +
                        coord_flip()
                    
                })
                
            })
            
        }
        
        showNotification("Generating additional plots",
                         duration = 25,
                         id = "m8")
        
        output$select_genes_add_plot_to_down_ui <- renderUI({
            
            filt_features <- req(input$selected_genes_for_feature_plot)
            
            div(class = "option-group",
                shinyWidgets::pickerInput(
                    inputId = "select_genes_add_plot_to_down",
                    label = "Select the genes that you want to download",
                    choices = sort(filt_features),
                    multiple = TRUE,
                    options = list(`actions-box` = TRUE)
                ))
            
        })
        
    })
    
    observeEvent(input$start_down_add_plots_tab1, {
        
        on.exit(removeNotification(id = "m8"), add = TRUE)
        
        shinyFeedback::feedbackWarning("select_genes_add_plot_to_down",
                                       is.null(input$select_genes_add_plot_to_down),
                                       "Please, select one or more genes.")
        
        genes <- req(input$select_genes_add_plot_to_down)
        
        withProgress(message = "Please wait, preparing the data for download.",
                     value = 0.5, {
                         
                         if (file.exists("./images") == F ) {
                             dir.create('./images')
                         }
                         
                         # Creates new folder to keep the plots organized. Adding date and time helps to avoid overwriting the data
                         path_new <- paste0("./images/one_sample_plots", format(Sys.time(),'_%Y-%m-%d__%H%M%S'))
                         
                         dir.create(path_new)
                         dir.create( paste0(path_new,"/feature_plots") )
                         dir.create( paste0(path_new,"/violin_plots") )
                         dir.create( paste0(path_new, "/dot_plots") )
                         
                         sc_data <- req( single_cell_data_reso_umap() )
                         assay_id <- "RNA"
                         
                         for( i in 1:length(genes) ){
                             
                             minimal <- min(sc_data[[assay_id]]@data[genes[i], ])
                             maximal <- max(sc_data[[assay_id]]@data[genes[i], ])
                             
                             p <- suppressMessages(Seurat::FeaturePlot(sc_data,
                                                                       cols = c("lightgrey", "red"),
                                                                       features = genes[i],
                                                                       slot = input$slot_selection_feature_plot,
                                                                       reduction = "umap") +
                                                       scale_colour_gradient2(limits=c(minimal, maximal),
                                                                              midpoint = maximal / 2,
                                                                              low = "gray80",
                                                                              mid = "gold",
                                                                              high = "red"))
                             
                             file <- paste0(path_new, "/feature_plots/", genes[i], ".", input$add_p_tab1_feat_format)
                             
                             ggplot2::ggsave(file,
                                             p,
                                             height=input$add_p_tab1_feat_height,
                                             width=input$add_p_tab1_feat_width,
                                             units="cm",
                                             dpi=as.numeric(input$add_p_tab1_feat_res))
                             
                             file <- paste0(path_new, "/violin_plots/", genes[i], ".", input$add_p_tab1_violin_format)
                             
                             p <- Seurat::VlnPlot(sc_data,
                                                  features = genes[i],
                                                  slot = input$slot_selection_feature_plot)
                             
                             ggplot2::ggsave(file,
                                             p,
                                             height=input$add_p_tab1_violin_height,
                                             width=input$add_p_tab1_violin_width,
                                             units="cm",
                                             dpi=as.numeric(input$add_p_tab1_violin_res))
                             
                             file <- paste0(path_new, "/dot_plots/", genes[i], ".", input$add_p_tab1_dot_format)
                             
                             p <- Seurat::DotPlot(sc_data,
                                                  features = genes[i],
                                                  cols = c("lightgrey", "red")) +
                                 scale_x_discrete(position = "top") +
                                 xlab("")+
                                 ylab("") +
                                 ggtitle("") +
                                 theme(axis.title.x=element_blank(),
                                       axis.title.y=element_blank(),
                                       axis.text.y = element_text( size = rel(1) ),
                                       axis.text.x = element_text( angle = 90) ) +
                                 coord_flip()
                             
                             ggplot2::ggsave(file,
                                             p,
                                             height=input$add_p_tab1_dot_height,
                                             width=input$add_p_tab1_dot_width,
                                             units="cm",
                                             dpi=as.numeric(input$add_p_tab1_dot_res))
                         }
                         
                     })
        
        showNotification("All plots were downloaded!",
                         duration = 15,
                         id = "")
        
    })
    
    #####################################
    ### Tab 2 - Integration pipeline ####
    #####################################
    
    output$load_integrated_ui <- renderUI({
        
        rds_list <- list.files('./RDS_files/', pattern = "*.rds")
        
        div(class = "option-group",
            shinyWidgets::pickerInput(
                inputId = "load_integrated",
                label = "Select the file containing the integrated data",
                choices = sort(rds_list),
                multiple = FALSE,
                options = list(`actions-box` = TRUE)
            ))
        
    })
    
    samples_list_integration <- reactive({
        
        req(input$samples_list_integration)
        
        ext <- tools::file_ext(input$samples_list_integration$name)
        config_input_file <- input$samples_list_integration$datapath
        
        test <- ext %in% c(
            'text/csv',
            'text/comma-separated-values',
            'text/tab-separated-values',
            'csv',
            'tsv')
        
        shinyFeedback::feedbackDanger("samples_list_integration",
                                      test == F,
                                      "Format is not supported! Upload a CSV or TSV file.")
        
        config_input <- switch(ext,
                               csv = vroom::vroom(config_input_file, delim = ","),
                               tsv = vroom::vroom(config_input_file, delim = "\t"),
                               validate("Invalid file; Please upload a .csv or .tsv file")
        )
        config_input <- as.data.frame(config_input)
        
        if ( ncol(config_input) >= 6 ) {
            shinyFeedback::feedbackSuccess("samples_list_integration",
                                           T, "")
        } else {
            shinyFeedback::feedbackDanger("samples_list_integration",
                                          T,
                                          "The config. file must have at least six columns!")
            
            validate( need( ncol(config_input) >= 6,
                            "The config. file must have at least six columns!",
                            "samples_list_integration") )
        }
        
        config_input
        
    })
    
    output$select_sample_tab2 <- renderUI({
        
        dir_list <- req( samples_list_integration() )
        dir_list <- unique(dir_list[, 1])
        
        div(class = "option-group",
            shinyWidgets::pickerInput(
                inputId = "sample_folder_tab2",
                label = "Select the samples to use",
                choices = sort(as.character(dir_list)),
                multiple = TRUE,
                options = list(`actions-box` = TRUE),
            ))
        
    })
    
    single_cell_data_reac_tab2 <- eventReactive(input$load_rds_file, {
        
        if ( input$integration_options == 1) { # Load file
            
            ##  the user needs to tell what normalization was used, since we should not scale the data when using SCTransform.
            shinyFeedback::feedbackWarning("load_rds_int_normalization",
                                           input$load_rds_int_normalization == "",
                                           #T,
                                           "Please select an option")
            validate(need(input$load_rds_int_normalization != "",
                          message = "",
                          label = "load_rds_int_normalization"))
            
            showNotification("Loading the integrated data",
                             id = "m9",
                             duration = NULL)
            
            sc_data <- readRDS( paste0("./RDS_files/", req(input$load_integrated)) )
            
            on.exit(removeNotification(id = "m9"), add = TRUE)
            
            sc_data
            
        } else if (input$integration_options == 0 ) { # new analysis
            
            shinyFeedback::feedbackWarning("sample_folder_tab2", is.null(input$sample_folder_tab2), "Required value")
            shinyFeedback::feedbackWarning("int_regex_mito", !shiny::isTruthy(input$int_regex_mito), "Required value")
            
            selected_samples <- req(input$sample_folder_tab2)
            config_csv <- req( samples_list_integration() )
            int_regex_mito <- req(input$int_regex_mito)
            
            # Select only the samples that the user specified.
            config_csv <- config_csv[config_csv[, 1] %in% selected_samples, ]
            
            project_name <- input$int_project_name
            
            if ( input$normaliz_method_tab2 == 0 ) { # LogNormalize
                
                shinyFeedback::feedbackWarning("sample_folder_tab2", is.null(input$sample_folder_tab2), "Required value")
                shinyFeedback::feedbackWarning("scale_factor_tab2", is.na(input$scale_factor_tab2), "Required value")
                shinyFeedback::feedbackWarning("n_of_var_genes_integration", is.na(input$n_of_var_genes_integration), "Required value")
                shinyFeedback::feedbackWarning("n_of_PCs_integration", is.na(input$n_of_PCs_integration), "Required value")
                
                scale_factor_tab2 <- req(input$scale_factor_tab2)
                n_of_var_genes_integration <- req(input$n_of_var_genes_integration)
                n_of_PCs_integration <- req(input$n_of_PCs_integration)
                
                showNotification("Integrating the data. Please wait, it can take a few minutes.",
                                 id = "m12",
                                 duration = NULL)
                
                sc_data <- integration_lognorm(config_csv = config_csv,
                                               project_name = project_name,
                                               int_regex_mito = int_regex_mito,
                                               scale_factor_tab2 = scale_factor_tab2,
                                               most_var_method_integration = input$most_var_method_tab2,
                                               n_of_var_genes_integration = n_of_var_genes_integration,
                                               n_of_PCs_integration = n_of_PCs_integration)
                
                on.exit(removeNotification(id = "m12"), add = TRUE)
                
                sc_data
                
            } else if ( input$normaliz_method_tab2 == 1 ) { # SCTransform
                
                showNotification("Integrating the data. Please wait, it can take a few minutes.",
                                 id = "m12",
                                 duration = NULL)
                
                sc_data <- integration_sctransform(config_csv = config_csv,
                                                   project_name = project_name,
                                                   int_regex_mito = int_regex_mito)
                
                on.exit(removeNotification(id = "m12"), add = TRUE)
                
                sc_data
                
            }
        }
        
        sc_data
        
    })
    
    output$download_int_data <- downloadHandler(
        
        filename = function() {
            paste("Integrated_datasets_without_clutering", ".rds", sep = "")
        },
        content = function(file) {
            
            withProgress(message = "Please wait, preparing the data for download.",
                         value = 0.5, {
                             
                             saveRDS( req( single_cell_data_reac_tab2() ), file)
                             
                         })
            
        }
        
    )
    
    output$VlnPlot_tab2 <- renderPlot({
        
        data_set <- req( single_cell_data_reac_tab2() )
        
        Seurat::VlnPlot(data_set,
                        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                        ncol = 3,
                        split.plot = F)
        
    })
    
    single_cell_data_filt_tab2 <- reactive({
        
        data_sc <- req( single_cell_data_reac_tab2() )
        
        test_cond <- if( !is.na(input$max_count_tab2) && !is.na(input$min_count_tab2) ) {
            
            shinyFeedback::feedbackWarning("max_count_tab2",
                                           input$min_count_tab2 > input$max_count_tab2,
                                           "No cells will be selected by appling this parameters!")
            
            shinyFeedback::feedbackWarning("min_count_tab2",
                                           input$min_count_tab2 > input$max_count_tab2,
                                           "No cells will be selected by appling this parameters!")
            validate(need(input$min_count_tab2 < input$max_count_tab2, "Error: No cells will be selected by appling this parameters!"))
            
        }
        
        if ( !is.na(input$min_count_tab2) ) {
            
            data_sc <- base::subset(data_sc,
                                    subset = nFeature_RNA > input$min_count_tab2)
            
        }
        
        if ( !is.na(input$max_count_tab2) ) {
            
            shinyFeedback::feedbackWarning("max_count_tab2",
                                           input$max_count_tab2 <= 0,
                                           "No cells will be selected by appling this parameters!")
            
            validate(need(input$max_count_tab2 > 0, "Error: No cells will be selected by appling this parameters!"))
            
            data_sc <- base::subset(data_sc,
                                    subset =  nFeature_RNA < input$max_count_tab2)
            
        }
        
        if ( !is.na(input$max_mito_perc_tab2) ) {
            
            shinyFeedback::feedbackWarning("max_mito_perc_tab2", input$max_mito_perc_tab2 <= 0, "No cells will be selected by appling this parameters!")
            validate(need(input$max_mito_perc_tab2 > 0, "Error: No cells will be selected by appling this parameters!"))
            
            data_sc <- base::subset(data_sc,
                                    subset =  percent.mt < input$max_mito_perc_tab2)
            
        }
        
        data_sc
        
    })
    
    observeEvent(input$run_vinplot_tab2, {
        
        output$VlnPlot_filt_tab2 <- renderPlot({
            
            Seurat::VlnPlot(req( single_cell_data_filt_tab2() ),
                            features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                            ncol = 3,
                            split.plot = F)
            
        })
        
        output$p5_down <- downloadHandler(
            
            filename = function() {
                paste("VlnPlot", ".", input$p5_format, sep = "")
            },
            content = function(file) {
                
                withProgress(message = "Please wait, preparing the data for download.",
                             value = 0.5, {
                                 
                                 shinyFeedback::feedbackWarning("p5_height", is.na(input$p5_height), "Required value")
                                 shinyFeedback::feedbackWarning("p5_width", is.na(input$p5_width), "Required value")
                                 
                                 height <- as.numeric( req( input$p5_height) )
                                 width <- as.numeric( req( input$p5_width) )
                                 res <- as.numeric(input$p5_res)
                                 
                                 p <- Seurat::VlnPlot(req( single_cell_data_filt_tab2() ),
                                                      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                                                      ncol = 3,
                                                      split.plot = F)
                                 
                                 ggplot2::ggsave(file,
                                                 p,
                                                 height=height,
                                                 width=width,
                                                 units="cm",
                                                 dpi=res)
                                 
                             })
            }
            
        )
        
    })
    
    single_cell_data_scaled_tab2 <- eventReactive(
        c(input$run_pca_tab2_1,
          input$run_pca_tab2_2), {
              
              single_cell_data_filt_tab2 <- req( single_cell_data_filt_tab2() )
              
              # If loading the data and normalization is SCTransform, skip the scaling.
              if (input$integration_options == 1 && input$load_rds_int_normalization == 1) {
                  
                  ## Same if running new analysis and normalization is SCtransform
              } else if (input$integration_options == 0 && input$normaliz_method_tab2 == 1) {
                  
              } else { # lognormalization
                  
                  showNotification("Scalling the data",
                                   duration = NULL,
                                   id = "tab2_m4")
                  
                  single_cell_data_filt_tab2 <- Seurat::ScaleData(single_cell_data_filt_tab2,
                                                                  verbose = T)
                  
                  on.exit(removeNotification(id = "tab2_m4"), add = TRUE)
                  
              }

              single_cell_data_filt_tab2
          })
    
    single_cell_data_pca_tab2 <- eventReactive(c(input$run_pca_tab2_1,
                                                 input$run_pca_tab2_2,
                                                 input$rerun_after_filtering_tab2), {
                                                     
                                                     if (input$filter_clusters_tab2 == 1) {
                                                         
                                                         sc_data <- req( single_cell_data_scaled_tab2_filtered() )
                                                         
                                                     } else if ( input$filter_clusters_tab2 == 0 ) {
                                                         
                                                         sc_data <- req( single_cell_data_scaled_tab2() )
                                                         
                                                     }
                                                     
                                                     showNotification("Running PCA",
                                                                      duration = NULL,
                                                                      id = "m5")
                                                     
                                                     sc_data <-  RunPCA(sc_data,
                                                                        verbose = T)
                                                     
                                                     on.exit(removeNotification(id = "m5"), add = TRUE)
                                                     
                                                     sc_data
                                                     
                                                 })
    
    ########################################################
    ## Section that allows filtering clusters of interest ##
    ########################################################
    output$cluster_list_tab2_ui <- renderUI ({
        
        clusters <- req( clusters_single_cell_data_reso_umap_tab2() )
        
        shinyWidgets::pickerInput(
            inputId = "cluster_list_tab2",
            label = "Choose clusters to select or exclude",
            choices = sort(clusters),
            multiple = TRUE,
            options = list(`actions-box` = TRUE)
        )
        
    })
    
    to_filter_tab2 <- reactive({
        
        shinyFeedback::feedbackWarning("cluster_list_tab2", is.null(input$cluster_list_tab2), "Required value")
        req(input$cluster_list_tab2)
        
        sc_data <- req( single_cell_data_clustered() )
        
        if ( input$filter_clusters_opt_tab2 == 0 ) { #"select"
            
            to_filter <- base::subset(sc_data,
                                      idents = as.numeric(input$cluster_list_tab2)
            )
            
            to_filter_ch <- to_filter@meta.data
            to_filter_ch <- base::rownames(to_filter_ch)
            
            to_filter_ch
            
        } else if ( input$filter_clusters_opt_tab2 == 1 ) { #"exclude"
            
            to_filter <- base::subset(sc_data,
                                      idents = as.numeric(input$cluster_list_tab2),
                                      invert = TRUE)
            
            to_filter_ch <- to_filter@meta.data
            to_filter_ch <- base::rownames(to_filter_ch)
            
            to_filter_ch
        }
        
    })
    
    single_cell_data_scaled_tab2_filtered <- eventReactive(input$rerun_after_filtering_tab2, {
        
        sc_data <- req( single_cell_data_filt_tab2() )
        cells_to_filter <- req( to_filter_tab2() )
        
        if ( input$filter_clusters_tab2 == 1 ) {
            
            sc_data <- base::subset(sc_data,
                                    cells = cells_to_filter)
            
        }
        
        if (input$integration_options == 1 && input$load_rds_int_normalization == 1) {
            
        } else if (input$integration_options == 0 && input$normaliz_method_tab2 == 1) {
            
        } else { # lognormalization
            
            showNotification("Scalling the data",
                             duration = NULL,
                             id = "tab2_m4")
            
            sc_data <- Seurat::ScaleData(sc_data, verbose = T)
            
            on.exit(removeNotification(id = "tab2_m4"), add = TRUE)
            
            sc_data
        }
        
        sc_data
    })
    
    ########################################################
    
    observeEvent(c(input$run_pca_tab2_1,
                   input$run_pca_tab2_2,
                   input$rerun_after_filtering_tab2), {
                       
                       output$n_of_PCAs_tab2 <- renderPlot({
                           
                           data_sc <- req( single_cell_data_pca_tab2() )
                           
                           Seurat::ElbowPlot(data_sc, ndims = 50, reduction = "pca")
                           
                       })
                       
                       output$p6_down <- downloadHandler(
                           
                           filename = function() {
                               paste("ElbowPlot", ".", input$p6_format, sep = "")
                           },
                           content = function(file) {
                               
                               withProgress(message = "Please wait, preparing the data for download.",
                                            value = 0.5, {
                                                
                                                shinyFeedback::feedbackWarning("p6_height", is.na(input$p6_height), "Required value")
                                                shinyFeedback::feedbackWarning("p6_width", is.na(input$p6_width), "Required value")
                                                
                                                height <- as.numeric( req( input$p6_height) )
                                                width <- as.numeric( req( input$p6_width) )
                                                res <- as.numeric(input$p6_res)
                                                
                                                data_sc <- req( single_cell_data_pca_tab2() )
                                                
                                                p <- Seurat::ElbowPlot(data_sc, ndims = 50, reduction = "pca")
                                                
                                                ggplot2::ggsave(file,
                                                                p,
                                                                height=height,
                                                                width=width,
                                                                units="cm",
                                                                dpi=res)
                                                
                                            })
                           }
                           
                       )
                       
                   })
    
    single_cell_data_clustered <- eventReactive(input$run_clustering_tab2, {
        
        shinyFeedback::feedbackWarning("n_of_PCs_tab2", is.na(input$n_of_PCs_tab2), "Required value")
        shinyFeedback::feedbackWarning("resolution_clust_tab2", is.na(input$resolution_clust_tab2), "Required value")
        
        req(input$n_of_PCs_tab2)
        req(input$resolution_clust_tab2)
        
        showNotification("Running the clustering step",
                         duration = NULL,
                         id = "m6")
        
        sc_data <- req( single_cell_data_pca_tab2() )
        
        sc_data <- Seurat::RunUMAP(sc_data,
                                   reduction = "pca",
                                   dims = 1:input$n_of_PCs_tab2)
        
        sc_data <- Seurat::FindNeighbors(sc_data,
                                         reduction = "pca",
                                         dims = 1:input$n_of_PCs_tab2)
        
        sc_data <- Seurat::FindClusters(sc_data,
                                        resolution = input$resolution_clust_tab2)
        
        sc_data
        
    })
    
    observeEvent(input$run_clustering_tab2, {
        
        output$umap_tab2 <- renderPlot({
            
            Seurat::DimPlot(req( single_cell_data_clustered() ), reduction = "umap", label = T, pt.size = .1)

        })
        
        output$umap_three_samples_comb <- renderPlot({
            
            Seurat::DimPlot(req( single_cell_data_clustered() ),
                            reduction = "umap",
                            label = F,
                            group.by = "treat",
                            pt.size = .1
            )
            
        })
        
        output$umap_three_samples <- renderPlot({
            
            Seurat::DimPlot(req( single_cell_data_clustered() ),
                            reduction = "umap",
                            label = T,
                            split.by = "treat", pt.size = .1
            )
        })
        
        output$p7_down <- downloadHandler(
            
            filename = function() {
                paste("Clustering_plot_", input$p7_down_opt, ".", input$p7_format, sep = "")
            },
            content = function(file) {
                
                withProgress(message = "Please wait, preparing the data for download.",
                             value = 0.5, {
                                 
                                 shinyFeedback::feedbackWarning("p7_height", is.na(input$p7_height), "Required value")
                                 shinyFeedback::feedbackWarning("p7_width", is.na(input$p7_width), "Required value")
                                 
                                 height <- as.numeric( req( input$p7_height) )
                                 width <- as.numeric( req( input$p7_width) )
                                 res <- as.numeric(input$p7_res)
                                 
                                 if ( input$p7_down_opt == "UMAP" ) {
                                     
                                     p <- Seurat::DimPlot(req( single_cell_data_clustered() ),
                                                          reduction = "umap",
                                                          label = T,
                                                          pt.size = .1)
                                     
                                 } else if ( input$p7_down_opt == "UMAP1" ) {
                                     
                                     p <- Seurat::DimPlot(req( single_cell_data_clustered() ),
                                                          reduction = "umap",
                                                          label = F,
                                                          group.by = "treat",
                                                          pt.size = .1
                                     )
                                     
                                 } else if ( input$p7_down_opt == "UMAP2" ) {
                                     
                                     p <- Seurat::DimPlot(req( single_cell_data_clustered() ),
                                                          reduction = "umap",
                                                          label = T,
                                                          split.by = "treat", pt.size = .1
                                     )
                                     
                                 }
                                 
                                 ggplot2::ggsave(file,
                                                 p,
                                                 height=height,
                                                 width=width,
                                                 units="cm",
                                                 dpi=res)
                                 
                             })
            }
            
        )
        
        output$cluster_size_tab2 <- renderPrint({
            
            sing_cell_data <- req( single_cell_data_clustered() )
            
            sc_meta <- as.data.frame(sing_cell_data[[]])
            sc_meta$cellcluster <- base::rownames(sc_meta)
            sc_meta <- sc_meta[, c("cellcluster", "seurat_clusters")]
            
            # returns the number of cells in each cluster
            cell_per_cluster <- sc_meta$seurat_clusters
            
            on.exit(removeNotification(id = "m6"), add = TRUE)
            
            table(sc_meta$seurat_clusters)
            
            
        })
        
    } )
    
    clusters_single_cell_data_reso_umap_tab2 <- reactive({
        
        sc_data <- req( single_cell_data_clustered() )
        as.numeric( unique( as.character(sc_data@meta.data$seurat_clusters ) ) )
        
    })
    
    treat_single_cell_data_reso_umap_tab2 <- reactive({
        
        sc_data <- req( single_cell_data_clustered() )
        unique( as.character(sc_data@meta.data$treat ) )
        
    })
 
    output$downloadRDS_tab2 <- downloadHandler(
        
        filename = function() {
            paste("Integrated_datasets_after_clutering", ".rds", sep = "")
        },
        content = function(file) {
            
            withProgress(message = "Please wait, preparing the data for download.",
                         value = 0.5, {
                             
                             saveRDS( req( single_cell_data_clustered() ), file)
                             
                         })
            
        }
        
    )
    
    output$find_markers_clust_id_tab2_ui <- renderUI ({
        
        clusters <- req( clusters_single_cell_data_reso_umap_tab2() )
        
        shinyWidgets::pickerInput(
            inputId = "find_markers_clust_id_tab2",
            label = "Select the cluster of interest",
            choices = sort(clusters),
            multiple = FALSE,
            options = list(`actions-box` = TRUE)
        )
        
    })
    
    output$find_markers_clust_ID1_tab2_ui <- renderUI ({
        
        clusters <- req( clusters_single_cell_data_reso_umap_tab2() )
        
        shinyWidgets::pickerInput(
            inputId = "find_markers_clust_ID1_tab2",
            label = "Select the cluster of interest",
            choices = sort(clusters),
            multiple = FALSE,
            options = list(`actions-box` = TRUE)
        )
        
    })
    
    output$find_markers_clust_ID2_tab2_ui <- renderUI ({
        
        clusters <- req( clusters_single_cell_data_reso_umap_tab2() )
        
        # Exclude the cluster already selected in the option 1, since the two cluster must be different
        clusters <- clusters[ !clusters == req(input$find_markers_clust_ID1_tab2)]
        
        shinyWidgets::pickerInput(
            inputId = "find_markers_clust_ID2_tab2",
            label = "Select the cluster(s) to compare",
            choices = sort(clusters),
            multiple = T,
            options = list(`actions-box` = TRUE),
            selected = sort(clusters)[1]
        )
        
    })
    
    output$find_markers_or_DE_tab2_cluster_ui <- renderUI ({
        
        clusters <- req( clusters_single_cell_data_reso_umap_tab2() )
        
        shinyWidgets::pickerInput(
            inputId = "find_markers_or_DE_tab2_cluster",
            label = "Select the cluster of interest",
            choices = sort(clusters),
            multiple = FALSE,
            options = list(`actions-box` = TRUE)
        )
        
    })
    
    output$find_markers_or_DE_tab2_treat1_ui <- renderUI({
        
        treat <- req( treat_single_cell_data_reso_umap_tab2() )
        
        shinyWidgets::pickerInput(
            inputId = "find_markers_or_DE_tab2_treat1",
            label = "Select the sample/treatment of interest",
            choices = treat,
            multiple = FALSE,
            options = list(`actions-box` = TRUE)
        )
        
    })
    
    output$find_markers_or_DE_tab2_treat2_ui  <- renderUI({
        
        treat <- req( treat_single_cell_data_reso_umap_tab2() )
        treat <- treat[!treat %in% req(input$find_markers_or_DE_tab2_treat1)]
        
        shinyWidgets::pickerInput(
            inputId = "find_markers_or_DE_tab2_treat2",
            label = "Select the sample/treatment of interest",
            choices = treat,
            multiple = T,
            options = list(`actions-box` = TRUE),
            selected = treat[1]
        )
        
    })
    
    # We change the assay here, so DE analysis and the visualization are performed in the "RNA" assay.
    single_cell_data_clustered_to_DE_vis <- reactive({
        
        sc_data <- single_cell_data_clustered()
        
        sc_data <- NormalizeData(sc_data,
                                 assay = "RNA",
                                 normalization.method = "LogNormalize",
                                 scale.factor = 10000)
        
        showNotification("Scalling the data of RNA assay",
                         duration = NULL,
                         id = "tab2_m4")
        
        all_genes <- rownames(sc_data@assays$RNA)
        sc_data <- Seurat::ScaleData(sc_data,
                                     assay = "RNA", 
                                     features = all_genes)
        
        on.exit(removeNotification(id = "tab2_m4"), add = TRUE)
        
        DefaultAssay(sc_data) <- "RNA"
        sc_data
    })
    
    # Identification of markers (tab 2)
    markers_tab2 <- eventReactive( input$run_ident_markers_tab2, {
        
        showNotification("Identifing markers or D.E. genes",
                         duration = NULL,
                         id = "tab2_n1")
        
        sc_data <- req( single_cell_data_clustered_to_DE_vis() )
        
        if (input$find_markers_or_DE_tab2 == 0 ) {
            
            showNotification("Identifing markers or D.E. genes",
                             duration = NULL,
                             id = "tab2_n1")
            
            if ( input$find_markers_tab2_opt == 0 ) {
                
                clusters <- as.numeric(as.character(unique(sc_data@meta.data$seurat_clusters)))
                markers_tab2 <- data.frame()
                
                for (i in clusters) {
                    
                    markers <- FindConservedMarkers(sc_data,
                                                    ident.1 = i,
                                                    grouping.var = "treat",
                                                    verbose = T,
                                                    assay = "RNA" )
                    
                    markers$cluster <- i
                    
                    markers_tab2 <- rbind(markers_tab2, markers)
                    
                }
                
                markers_tab2
                
            } else if ( input$find_markers_tab2_opt == 1 ) {
                
                markers_tab2 <- FindConservedMarkers(sc_data,
                                                     ident.1 = req(input$find_markers_clust_id_tab2),
                                                     grouping.var = "treat",
                                                     assay = "RNA" )
                
                markers_tab2$cluster <- input$find_markers_clust_id_tab2
                
            } else if ( input$find_markers_tab2_opt == 2 ) {
                
                markers_tab2 <- FindConservedMarkers(sc_data,
                                                     ident.1 = req(input$find_markers_clust_ID1_tab2),
                                                     ident.2 = req(input$find_markers_clust_ID2_tab2),
                                                     grouping.var = "treat",
                                                     assay = "RNA" )
                
                markers_tab2$cluster <- paste0(input$find_markers_clust_ID1_tab2,
                                               "_vs_",
                                               input$find_markers_clust_ID2_tab2)
                
            }
            
        } else if (input$find_markers_or_DE_tab2 == 1 ) {
            
            showNotification("Identifing markers or D.E. genes",
                             duration = NULL,
                             id = "tab2_n1")
            
            shinyFeedback::feedbackWarning("find_markers_or_DE_tab2_pvalue",
                                           is.na(input$find_markers_or_DE_tab2_pvalue),
                                           "Required value")
            req(input$find_markers_or_DE_tab2_pvalue)
            
            sc_data$celltype.treat <- paste(Idents(sc_data), sc_data$treat, sep = "_")
            Idents(sc_data) <- "celltype.treat"
            
            markers_tab2 <- FindMarkers(sc_data,
                                        ident.1 = paste(req(input$find_markers_or_DE_tab2_cluster),
                                                        req(input$find_markers_or_DE_tab2_treat1),
                                                        sep = "_"),
                                        ident.2 = paste(req(input$find_markers_or_DE_tab2_cluster),
                                                        req(input$find_markers_or_DE_tab2_treat2),
                                                        sep = "_"),
                                        test.use = input$find_markers_tab2_test.use,
                                        assay = "RNA",
                                        verbose = T,
            )
            
            markers_tab2 <- markers_tab2[markers_tab2$p_val_adj < input$find_markers_or_DE_tab2_pvalue, ]
            
            markers_tab2$comparison <- paste("Cluster", input$find_markers_or_DE_tab2_cluster,
                                             input$find_markers_or_DE_tab2_treat1,
                                             "vs",
                                             input$find_markers_or_DE_tab2_treat2,
                                             sep = "_")
            
        }
        
        markers_tab2$geneID <- base::rownames(markers_tab2)
        markers_tab2
        
    })
    
    output$markers_tab2_react <- renderReactable({
        
        markers_tab2 <- req( markers_tab2() )
        markers_tab2 <- markers_tab2[ , c( ncol(markers_tab2), (ncol(markers_tab2)-1), 1: (ncol(markers_tab2) -2) ) ]
        
        on.exit(removeNotification(id = "tab2_n1"), add = TRUE)
        
        my_reactable(markers_tab2)
        
    })
    
    output$download_markers_tab2 <- downloadHandler(
        
        filename = function() {
            
            if (input$find_markers_or_DE_tab2 == 0 ) {
                
                if ( input$find_markers_tab2_opt == 0 ) {
                    
                    paste("list_of_markers_int.dataset_all_clusters", ".csv", sep = "")
                    
                } else if (input$find_markers_tab2_opt == 1) {
                    
                    paste("list_of_markers_int.dataset", "cluster", input$find_markers_clust_id_tab2, ".csv", sep = "")
                    
                } else if (input$find_markers_tab2_opt == 2) {
                    
                    paste("list_of_markers_int.dataset", "cluster", input$find_markers_clust_ID1_tab2, "_vs_" , input$find_markers_clust_ID2_tab2, ".csv", sep = "")
                    
                }
                
            } else if (input$find_markers_or_DE_tab2 == 1 ) {
                
                paste("list_of_markers_int.dataset_DEGs_",
                      input$find_markers_or_DE_tab2_treat1,
                      "_vs_",
                      input$find_markers_or_DE_tab2_treat2,
                      "_cluster_",
                      input$find_markers_or_DE_tab2_cluster,
                      ".csv", sep = "")
                
            }
        },
        content = function(file) {
            
            withProgress(message = "Please wait, preparing the data for download.",
                         value = 0.5, {
                             
                             markers_tab2 <- req( markers_tab2() )
                             markers_tab2 <- markers_tab2[ , c( ncol(markers_tab2), (ncol(markers_tab2)-1), 1: (ncol(markers_tab2) -2) ) ]
                             
                             write.csv(markers_tab2, file, row.names = FALSE)
                             
                         })
        }
        
    )

    features_tab2 <- reactive({
        
        req(input$markers_list_tab2)
        req(input$markers_list_header_opt_tab2)
        
        ext <- tools::file_ext(input$markers_list_tab2$name)
        markers_list_file <- input$markers_list_tab2$datapath
        
        if(input$markers_list_header_opt_tab2 == "") {
            shinyFeedback::feedbackWarning("markers_list_header_opt_tab2",
                                           TRUE,
                                           "Please, inform if the file has a header.")
        }
        
        read_file_markers("tab2_readfile",
                          markers_list_file = markers_list_file,
                          feed_ID ="markers_list_tab2",
                          ext = ext,
                          header_opt = input$markers_list_header_opt_tab2)
    })
    
    observeEvent(input$load_markers_tab2, {
        
        output$marker_group_selec_tab2 = renderUI({
            
            pickerInput_markers_group("features_group_tab2", genes = req( features_tab2() ) )
            
        })
        
        output$marker_genes_selec_tab2 = renderUI({
            
            features_group <- req(input$features_group_tab2)
            id_choice <- req(input$genes_ids_tab2)
            
            pickerInput_markers_genes("selected_genes_tab2",
                                      genes = req( features_tab2() ),
                                      features_group = features_group,
                                      id_choice = id_choice)
        })
        
    })
    
    filt_features_tab2 <- eventReactive(input$run_heatmap_tab2, {
        
        features_f <- req( features_tab2() )
        features_group <- req(input$features_group_tab2)
        
        features_f <- dplyr::filter(features_f,
                                    Group %in% features_group)
        
        if (input$filter_genes_q_tab2 == 0) {
            
            shinyFeedback::feedbackWarning("selected_genes_tab2", is.null(input$selected_genes_tab2), "Please, select one or more genes")
            req(input$selected_genes_tab2)
            
            if (input$genes_ids_tab2 == "ID") {
                
                features_f <- features_f[features_f$GeneID %in% req(input$selected_genes_tab2), ]
                
            } else if (input$genes_ids_tab2 == "name") {
                
                features_f <- features_f[features_f$Name %in% req(input$selected_genes_tab2), ]
                
            }
            
        }
        
        features_f
        
    })
    
    sc_data_av_react_tab2 <- eventReactive(input$run_heatmap_tab2, {
        
        sc_data_av_tab2  <- Seurat::AverageExpression(req( single_cell_data_clustered_to_DE_vis() ),
                                                      assays = "RNA",
                                                      slot = input$slot_selection_heatmap_tab2)
        
        sc_data_av_tab2 <- as.matrix(sc_data_av_tab2[[1]])
        
        on.exit(removeNotification(id = "m7"), add = TRUE)
        sc_data_av_tab2
        
    })
    
    observeEvent(input$run_heatmap_tab2, {
        
        heatmap_n_genes_tab2 <- reactive({
            
            req(filt_features_tab2())
            
            heatmap_genes(features = req( filt_features_tab2() ),
                          sc_data_av = req( sc_data_av_react_tab2()) )
            
        })
        
        heatmap_Height_tab2 <- reactive( 150 + ( 20 * req( heatmap_n_genes_tab2() ) ) )
        
        heat_map_prep_tab2 <- reactive({
            
            features <- req(filt_features_tab2())
            sc_data_av <- req( sc_data_av_react_tab2() )
            features_selec <- as.data.frame(unique(features$GeneID))
            
            if (nrow(features) == 1) {
                
                sc_data_av_feat <- sc_data_av[base::rownames(sc_data_av) %in% features_selec[, 1], ]
                sc_data_av_feat <- as.data.frame(t(sc_data_av_feat))
                base::rownames(sc_data_av_feat) <- features_selec[1, 1]
                
            } else {
                
                sc_data_av_feat <- sc_data_av[base::rownames(sc_data_av) %in% features_selec[, 1], ]
                
            }
            
            sc_data_av_feat
            
        })
        
        output$heat_map_tab2 <- renderPlot({
            
            heat_map_prep_tab2 <- req( heat_map_prep_tab2() )
            
            ComplexHeatmap::Heatmap(heat_map_prep_tab2, border = TRUE,
                                    rect_gp = gpar(col = "white", lwd = 2),
                                    column_title = "Clusters",
                                    column_title_side = "bottom",
                                    name = "Expression",
                                    show_row_dend = T)
            
        }, height = req( heatmap_Height_tab2()) )
        
        output$heat_map_ui_tab2 <- renderUI({
            
            plotOutput("heat_map_tab2", height = req( heatmap_Height_tab2()) )
            
        })
        
        output$marker_to_feature_plot_tab2 = renderUI({
            
            sc_data_av_react <- req( sc_data_av_react_tab2() )
            genes_names <- req( filt_features_tab2() )
            
            genes_names <- genes_names[genes_names$GeneID %in% base::rownames(sc_data_av_react), ]
            
            if (input$genes_ids_tab2 == "ID") {
                
                genes_names <- unique(genes_names$GeneID)
                
            } else if (input$genes_ids_tab2 == "name") {
                
                genes_names <- unique(genes_names$Name)
                
            }
            
            div(class = "option-group",
                shinyWidgets::pickerInput(
                    inputId = "selected_genes_for_feature_plot_tab2",
                    label = "Select the genes to feature plots",
                    choices = sort(as.character(genes_names)),
                    multiple = TRUE,
                    options = list(`actions-box` = TRUE)
                ))
            
        })
        
        output$p8_down <- downloadHandler(
            
            filename = function() {
                paste("Heatmap", ".", input$p8_format, sep = "")
            },
            content = function(file) {
                
                withProgress(message = "Please wait, preparing the data for download.",
                             value = 0.5, {
                                 
                                 shinyFeedback::feedbackWarning("p8_height", is.na(input$p8_height), "Required value")
                                 shinyFeedback::feedbackWarning("p8_width", is.na(input$p8_width), "Required value")
                                 
                                 height <- as.numeric( req( input$p8_height) )
                                 width <- as.numeric( req( input$p8_width) )
                                 res <- as.numeric(input$p8_res)
                                 
                                 heat_map_prep_tab2 <- req( heat_map_prep_tab2() )
                                 
                                 p <- ComplexHeatmap::Heatmap(heat_map_prep_tab2, border = TRUE,
                                                              rect_gp = gpar(col = "white", lwd = 2),
                                                              column_title = "Clusters",
                                                              column_title_side = "bottom",
                                                              name = "Expression",
                                                              show_row_dend = T)
                                 
                                 gb = grid.grabExpr(draw(p))
                                 
                                 ggplot2::ggsave(file,
                                                 gb,
                                                 height=height,
                                                 width=width,
                                                 units="cm",
                                                 dpi=res)
                                 
                             })
            }
            
        )
        
    })
    
    features_selec_tab2 <- eventReactive(input$run_feature_plot_tab2, {
        
        shinyFeedback::feedbackWarning("selected_genes_for_feature_plot_tab2",
                                       is.null(input$selected_genes_for_feature_plot_tab2),
                                       "Please, select one or more genes")
        req(input$selected_genes_for_feature_plot_tab2)
        
        features <- req( filt_features_tab2() )
        
        if (input$genes_ids_tab2 == "ID") {
            
            features_f <- dplyr::filter(features,
                                        GeneID %in% input$selected_genes_for_feature_plot_tab2)
            
        } else if (input$genes_ids_tab2 == "name") {
            
            features_f <- dplyr::filter(features,
                                        Name %in% input$selected_genes_for_feature_plot_tab2)
            
        }
        
        as.character(unique(features_f$GeneID))
        
    })
    
    observeEvent(input$run_feature_plot_tab2, {
        
        output$feature_plot_tab2 <- renderUI({
            
            feat_length <- length(req( features_selec_tab2()) )
            
            plot_output_list <- lapply(1:feat_length, function(i) {
                plotname <- paste("plot1_tab2", i, sep="")
                plotOutput(plotname, height = 300)
            })
            
            # Convert the list to a tagList - this is necessary for the list of items
            # to display properly.
            do.call(tagList, plot_output_list)
        })
        
        output$umap2_tab2 <- renderUI({
            
            feat_length <- length(req( features_selec_tab2()) )
            
            plot_output_list <- lapply(1:feat_length, function(i) {
                plotname <- paste("plot3_tab2", i, sep="")
                plotOutput(plotname, height = 300)
            })
            
            do.call(tagList, plot_output_list)
        })
        
        output$run_vln_plot_tab2 <- renderUI({
            
            feat_length <- length(req( features_selec_tab2()) )
            
            plot_output_list <- lapply(1:feat_length, function(i) {
                plotname <- paste("plot4_tab2", i, sep="")
                plotOutput(plotname, height = 300)
            })
            
            do.call(tagList, plot_output_list)
        })
        
        output$run_dot_plot_tab2 <- renderUI({
            
            feat_length <- length(req( features_selec_tab2()) )
            
            plot_output_list <- lapply(1:feat_length, function(i) {
                plotname <- paste("plot5_tab2", i, sep="")
                plotOutput(plotname, height = 300)
            })
            
            do.call(tagList, plot_output_list)
        })
        
        sc_data_tab2 <- req( single_cell_data_clustered_to_DE_vis() )
        assay_id <- "RNA"
        
        for (i in 1:length(req( features_selec_tab2())) ) {
            
            features <- req( features_selec_tab2() )
            local({
                
                my_i <- i
                
                plotname <- paste("plot1_tab2", my_i, sep="")
                output[[plotname]] <- renderPlot({
                    
                    p_list <- FeaturePlotSingle(sc_data_tab2,
                                                feature = features[my_i],
                                                metadata_column = "treat",
                                                assay_id = assay_id,
                                                pt.size = 0.2,
                                                order = TRUE,
                                                reduction = "umap",
                                                label = T)
                    
                    groups <- unique(sc_data_tab2@meta.data[, "treat"])
                    
                    if ( length(groups) > 3 ) {
                        
                        wrap_plots(p_list,
                                   guides = 'collect',
                                   ncol = length(groups) ) +
                            plot_annotation(title = paste("Gene name:",
                                                          features[my_i]),
                                            theme = theme(plot.title = element_text(size = 16,
                                                                                    face = "bold",
                                                                                    hjust = 0.5)) )
                    } else {
                        
                        wrap_plots(p_list,
                                   guides = 'collect',
                                   ncol = length(groups),
                                   ceiling( (length(groups) / 3  + 1 ) ) ) +
                            plot_annotation(title = paste("Gene name:",
                                                          features[my_i]),
                                            theme = theme(plot.title = element_text(size = 16,
                                                                                    face = "bold",
                                                                                    hjust = 0.5)) )
                    }
                    
                })
                
                plotname <- paste("plot3_tab2", my_i, sep="")
                output[[plotname]] <- renderPlot({
                    
                    Seurat::DimPlot(req( single_cell_data_clustered_to_DE_vis() ), reduction = "umap", label = T, pt.size = .1)
                    
                })
                
                plotname <- paste("plot4_tab2", my_i, sep="")
                output[[plotname]] <- renderPlot({
                    
                    Seurat::VlnPlot(sc_data_tab2,
                                    features = features[my_i],
                                    split.by = "treat",
                                    assay = "RNA"
                    )
                    
                })
                
                plotname <- paste("plot5_tab2", my_i, sep="")
                output[[plotname]] <- renderPlot({
                    
                    Seurat::DotPlot(sc_data_tab2,
                                    features = features[my_i],
                                    cols = c("lightgrey", "red"),
                                    split.by = "treat",
                                    assay = "RNA" ) +
                        scale_x_discrete(position = "top") +
                        xlab("")+
                        ylab("") +
                        ggtitle("") +
                        theme(legend.position = "",
                              axis.title.x=element_blank(),
                              axis.title.y=element_blank(),
                              axis.text.y = element_text( size = rel(1) ),
                              axis.text.x = element_text( angle = 90) ) +
                        coord_flip()
                    
                })
                
            })
            
        }
        showNotification("Generating additional plots",
                         duration = 30,
                         id = "m8")
        
        output$select_genes_add_plot_to_down_tab2_ui <- renderUI({
            
            filt_features <- req(input$selected_genes_for_feature_plot_tab2)
            
            div(class = "option-group",
                shinyWidgets::pickerInput(
                    inputId = "select_genes_add_plot_to_down_tab2",
                    label = "Select the genes that you want to download",
                    choices = sort(filt_features),
                    multiple = TRUE,
                    options = list(`actions-box` = TRUE)
                ))
            
        })
    })
    
    observeEvent(input$start_down_add_plots_tab2, {
        
        on.exit(removeNotification(id = "m8"), add = TRUE)
        
        shinyFeedback::feedbackWarning("select_genes_add_plot_to_down_tab2",
                                       is.null(input$select_genes_add_plot_to_down_tab2),
                                       "Please, select one or more genes")
        genes <- req(input$select_genes_add_plot_to_down_tab2)
        
        withProgress(message = "Please wait, preparing the data for download.",
                     value = 0.5, {
                         
                         if (file.exists("./images") == F ) {
                             dir.create('./images')
                         }
                         
                         # Creates new folder to keep the plots organized. Adding date and time helps to avoid overwriting the data
                         path_new <- paste0("./images/integrated_sample_plots", format(Sys.time(),'_%Y-%m-%d__%H%M%S'))
                         
                         dir.create(path_new)
                         dir.create( paste0(path_new,"/feature_plots_combined_samples") )
                         dir.create( paste0(path_new,"/violin_plots_combined_samples") )
                         dir.create( paste0(path_new, "/dot_plots_combined_samples") )
                         
                         dir.create( paste0(path_new,"/feature_plots") )
                         dir.create( paste0(path_new,"/violin_plots") )
                         dir.create( paste0(path_new, "/dot_plots") )
                         
                         sc_data <- req( single_cell_data_clustered_to_DE_vis() )
                         assay_id <- "RNA"
                         
                         for( i in 1:length(genes) ){
                             
                             minimal <- min( sc_data[[ assay_id ]]@data[ genes[i], ] )
                             maximal <- max( sc_data[[ assay_id ]]@data[ genes[i], ] )
                             
                             p <- suppressMessages(Seurat::FeaturePlot(sc_data,
                                                                       features = genes[i],
                                                                       slot = input$slot_selection_feature_plot_tab2,
                                                                       reduction = "umap") +
                                                       scale_colour_gradient2(limits=c(minimal, maximal),
                                                                              midpoint = maximal / 2,
                                                                              low = "gray80",
                                                                              mid = "gold",
                                                                              high = "red"))
                             
                             file <- paste0(path_new, "/feature_plots_combined_samples/", genes[i], ".", req(input$add_p_tab2_feat_format) )
                             
                             ggplot2::ggsave(file,
                                             p,
                                             height= req(input$add_p_tab2_feat_height),
                                             width=req(input$add_p_tab2_feat_width),
                                             units="cm",
                                             dpi=as.numeric(req(input$add_p_tab2_feat_res)))
                             
                             file <- paste0(path_new, "/violin_plots_combined_samples/", genes[i], ".", req(input$add_p_tab2_violin_format) )
                             
                             p <- Seurat::VlnPlot(sc_data,
                                                  features = genes[i],
                                                  slot = req(input$slot_selection_feature_plot_tab2) )
                             
                             ggplot2::ggsave(file,
                                             p,
                                             height=req(input$add_p_tab2_violin_height),
                                             width=req(input$add_p_tab2_violin_width),
                                             units="cm",
                                             dpi=as.numeric( req(input$add_p_tab2_violin_res) ))
                             
                             file <- paste0(path_new, "/dot_plots_combined_samples/", genes[i], ".", req(input$add_p_tab2_dot_format) )
                             
                             p <- Seurat::DotPlot(sc_data,
                                                  features = genes[i],
                                                  cols = c("lightgrey", "red")) +
                                 scale_x_discrete(position = "top") +
                                 xlab("")+
                                 ylab("") +
                                 ggtitle("") +
                                 theme(legend.position = "",
                                       axis.title.x=element_blank(),
                                       axis.title.y=element_blank(),
                                       axis.text.y = element_text( size = rel(1) ),
                                       axis.text.x = element_text( angle = 90) ) +
                                 coord_flip()
                             
                             ggplot2::ggsave(file,
                                             p,
                                             height= req(input$add_p_tab2_dot_height),
                                             width= req(input$add_p_tab2_dot_width),
                                             units="cm",
                                             dpi=as.numeric( req(input$add_p_tab2_dot_res) ))
                             
                             p_list <- FeaturePlotSingle(sc_data,
                                                         feature = genes[i],
                                                         metadata_column = "treat",
                                                         assay_id = assay_id,
                                                         pt.size = 0.2,
                                                         order = TRUE,
                                                         reduction = "umap",
                                                         label = T)
                             
                             groups <- unique(sc_data@meta.data[, "treat"])
                             
                             if ( length(groups) > 3 ) {
                                 
                                 p2 <-  wrap_plots(p_list,
                                                   guides = 'collect',
                                                   ncol = length(groups) ) +
                                     plot_annotation(title = paste("Gene name:",
                                                                   genes[i]),
                                                     theme = theme(plot.title = element_text(size = 16,
                                                                                             face = "bold",
                                                                                             hjust = 0.5)) )
                             } else {
                                 
                                 p2 <- wrap_plots(p_list,
                                                  guides = 'collect',
                                                  ncol = length(groups),
                                                  ceiling( (length(groups) / 3  + 1 ) ) ) +
                                     plot_annotation(title = paste("Gene name:",
                                                                   genes[i]),
                                                     theme = theme(plot.title = element_text(size = 16,
                                                                                             face = "bold",
                                                                                             hjust = 0.5)) )
                             }
                             
                             file <- paste0(path_new, "/feature_plots/", genes[i], ".", req(input$add_p_tab2_feat_format) )
                             
                             ggplot2::ggsave(file,
                                             p2,
                                             height = req(input$add_p_tab2_feat_height),
                                             width = (req(input$add_p_tab2_feat_width) * length(groups)),
                                             units = "cm",
                                             dpi = as.numeric(req(input$add_p_tab2_feat_res)))
                             
                             file <- paste0(path_new, "/violin_plots/", genes[i], ".", req(input$add_p_tab2_violin_format) )
                             
                             p2 <- Seurat::VlnPlot(sc_data,
                                                   features = genes[i],
                                                   split.by = "treat",
                                                   assay = "RNA",
                                                   slot = req(input$slot_selection_feature_plot_tab2) )
                             
                             ggplot2::ggsave(file,
                                             p2,
                                             height = req(input$add_p_tab2_violin_height),
                                             width = ( req(input$add_p_tab2_violin_width) * length(groups) ),
                                             units = "cm",
                                             dpi = as.numeric( req(input$add_p_tab2_violin_res) ))

                             file <- paste0(path_new, "/dot_plots/", genes[i], ".", req(input$add_p_tab2_dot_format) )
                             
                             p2 <- Seurat::DotPlot(sc_data,
                                                   features = genes[i],
                                                   cols = c("lightgrey", "red"),
                                                   split.by = "treat",
                                                   assay = "RNA" ) +
                                 scale_x_discrete(position = "top") +
                                 xlab("")+
                                 ylab("") +
                                 ggtitle("") +
                                 theme(legend.position = "",
                                       axis.title.x=element_blank(),
                                       axis.title.y=element_blank(),
                                       axis.text.y = element_text( size = rel(1) ),
                                       axis.text.x = element_text( angle = 90) ) +
                                 coord_flip()
                             
                             ggplot2::ggsave(file,
                                             p2,
                                             height= req(input$add_p_tab2_dot_height),
                                             width= ( req(input$add_p_tab2_dot_width) * length(groups) ),
                                             units="cm",
                                             dpi=as.numeric( req(input$add_p_tab2_dot_res) ))
                         }
                         
                     })
        
        showNotification("All plots were downloaded!",
                         duration = 15,
                         id = "")
        
    })
    
    #####################################
    ### Tab 3 - trajectory inference ####
    #####################################
    
    output$load_integrated_ui_tab3 <- renderUI({
        
        rds_list <- list.files('./RDS_files/', pattern = "*.rds")
        
        div(class = "option-group",
            shinyWidgets::pickerInput(
                inputId = "load_integrated_tab3",
                label = "Select the file containing the clustered data",
                choices = sort(rds_list),
                multiple = FALSE,
                options = list(`actions-box` = TRUE)
            ))
    })
    
    sc_data_traj_inf <- eventReactive(input$run_ti_model, {
        
        #if ( input$rds_location_tab3 == 0 ) {
        
        showNotification("Loading the clustered data",
                         id = "tab3_m1",
                         duration = NULL)
        
        load_rds <- req(input$load_integrated_tab3)
        sc_data <- readRDS( paste0("./RDS_files/", load_rds) )
        
        # } else if ( input$rds_location_tab3 == 1 ) {
        #
        #     ext <- tools::file_ext(input$file_input_rds_tab3$name)
        #     config_input_file <- input$file_input_rds_tab3$datapath
        #
        #     req(input$file_input_rds_tab3$datapath)
        #
        #     showNotification("Loading the clustered data",
        #                      id = "tab3_m1",
        #                      duration = NULL)
        #
        #     sc_data <- readRDS(config_input_file)
        #
        # }
        
        
        on.exit(removeNotification(id = "tab3_m1"), add = TRUE)
        
        validate(need("seurat_clusters" %in% colnames(sc_data@meta.data),
                      "Clustering was not detected. Was the data clustered in the other tabs?", ""))
        sc_data
        
    } )
    
    expressed_genes <- reactive({
        
        req(sc_data_traj_inf())
        
        sce <- as.SingleCellExperiment(sc_data_traj_inf())
        sce@assays@data@listData$counts@Dimnames[[1]]
        
    })
    
    # extract the expression matrix from the Seurat object
    object_expression <- reactive({
        
        req(sc_data_traj_inf())
        sing_cell_data <- sc_data_traj_inf()
        
        Matrix::t(as(as.matrix(sing_cell_data@assays$RNA@data), 'sparseMatrix'))
        
    })
    
    object_counts <- reactive({
        
        req(sc_data_traj_inf())
        sing_cell_data <- sc_data_traj_inf()
        
        Matrix::t(as(as.matrix(sing_cell_data@assays$RNA@counts), 'sparseMatrix'))
        
    })
    
    dataset_inf <- reactive({
        
        wrap_expression(
            counts = req( object_counts() ),
            expression = req( object_expression() )
        )
        
    })
    
    # Extracts meta data to fill the dyno object
    sc_meta <- reactive({
        
        req( sc_data_traj_inf() )
        sing_cell_data <- sc_data_traj_inf()
        
        sing_cell_data[[]] %>%
            mutate(cells = rownames(.))
        
    })
    
    ## Clustering
    sc_meta_cluster <- reactive({
        
        sc_meta <- req( sc_meta() )
        
        data.frame(cell_id = sc_meta$cells,
                   group_id = sc_meta$seurat_clusters)
        
    })
    
    sc_cells_time_vec <- reactive ({
        
        sc_meta <- req( sc_meta() )
        
        sc_cells_time_vec <- as.character(sc_meta$treat)
        names(sc_cells_time_vec) <- sc_meta$cells
        
        sc_cells_time_vec
    })
    
    output$ti_methods_list_ui <- renderUI ({
        
        ti_methods_list <- dynmethods::methods$method_id
        
        div(class = "option-group",
            shinyWidgets::pickerInput(
                inputId = "ti_methods_list",
                label = "Select the dynverse method to execute",
                choices = sort(as.character(ti_methods_list)),
                multiple = FALSE,
                selected = "slingshot",
                options = list(`actions-box` = TRUE)
            ))
        
    })
    
    model <- eventReactive(input$run_ti_model, { # change it to run model
        
        
        req(sc_data_traj_inf())
        sing_cell_data <- sc_data_traj_inf()
        
        showNotification("Running trajectory inference model",
                         id = "tab3_m2",
                         duration = NULL)
        
        ## Extracts the dimension reduction info
        dimred <- sing_cell_data@reductions$pca@cell.embeddings
        
        sc_meta <- req( sc_meta() )
        sc_meta_cluster <-  req( sc_meta_cluster() )
        object_expression <- req( object_expression() )
        object_counts <- req( object_counts() )
        
        # test if the user set the initial cluster
        if ( !is.na(input$traj_init_clusters) ) {
            
            start_cells <- sc_meta[sc_meta$seurat_clusters == input$traj_init_clusters, ]
            
            start_cells <- as.character(start_cells$cells)
            
        } else {
            
            start_cells <- NULL
            
        }
        
        # test if the user set the end cluster
        if (!is.na(input$traj_end_clusters)) {
            
            end_cells <- sc_meta[sc_meta$seurat_clusters == input$traj_end_clusters, ]
            
            end_cells <- as.character(end_cells$cells)
            
        } else {
            
            end_cells <- NULL
            
        }
        
        if ( input$ti_select_models == 0 ) { # slingshot locally
            
            if (input$ti_sample_number == 1) { # there are multiple samples
                
                priors <- list(
                    start_id = start_cells,
                    end_id = end_cells,
                    groups_id = sc_meta_cluster,
                    dimred = dimred
                )
                
            } else if (input$ti_sample_number == 0) { # there are only one sample
                
                priors <- list(
                    start_id = start_cells,
                    end_id = end_cells,
                    groups_id = sc_meta_cluster,
                    dimred = dimred
                )
                
            }
            
            sds <- run_fun_slingshot(expression = object_expression, priors = priors, verbose = T)
            
            on.exit(removeNotification(id = "tab3_m2"), add = TRUE)
            
        } else if ( input$ti_select_models == 1 ) { # dynverse models - depends on docker
            
            dataset <- wrap_expression(
                counts = object_counts,
                expression = object_expression
            )
            
            # fetch newest version of the method
            method_id <- paste0("dynverse/ti_", req(input$ti_methods_list), ":latest")
            methods_selected <- create_ti_method_container(method_id)
            
            if ( input$ti_sample_number == 1 ) { # there are multiple samples
                
                
                dataset <- add_prior_information(
                    
                    dataset,
                    start_id = start_cells,
                    end_id = end_cells,
                    groups_id = sc_meta_cluster,
                    dimred = dimred
                )
                
                sds <- infer_trajectory(dataset,
                                        req( methods_selected() ),
                                        verbose = T,
                                        give_priors = c("start_id",
                                                        "end_id",
                                                        "groups_id",
                                                        "dimred"))
                
            } else if ( input$ti_sample_number == 0 ) { # there only one samples
                
                dataset <- add_prior_information(
                    
                    dataset,
                    start_id = start_cells,
                    end_id = end_cells,
                    groups_id = sc_meta_cluster,
                    dimred = dimred
                )
                
                sds <- infer_trajectory(dataset,
                                        req( methods_selected() ),
                                        verbose = T,
                                        give_priors = c("start_id",
                                                        "end_id",
                                                        "groups_id",
                                                        "dimred"))
                
            }
            
        }
        
        sds
        
    })
    
    output$ti_order <- renderPlot({
        
        model <- req( model() )
        sc_meta_cluster <- req( sc_meta_cluster() )
        
        if ( input$ti_sample_number == 0 ) {
            
            plot_dimred(model,
                        grouping = sc_meta_cluster,
                        color_density = "grouping") +
                ggtitle("Cell grouping")#,
            
        } else if (input$ti_sample_number == 1) {
            
            if ( input$ti_graphs_color_choice == 1 ) {
                
                plot_dimred(model,
                            grouping = sc_meta_cluster,
                            color_density = "grouping") +
                    ggtitle("Cell grouping")#,
                
            } else if (input$ti_graphs_color_choice == 0 ) {
                
                sc_cells_time_vec <- req( sc_cells_time_vec() )
                
                plot_dimred(model,
                            grouping = sc_cells_time_vec,
                            color_density = "grouping") +
                    ggtitle("Cell grouping")#,
                
            }
        }
        
    })
    
    output$ti_traject <- renderPlot({
        
        model <- req( model() )
        sc_meta_cluster <- req( sc_meta_cluster() )
        
        if ( input$ti_sample_number == 0 ) {
            
            plot_dendro(model,
                        grouping = sc_meta_cluster) +
                ggtitle("Trajectory")#,
            
        } else if (input$ti_sample_number == 1) {
            
            if ( input$ti_graphs_color_choice == 1 ) {
                
                plot_dendro(model,
                            grouping = sc_meta_cluster) +
                    ggtitle("Trajectory")#,
                
            } else if (input$ti_graphs_color_choice == 0 ) {
                sc_cells_time_vec <- req( sc_cells_time_vec() )
                
                plot_dendro(model,
                            grouping = sc_cells_time_vec) +
                    ggtitle("Trajectory")#,
                
            }
        }
        
        
    })
    
    output$ti_graph <- renderPlot({
        
        model <- req( model() )
        sc_meta_cluster <- req( sc_meta_cluster() )
        
        if ( input$ti_sample_number == 0 ) {
            
            plot_graph(model,
                       grouping = sc_meta_cluster,
                       expression_source = dataset) +
                ggtitle("Trajectory represented as a graph")
            
        } else if (input$ti_sample_number == 1) {
            
            if ( input$ti_graphs_color_choice == 1 ) {
                
                plot_graph(model,
                           grouping = sc_meta_cluster,
                           expression_source = dataset) +
                    ggtitle("Trajectory represented as a graph")
                
            } else if (input$ti_graphs_color_choice == 0 ) {
                sc_cells_time_vec <- req( sc_cells_time_vec() )
                
                plot_graph(model,
                           grouping = sc_cells_time_vec,
                           expression_source = dataset) +
                    ggtitle("Trajectory represented as a graph")
                
            }
        }
        
    })
    
    output$p9_down <- downloadHandler(
        
        filename = function() {
            
            if ( input$p9_down_opt == 0 ) {
                
                paste("Trajectory_Dimension_Reduction", ".", input$p9_format, sep = "")
                
            } else if (input$p9_down_opt == 1) {
                
                paste("Trajectory_Dendrogram", ".", input$p9_format, sep = "")
                
            } else if (input$p9_down_opt == 2) {
                
                paste("Trajectory_Graph", ".", input$p9_format, sep = "")
                
            }
            
        },
        content = function(file) {
            
            withProgress(message = "Please wait, preparing the data for download.",
                         value = 0.5, {
                             
                             shinyFeedback::feedbackWarning("p9_height", is.na(input$p9_height), "Required value")
                             shinyFeedback::feedbackWarning("p9_width", is.na(input$p9_width), "Required value")
                             
                             height <- as.numeric( req( input$p9_height) )
                             width <- as.numeric( req( input$p9_width) )
                             res <- as.numeric(input$p9_res)
                             
                             if ( input$p9_down_opt == 0 ) {
                                 
                                 model <- req( model() )
                                 sc_meta_cluster <- req( sc_meta_cluster() )
                                 
                                 if ( input$ti_sample_number == 0 ) {
                                     
                                     p <- plot_dimred(model,
                                                      grouping = sc_meta_cluster,
                                                      color_density = "grouping") +
                                         ggtitle("Cell grouping")#,
                                     
                                 } else if (input$ti_sample_number == 1) {
                                     
                                     if ( input$ti_graphs_color_choice == 1 ) {
                                         
                                         p <-  plot_dimred(model,
                                                           grouping = sc_meta_cluster,
                                                           color_density = "grouping") +
                                             ggtitle("Cell grouping")#,
                                         
                                     } else if (input$ti_graphs_color_choice == 0 ) {
                                         
                                         sc_cells_time_vec <- req( sc_cells_time_vec() )
                                         
                                         p <-  plot_dimred(model,
                                                           grouping = sc_cells_time_vec,
                                                           color_density = "grouping") +
                                             ggtitle("Cell grouping")#,
                                         
                                     }
                                 }
                                 
                             } else if ( input$p9_down_opt == 1 ) {
                                 
                                 model <- req( model() )
                                 sc_meta_cluster <- req( sc_meta_cluster() )
                                 
                                 if ( input$ti_sample_number == 0 ) {
                                     
                                     p <- plot_dendro(model,
                                                      grouping = sc_meta_cluster) +
                                         ggtitle("Trajectory")#,
                                     
                                 } else if (input$ti_sample_number == 1) {
                                     
                                     if ( input$ti_graphs_color_choice == 1 ) {
                                         
                                         p <- plot_dendro(model,
                                                          grouping = sc_meta_cluster) +
                                             ggtitle("Trajectory")#,
                                         
                                     } else if (input$ti_graphs_color_choice == 0 ) {
                                         sc_cells_time_vec <- req( sc_cells_time_vec() )
                                         
                                         p <-  plot_dendro(model,
                                                           grouping = sc_cells_time_vec) +
                                             ggtitle("Trajectory")#,
                                         
                                     }
                                 }
                                 
                             } else if ( input$p9_down_opt == 2 ) {
                                 
                                 model <- req( model() )
                                 sc_meta_cluster <- req( sc_meta_cluster() )
                                 
                                 if ( input$ti_sample_number == 0 ) {
                                     
                                     p <- plot_graph(model,
                                                     grouping = sc_meta_cluster,
                                                     expression_source = dataset) +
                                         ggtitle("Trajectory represented as a graph")
                                     
                                 } else if (input$ti_sample_number == 1) {
                                     
                                     if ( input$ti_graphs_color_choice == 1 ) {
                                         
                                         p <- plot_graph(model,
                                                         grouping = sc_meta_cluster,
                                                         expression_source = dataset) +
                                             ggtitle("Trajectory represented as a graph")
                                         
                                     } else if (input$ti_graphs_color_choice == 0 ) {
                                         
                                         sc_cells_time_vec <- req( sc_cells_time_vec() )
                                         
                                         p <- plot_graph(model,
                                                         grouping = sc_cells_time_vec,
                                                         expression_source = dataset) +
                                             ggtitle("Trajectory represented as a graph")
                                         
                                     }
                                 }
                                 
                             }
                             
                             ggplot2::ggsave(file,
                                             p,
                                             height=height,
                                             width=width,
                                             units="cm",
                                             dpi=res)
                             
                         })
        }
        
    )
    
    ## Print the lineages that are been show in the plots
    
    output$lineages <- renderPrint({
        
        sds <- req( model() )
        sds$lineages
        
    })
    
    ## Gene expression
    features_tab3 <- reactive({
        
        req(input$markers_list_tab3)
        ext <- tools::file_ext(input$markers_list_tab3$name)
        markers_list_file <- input$markers_list_tab3$datapath
        
        if(input$markers_list_header_opt_tab3 == "") {
            shinyFeedback::feedbackWarning("markers_list_header_opt_tab3",
                                           TRUE,
                                           "Please, inform if the file has a header.")
        }
        req(input$markers_list_header_opt_tab3)
        
        read_file_markers("tab3_readfile",
                          markers_list_file = markers_list_file,
                          feed_ID ="markers_list_tab3",
                          ext = ext,
                          header_opt = input$markers_list_header_opt_tab3)
    })
    
    # Load the file and offers the parameters for heatmap
    observeEvent(input$load_markers_tab3, {
        
        # Painel that will apper after loading the list of markers containing the filter options
        output$marker_group_selec_tab3 = renderUI({
            
            pickerInput_markers_group("features_group_tab3", genes = req( features_tab3()) )
            
        })
        
        output$marker_genes_selec_tab3 = renderUI({
            
            features_group <- req(input$features_group_tab3)
            id_choice <- req(input$genes_ids_tab3)
            
            pickerInput_markers_genes("selected_genes_tab3",
                                      genes = req( features_tab3() ),
                                      features_group = features_group,
                                      id_choice = id_choice)
        })
        
    })
    
    # Filter accordly with the parameters selected by the user
    filt_features_tab3 <- eventReactive(input$run_heatmap_tab3, {
        
        features_f <- req( features_tab3() )
        features_group <- req(input$features_group_tab3)
        
        features_f <- dplyr::filter(features_f,
                                    Group %in% features_group)
        
        if (input$filter_genes_q_tab3 == 0) {
            
            shinyFeedback::feedbackWarning("selected_genes_tab3",
                                           is.null(input$selected_genes_tab3),
                                           "Please, select one or more genes")
            req(input$selected_genes_tab3)
            
            if (input$genes_ids_tab3 == "ID") {
                
                features_f <- features_f[features_f$GeneID %in% input$selected_genes_tab3, ]
                
            } else if (input$genes_ids_tab3 == "name") {
                
                features_f <- features_f[features_f$Name %in% input$selected_genes_tab3, ]
                
            }
            
        }
        
        features_f
        
    })
    
    dynverse_genes_list <- eventReactive(input$dynverse_def_imp_genes, {
        
        showNotification("Defining the most relevant genes",
                         id = "tab3_m5",
                         duration = NULL)
        
        if ( input$dynverse_opt == 0 ) { # global
            
            feat_importances <- dynfeature::calculate_overall_feature_importance(req( model() ),
                                                                                 expression_source = req( dataset_inf()))
            
            
        } else if ( input$dynverse_opt == 1 ) { # branch/lineage
            
            feat_importances <- dynfeature::calculate_branch_feature_importance(req( model() ),
                                                                                expression_source = req( dataset_inf()) )
            
        } else if ( input$dynverse_opt == 2 ) { # bifurcation
            
            if ( is.na(input$branching_milestone) ) {
                
                feat_importances <- calculate_branching_point_feature_importance(req( model() ),
                                                                                 expression_source = req( dataset_inf()) )
                
            } else {
                
                feat_importances <- calculate_branching_point_feature_importance(req( model() ),
                                                                                 milestones_oi = input$branching_milestone,
                                                                                 expression_source = req( dataset_inf()))
                
            }
            
        }
        
        on.exit(removeNotification(id = "tab3_m5"), add = TRUE)
        
        showNotification("Select the number of genes to show in the heatmap and trajectory plots",
                         id = "tab3_m6",
                         duration = 30)
        
        feat_importances
        
    })
    
    dynverse_genes_list_filt <- eventReactive( c(input$dynverse_def_imp_genes, input$run_heatmap_tab3_dynverse), {
        
        shinyFeedback::feedbackWarning("dynverse_n_genes",
                                       is.na(input$dynverse_n_genes),
                                       "Required value")
        req(input$dynverse_n_genes)
        
        features_imp <- req( dynverse_genes_list() )
        
        if ( input$dynverse_opt == 0 ) { # global
            
            features <- features_imp %>%
                top_n(input$dynverse_n_genes, importance)
            
        } else if ( input$dynverse_opt == 1 ) { # branch/lineage
            
            if ( is.na(req(input$dynverse_branch_from)) ) {
                
                features <- features_imp %>%
                    filter(to == req(input$dynverse_branch_to)) %>%
                    top_n(input$dynverse_n_genes, importance)
                
            } else if ( is.na(req(input$dynverse_branch_to)) ) {
                
                features <- features_imp %>%
                    filter( from  == req(input$dynverse_branch_from)) %>%
                    top_n(input$dynverse_n_genes, importance)
                
            } else {
                
                features <- features_imp %>%
                    filter( from  == req(input$dynverse_branch_from) & to == req(input$dynverse_branch_to)) %>%
                    top_n(input$dynverse_n_genes, importance)
                
            }
            
        } else if ( input$dynverse_opt == 2 ) { # bifurcation
            
            features <- features_imp %>%
                top_n(input$dynverse_n_genes, importance)
            
        }
        
        # features <- features$b
        
        features
    })
    
    ##################################################
    ### Expression plots - Using markers as input ####
    ##################################################
    
    observeEvent( c(input$dynverse_def_imp_genes, input$run_heatmap_tab3), {
        
        # This calculate the height of the plot based on the n of genes
        heatmap_n_genes_tab3 <- reactive({
            
            # Test what input to use for the heatmap
            if (input$ti_expre_opt == 0) {
                
                features <- req( filt_features_tab3())
                features <- unique(features$GeneID)
                
            } else if (input$ti_expre_opt == 1) {
                
                features <- req( dynverse_genes_list_filt() )
                features <- features$feature_id
                
            } else if (input$ti_expre_opt == 2) {
                
                features <- req( tradseq_genes_list() )
                
            }
            
            as.numeric(length(features))
        })
        heatmap_Height_tab3 <- reactive( 250 + ( 20 * req( heatmap_n_genes_tab3() ) ) )
        
        # Generates the heatmap plot
        output$heat_map_tab3 <- renderPlot({
            
            showNotification("Generating the heatmap",
                             id = "tab3_m4",
                             duration = NULL)
            
            # Test what input to use for the heatmap
            if (input$ti_expre_opt == 0) {
                
                features <- req( filt_features_tab3() )
                features <- unique(features$GeneID)
                
            } else if (input$ti_expre_opt == 1) {
                
                features <- req( dynverse_genes_list_filt() )
                features <- features$feature_id
                
            } else if (input$ti_expre_opt == 2) {
                
                features <- req( tradseq_genes_list() )
                
            }
            
            features <- as.character(features)
            
            ### Drawing heat map ###
            heatmap <- plot_heatmap(
                req( model() ),
                expression_source = req( dataset_inf() ),
                grouping = req( sc_meta_cluster() ),
                features_oi = features
            )
            
            on.exit(removeNotification(id = "tab3_m4"), add = TRUE)
            
            heatmap
            
        }, height = req( heatmap_Height_tab3() ) )
        
        output$heat_map_ui_tab3 <- renderUI({
            
            plotOutput( "heat_map_tab3", height = req( heatmap_Height_tab3() ))
            
        })
        
        output$marker_to_feature_plot_tab3 = renderUI({
            
            # Test what input to use for the heatmap
            if (input$ti_expre_opt == 0) {
                
                features <- req( filt_features_tab3() )
                features <- unique(features$GeneID)
                
            } else if (input$ti_expre_opt == 1) {
                
                features <- req( dynverse_genes_list_filt() )
                features <- features$feature_id
                
            } else if (input$ti_expre_opt == 2) {
                
                features <- req( tradseq_genes_list() )
                
            }
            
            div(class = "option-group",
                shinyWidgets::pickerInput(
                    inputId = "selected_genes_for_feature_plot_tab3",
                    label = "Select the genes to feature plots",
                    choices = sort(as.character(features)),
                    multiple = TRUE,
                    options = list(`actions-box` = TRUE)
                ))
            
        })
        
        output$p10_down <- downloadHandler(
            
            filename = function() {
                paste("Heatmap_int.dataset", ".", input$p10_format, sep = "")
            },
            content = function(file) {
                
                withProgress(message = "Please wait, preparing the data for download.",
                             value = 0.5, {
                                 
                                 shinyFeedback::feedbackWarning("p10_height", is.na(input$p10_height), "Required value")
                                 shinyFeedback::feedbackWarning("p10_width", is.na(input$p10_width), "Required value")
                                 
                                 height <- as.numeric( req( input$p10_height) )
                                 width <- as.numeric( req( input$p10_width) )
                                 res <- as.numeric(input$p10_res)
                                 
                                 # Test what input to use for the heatmap
                                 if (input$ti_expre_opt == 0) {
                                     
                                     features <- req( filt_features_tab3() )
                                     features <- unique(features$GeneID)
                                     
                                 } else if (input$ti_expre_opt == 1) {
                                     
                                     features <- req( dynverse_genes_list_filt() )
                                     features <- features$feature_id
                                     
                                 } else if (input$ti_expre_opt == 2) {
                                     
                                     features <- req( tradseq_genes_list() )
                                     
                                 }
                                 
                                 #    features <- paste0("gene:", as.character(features))
                                 features <- as.character(features)
                                 
                                 ### Drawing heat map ###
                                 p <- plot_heatmap(
                                     req( model() ),
                                     expression_source = req( dataset_inf() ),
                                     grouping = req( sc_meta_cluster() ),
                                     features_oi = features
                                 )
                                 
                                 ggplot2::ggsave(file,
                                                 p,
                                                 height=height,
                                                 width=width,
                                                 units="cm",
                                                 dpi=res)
                                 
                             })
            }
            
        )
        
    })
    
    #######3 >>>>>>>>>sc_meta_cluster
    output$dynverse_branch_from_ui <- renderUI({
        
        sc_meta_cluster <- req( sc_meta_cluster() )
        clust <- unique(as.character(sc_meta_cluster$group_id))
        
        div(class = "option-group",
            shinyWidgets::pickerInput(
                inputId = "dynverse_branch_from",
                label = "Set the start branch (from)",
                choices = sort(as.character(clust)),
                multiple = F,
                options = list(`actions-box` = TRUE)
            ))
        
    })
    
    output$dynverse_branch_to_ui <- renderUI({
        
        sc_meta_cluster <- req( sc_meta_cluster() )
        
        clust <- unique(as.character(sc_meta_cluster$group_id))
        clust <- clust[!clust %in% req(input$dynverse_branch_from)]
        
        div(class = "option-group",
            shinyWidgets::pickerInput(
                inputId = "dynverse_branch_to",
                label = "Set the start branch (to)",
                choices = sort(as.character(clust)),
                multiple = F,
                options = list(`actions-box` = TRUE)
            ))
        
    })
    
    observeEvent(input$run_feature_plot_tab3, {
        
        showNotification("Generating the expression plots",
                         id = "tab3_m7",
                         duration = 30)
        
        output$ti_order_express <- renderUI({
            
            shinyFeedback::feedbackWarning("selected_genes_for_feature_plot_tab3",
                                           is.null(input$selected_genes_for_feature_plot_tab3),
                                           "Please, select one or more genes")
            req(input$selected_genes_for_feature_plot_tab3)
            
            feat_length <- length(req(input$selected_genes_for_feature_plot_tab3))
            
            plot_output_list <- lapply(1:feat_length, function(i) {
                plotname <- paste("ti_order_exp", i, sep="")
                plotOutput(plotname, height = 300)
            })
            
            # Convert the list to a tagList - this is necessary for the list of items
            # to display properly.
            do.call(tagList, plot_output_list)
        })
        
        output$ti_traject_express <- renderUI({
            
            shinyFeedback::feedbackWarning("selected_genes_for_feature_plot_tab3",
                                           is.null(input$selected_genes_for_feature_plot_tab3),
                                           "Please, select one or more genes")
            req(input$selected_genes_for_feature_plot_tab3)
            
            feat_length <- length(req(input$selected_genes_for_feature_plot_tab3))
            
            plot_output_list <- lapply(1:feat_length, function(i) {
                plotname <- paste("ti_traject_expr", i, sep="")
                plotOutput(plotname, height = 300)
            })
            
            # Convert the list to a tagList - this is necessary for the list of items
            # to display properly.
            do.call(tagList, plot_output_list)
        })
        
        output$ti_graph_express <- renderUI({
            
            shinyFeedback::feedbackWarning("selected_genes_for_feature_plot_tab3",
                                           is.null(input$selected_genes_for_feature_plot_tab3),
                                           "Please, select one or more genes")
            req(input$selected_genes_for_feature_plot_tab3)
            
            feat_length <- length(req(input$selected_genes_for_feature_plot_tab3))
            
            plot_output_list <- lapply(1:feat_length, function(i) {
                plotname <- paste("ti_graph_expr", i, sep="")
                plotOutput(plotname, height = 300)
            })
            
            # Convert the list to a tagList - this is necessary for the list of items
            # to display properly.
            do.call(tagList, plot_output_list)
        })
        
        features_to_plot <- req(input$selected_genes_for_feature_plot_tab3)
        
        for ( i in 1:length( features_to_plot ) ) {
            
            model <- req( model() )
            
            local({
                
                my_i <- i
                
                plotname <- paste("ti_order_exp", my_i, sep="")
                output[[plotname]] <- renderPlot({
                    
                    plot_dimred(model,
                                feature_oi = features_to_plot[my_i],
                                expression_source = req( dataset_inf())) +
                        ggtitle("Cell grouping")
                    
                })
                
                plotname <- paste("ti_traject_expr", my_i, sep="")
                output[[plotname]] <- renderPlot({
                    
                    plot_dendro(model,
                                feature_oi = features_to_plot[my_i],
                                expression_source = req( dataset_inf())) +
                        ggtitle("Trajectory")#,
                    
                })
                
                plotname <- paste("ti_graph_expr", my_i, sep="")
                output[[plotname]] <- renderPlot({
                    
                    plot_graph(model,
                               feature_oi = features_to_plot[my_i],
                               expression_source = req( dataset_inf()) ) +
                        ggtitle("Trajectory represented as a graph")#,
                    
                })
                
            })
            
        }
        
    })
    
    
    output$select_genes_add_plot_to_down_tab3_ui <- renderUI({
        shinyFeedback::feedbackWarning("selected_genes_for_feature_plot_tab3",
                                       is.null(input$selected_genes_for_feature_plot_tab3),
                                       "Please, select one or more genes")
        req(input$selected_genes_for_feature_plot_tab3)
        
        filt_features <- req(input$selected_genes_for_feature_plot_tab3)
        
        div(class = "option-group",
            shinyWidgets::pickerInput(
                inputId = "select_genes_add_plot_to_down_tab3",
                label = "Select the genes that you want to download",
                choices = sort(as.character(filt_features)),
                multiple = TRUE,
                options = list(`actions-box` = TRUE)
            ))
        
    })
    
    observeEvent(input$start_down_add_plots_tab3, {
        
        shinyFeedback::feedbackWarning("select_genes_add_plot_to_down_tab3",
                                       is.null(input$select_genes_add_plot_to_down_tab3),
                                       "Please, select one or more genes")
        genes <- req(input$select_genes_add_plot_to_down_tab3)
        
        withProgress(message = "Please wait, preparing the data for download.",
                     value = 0.5, {
                         
                         
                         if (file.exists("./images") == F ) {
                             dir.create('./images')
                         }
                         
                         # Creates new folder to keep the plots organized. Adding date and time helps to avoid overwriting the data
                         path_new <- paste0("./images/trajectories_plots", format(Sys.time(),'_%Y-%m-%d__%H%M%S'))
                         
                         dir.create(path_new)
                         dir.create( paste0(path_new,"/Dimension_reduction") )
                         dir.create( paste0(path_new,"/Dendrogram") )
                         dir.create( paste0(path_new, "/Graph") )
                         
                         model <- req( model() )
                         
                         for( i in 1:length(genes) ){
                             
                             # Saves the feature plots
                             
                             p <- plot_dimred(model,
                                              feature_oi = genes[i],
                                              expression_source = req( dataset_inf()) ) +
                                 ggtitle("Cell grouping")
                             
                             file <- paste0(path_new, "/Dimension_reduction/", genes[i], ".", input$add_p_tab3_feat_format)
                             
                             ggplot2::ggsave(file,
                                             p,
                                             height=input$add_p_tab3_feat_height,
                                             width=input$add_p_tab3_feat_width,
                                             units="cm",
                                             dpi=as.numeric(input$add_p_tab3_feat_res))
                             
                             # Saves the violin plots
                             
                             file <- paste0(path_new, "/Dendrogram/", genes[i], ".", input$add_p_tab3_violin_format)
                             
                             p <- plot_dendro(model,
                                              feature_oi = genes[i],
                                              expression_source = req( dataset_inf()) )+
                                 ggtitle("Trajectory")
                             
                             ggplot2::ggsave(file,
                                             p,
                                             height=input$add_p_tab3_violin_height,
                                             width=input$add_p_tab3_violin_width,
                                             units="cm",
                                             dpi=as.numeric(input$add_p_tab3_violin_res))
                             
                             # Saves the dot plots
                             
                             file <- paste0(path_new, "/Graph/", genes[i], ".", input$add_p_tab3_dot_format)
                             
                             p <- plot_graph(model,
                                             feature_oi = genes[i],
                                             expression_source = req( dataset_inf())) +
                                 ggtitle("Trajectory represented as a graph")#,
                             
                             ggplot2::ggsave(file,
                                             p,
                                             height=input$add_p_tab3_dot_height,
                                             width=input$add_p_tab3_dot_width,
                                             units="cm",
                                             dpi=as.numeric(input$add_p_tab3_dot_res))
                             
                             
                         }
                         
                     })
        
        showNotification("All plots were downloaded!",
                         duration = 15,
                         id = "")
        
    })
    
    #################################################################
    ### Expression plots - Using dynverse most relevant features ####
    #################################################################
    
    observeEvent(c(input$run_heatmap_tab3_dynverse, input$dynverse_def_imp_genes), {
        
        # This calculate the height of the plot based on the n of genes
        heatmap_n_genes_tab3_dynv <- reactive({
            
            # Test what input to use for the heatmap
            if (input$ti_expre_opt == 0) {
                
                features <- req( filt_features_tab3() )
                features <- unique(features$GeneID)
                
            } else if (input$ti_expre_opt == 1) {
                
                features <- req( dynverse_genes_list_filt() )
                features <- features$feature_id
                
            }
            
            as.numeric(length(features))
        })
        heatmap_Height_tab3_dynv <- reactive( 250 + ( 20 * req( heatmap_n_genes_tab3_dynv() ) ) )
        
        # Generates the heatmap plot
        output$heat_map_tab3_dynv <- renderPlot({
            
            showNotification("Generating the heatmap",
                             id = "tab3_m8",
                             duration = NULL)
            
            # Test what input to use for the heatmap
            if (input$ti_expre_opt == 0) {
                
                features <- req( filt_features_tab3() )
                features <- unique(features$GeneID)
                
            } else if (input$ti_expre_opt == 1) {
                
                features <- req( dynverse_genes_list_filt() )
                features <- features$feature_id
                
            }
            
            features <- as.character(features)
            
            ### Drawing heat map ###
            heatmap <- plot_heatmap(
                req(  model() ),
                expression_source = req( dataset_inf() ),
                grouping = req( sc_meta_cluster() ),
                features_oi = features
            )
            
            on.exit(removeNotification(id = "tab3_m8"), add = TRUE)
            
            heatmap
            
        }, height = heatmap_Height_tab3_dynv() )
        
        output$heat_map_ui_tab3_dynv <- renderUI({
            
            plotOutput( "heat_map_tab3_dynv", height = heatmap_Height_tab3_dynv() )
            
        })
        
        output$marker_to_feature_plot_tab3_dynv = renderUI({
            
            # Test what input to use for the heatmap
            if (input$ti_expre_opt == 0) {
                
                features <- filt_features_tab3()
                features <- features$GeneID
                
            } else if (input$ti_expre_opt == 1) {
                
                features <- dynverse_genes_list_filt()
                features <- features$feature_id
                
            } else if (input$ti_expre_opt == 2) {
                
                features <- tradseq_genes_list()
                
            }
            
            div(class = "option-group",
                shinyWidgets::pickerInput(
                    inputId = "selected_genes_for_feature_plot_tab3_dynv",
                    label = "Select the genes to feature plots",
                    choices = sort(as.character(features)),
                    multiple = TRUE,
                    options = list(`actions-box` = TRUE)
                ))
            
        })
        
        output$p11_down <- downloadHandler(
            
            filename = function() {
                paste("Heatmap_int.dataset_most_important_features", ".", input$p11_format, sep = "")
            },
            content = function(file) {
                
                withProgress(message = "Please wait, preparing the data for download.",
                             value = 0.5, {
                                 
                                 shinyFeedback::feedbackWarning("p11_height", is.na(input$p11_height), "Required value")
                                 shinyFeedback::feedbackWarning("p11_width", is.na(input$p11_width), "Required value")
                                 
                                 height <- as.numeric( req( input$p11_height) )
                                 width <- as.numeric( req( input$p11_width) )
                                 res <- as.numeric(input$p11_res)
                                 
                                 # Test what input to use for the heatmap
                                 if (input$ti_expre_opt == 0) {
                                     
                                     features <- filt_features_tab3()
                                     features <- unique(features$GeneID)
                                     
                                 } else if (input$ti_expre_opt == 1) {
                                     
                                     features <- dynverse_genes_list_filt()
                                     features <- features$feature_id
                                     
                                 }
                                 
                                 features <- as.character(features)
                                 
                                 ### Drawing heat map ###
                                 p <- plot_heatmap(
                                     model(),
                                     expression_source = dataset_inf(),
                                     grouping = sc_meta_cluster(),
                                     features_oi = features
                                 )
                                 
                                 ggplot2::ggsave(file,
                                                 p,
                                                 height=height,
                                                 width=width,
                                                 units="cm",
                                                 dpi=res)
                                 
                             })
            }
            
        )
        
    })
    
    observeEvent(input$run_feature_plot_tab3_dynv, {
        
        showNotification("Generating the expression plots",
                         id = "tab3_m8",
                         duration = 30)
        
        output$ti_order_express_dynv <- renderUI({
            
            shinyFeedback::feedbackWarning("selected_genes_for_feature_plot_tab3_dynv",
                                           is.null(input$selected_genes_for_feature_plot_tab3_dynv),
                                           "Please, select one or more genes")
            req(input$selected_genes_for_feature_plot_tab3_dynv)
            
            feat_length <- length(input$selected_genes_for_feature_plot_tab3_dynv)
            
            plot_output_list <- lapply(1:feat_length, function(i) {
                plotname <- paste("ti_order_exp_dynv", i, sep="")
                plotOutput(plotname, height = 300)
            })
            
            # Convert the list to a tagList - this is necessary for the list of items
            # to display properly.
            do.call(tagList, plot_output_list)
        })
        
        output$ti_traject_express_dynv <- renderUI({
            
            shinyFeedback::feedbackWarning("selected_genes_for_feature_plot_tab3_dynv",
                                           is.null(input$selected_genes_for_feature_plot_tab3_dynv),
                                           "Please, select one or more genes")
            req(input$selected_genes_for_feature_plot_tab3_dynv)
            
            feat_length <- length(input$selected_genes_for_feature_plot_tab3_dynv)
            
            plot_output_list <- lapply(1:feat_length, function(i) {
                plotname <- paste("ti_traject_expr_dynv", i, sep="")
                plotOutput(plotname, height = 300)
            })
            
            # Convert the list to a tagList - this is necessary for the list of items
            # to display properly.
            do.call(tagList, plot_output_list)
        })
        
        output$ti_graph_express_dynv <- renderUI({
            
            shinyFeedback::feedbackWarning("selected_genes_for_feature_plot_tab3_dynv",
                                           is.null(input$selected_genes_for_feature_plot_tab3_dynv),
                                           "Please, select one or more genes")
            req(input$selected_genes_for_feature_plot_tab3_dynv)
            
            feat_length <- length(input$selected_genes_for_feature_plot_tab3_dynv)
            
            plot_output_list <- lapply(1:feat_length, function(i) {
                plotname <- paste("ti_graph_expr_dynv", i, sep="")
                plotOutput(plotname, height = 300)
            })
            
            # Convert the list to a tagList - this is necessary for the list of items
            # to display properly.
            do.call(tagList, plot_output_list)
        })
        
        features_to_plot <- input$selected_genes_for_feature_plot_tab3_dynv
        
        for (i in 1:length( features_to_plot ) ) {
            
            model <- model()
            
            local({
                
                my_i <- i
                
                plotname <- paste("ti_order_exp_dynv", my_i, sep="")
                output[[plotname]] <- renderPlot({
                    
                    plot_dimred(model,
                                feature_oi = features_to_plot[my_i],
                                expression_source = dataset_inf()) +
                        ggtitle("Cell grouping")
                    
                })
                
                plotname <- paste("ti_traject_expr_dynv", my_i, sep="")
                output[[plotname]] <- renderPlot({
                    
                    plot_dendro(model,
                                feature_oi = features_to_plot[my_i],
                                expression_source = dataset_inf()) +
                        ggtitle("Trajectory")#,
                    
                })
                
                plotname <- paste("ti_graph_expr_dynv", my_i, sep="")
                output[[plotname]] <- renderPlot({
                    
                    plot_graph(model,
                               feature_oi = features_to_plot[my_i],
                               expression_source = dataset_inf()) +
                        ggtitle("Trajectory represented as a graph")#,
                    
                })
                
            })
            
        }
        
    })
    
    output$select_genes_add_plot_to_down_tab3_ui_dynv <- renderUI({
        
        shinyFeedback::feedbackWarning("selected_genes_for_feature_plot_tab3_dynv",
                                       is.null(input$selected_genes_for_feature_plot_tab3_dynv),
                                       "Please, select one or more genes")
        req(input$selected_genes_for_feature_plot_tab3_dynv)
        
        filt_features <- input$selected_genes_for_feature_plot_tab3_dynv
        
        div(class = "option-group",
            shinyWidgets::pickerInput(
                inputId = "select_genes_add_plot_to_down_tab3_dynv",
                label = "Select the genes that you want to download",
                choices = sort(as.character(filt_features)),
                multiple = TRUE,
                options = list(`actions-box` = TRUE)
            ))
        
    })
    
    observeEvent(input$start_down_add_plots_tab3_dynv, {
        
        shinyFeedback::feedbackWarning("select_genes_add_plot_to_down_tab3_dynv",
                                       is.null(input$select_genes_add_plot_to_down_tab3_dynv),
                                       "Please, select one or more genes")
        req(input$select_genes_add_plot_to_down_tab3_dynv)
        
        withProgress(message = "Please wait, preparing the data for download.",
                     value = 0.5, {
                         
                         if (file.exists("./images") == F ) {
                             dir.create('./images')
                         }
                         
                         # Creates new folder to keep the plots organized. Adding date and time helps to avoid overwriting the data
                         path_new <- paste0("./images/trajectories_plots", format(Sys.time(),'_%Y-%m-%d__%H%M%S'))
                         
                         dir.create(path_new)
                         dir.create( paste0(path_new,"/Dimension_reduction") )
                         dir.create( paste0(path_new,"/Dendrogram") )
                         dir.create( paste0(path_new, "/Graph") )
                         
                         model <- model()
                         
                         
                         
                         for( i in 1:length(genes) ){
                             
                             # Saves the feature plots
                             
                             p <- plot_dimred(model,
                                              feature_oi = genes[i],
                                              expression_source = dataset_inf()) +
                                 ggtitle("Cell grouping")
                             
                             file <- paste0(path_new, "/Dimension_reduction/", genes[i], ".", input$add_p_tab3_feat_format_dynv)
                             
                             ggplot2::ggsave(file,
                                             p,
                                             height=input$add_p_tab3_feat_height_dynv,
                                             width=input$add_p_tab3_feat_width_dynv,
                                             units="cm",
                                             dpi=as.numeric(input$add_p_tab3_feat_res_dynv))
                             
                             # Saves the violin plots
                             
                             file <- paste0(path_new, "/Dendrogram/", genes[i], ".", input$add_p_tab3_violin_format_dynv)
                             
                             p <- plot_dendro(model,
                                              feature_oi = genes[i],
                                              expression_source = dataset_inf()) +
                                 ggtitle("Trajectory")
                             
                             ggplot2::ggsave(file,
                                             p,
                                             height=input$add_p_tab3_violin_height_dynv,
                                             width=input$add_p_tab3_violin_width_dynv,
                                             units="cm",
                                             dpi=as.numeric(input$add_p_tab3_violin_res_dynv))
                             
                             # Saves the dot plots
                             
                             file <- paste0(path_new, "/Graph/", genes[i], ".", input$add_p_tab3_dot_format_dynv)
                             
                             p <- plot_graph(model,
                                             feature_oi = genes[i],
                                             expression_source = dataset_inf()) +
                                 ggtitle("Trajectory represented as a graph")#,
                             
                             ggplot2::ggsave(file,
                                             p,
                                             height=input$add_p_tab3_dot_height_dynv,
                                             width=input$add_p_tab3_dot_width_dynv,
                                             units="cm",
                                             dpi=as.numeric(input$add_p_tab3_dot_res_dynv
                                                            
                                             ))
                             
                         }
                         
                     })
        
        showNotification("All plots were downloaded!",
                         duration = 15,
                         id = "")
        
    })
    
    output$download_dynverse_genes <- downloadHandler(
        
        filename = function() {
            paste("List_of_most_important_genes_within_trajectory", ".csv", sep = "")
        },
        content = function(file) {
            
            withProgress(message = "Please wait, preparing the data for download.",
                         value = 0.5, {
                             
                             write.csv(dynverse_genes_list(), file, row.names = FALSE)
                             
                         })
            
        }
        
    )
    
    output$download_dynverse_genes_filt <- downloadHandler(
        
        filename = function() {
            paste("_Filtered_list_of_most_important_genes_within_trajectory", ".csv", sep = "")
        },
        content = function(file) {
            
            withProgress(message = "Please wait, preparing the data for download.",
                         value = 0.5, {
                             
                             write.csv(dynverse_genes_list_filt(), file, row.names = FALSE)
                             
                         })
        }
    )
    
    ############################
    ### Annotation - BioMart ###
    ############################
    
    ##############################
    ### Load gene ids from csv ###
    ##############################
    genes_list <- reactive({
        req(input$genescsv) # Wait for file
        as.list(read.csv(input$genescsv$datapath, stringsAsFactors = FALSE)[,1])
    })
    
    #########################################
    ### Create a database selection input ###
    #########################################
    output$dbselection <- renderUI({
        shinyWidgets::pickerInput("dbselection",
                                  "Select BioMart db:",
                                  choices  = c('', 'PHYTOZOME', 'plants_mart', 'ENSEMBL_MART_ENSEMBL', 'ENSEMBL_MART_MOUSE',
                                               'ENSEMBL_MART_SNP', 'ENSEMBL_MART_FUNCGEN'),
                                  multiple = FALSE,
                                  #selected = "plants_mart",
                                  options  = pickerOptions(actionsBox = TRUE))
    })
    
    #########################################
    ### Create the selected mart variable ###
    #########################################
    selectedMart <- reactive({
        req(input$dbselection != '') # wait for a database to be selected
        
        # check if database is phytozome or ensembl
        # they have different hosts
        if (as.character(req(input$dbselection) == 'PHYTOZOME')) {
            phytozome_mart
        } else if (as.character(req(input$dbselection) == 'plants_mart')) {
            #useMart(as.character(req(input$dbselection)), host = "plants.ensembl.org", ensemblRedirect = FALSE)
            useMart(as.character(req(input$dbselection)), host = "plants.ensembl.org")
        } else {
            #useMart(as.character(req(input$dbselection)), host = "www.ensembl.org", ensemblRedirect = FALSE)
            useMart(as.character(req(input$dbselection)), host = "www.ensembl.org")
        }
    })
    
    #######################################
    ### Show Biomart available datasets ###
    ### (for the selected database)     ###
    #######################################
    # output$biomartdbs <- renderPrint({
    #     if (length(input$dbselection) != 0) {
    #         print(
    #             head(as.data.frame(
    #                 listDatasets(selectedMart())
    #             ), 5)
    #         )
    #     } else {
    #         print(
    #             "Please select a BioMart database!"
    #         )
    #     }
    # })
    
    ######################################
    ### Get available biomaRt datasets ###
    ######################################
    biomartdatasets <- reactive({
        req(input$dbselection != '') # wait for a database to be selected
        withProgress(
            message = "Getting available datasets in the database.",
            value   = 0.5,
            {
                retry(listDatasets(selectedMart())$dataset) # try connection 5 times
            }
        )
    })
    
    ########################################
    ### Create a dataset selection input ###
    ########################################
    output$datasetselection <- renderUI({
        req(input$dbselection, biomartdatasets()) # wait for db and dataset to be selected
        shinyWidgets::pickerInput("datasetselection",
                                  "Select BioMart dataset:",
                                  choices  = c('', biomartdatasets()),
                                  multiple = FALSE,
                                  #selected = "athaliana_eg_gene",
                                  options  = pickerOptions(actionsBox = TRUE))
    })
    
    #####################################
    ### Get available biomaRt filters ###
    #####################################
    datasetfilters <- reactive({
        if (input$datasetselection != '') {
            # This if statement makes the server wait for a dataset
            # to be chosen in order to search for available filters.
            withProgress(
                message = "Getting available filters in the selected dataset.",
                value   = 0.5,
                {
                    retry(
                        listFilters(useDataset(as.character(input$datasetselection), mart = selectedMart()))$name
                    ) # try connection 5 times
                }
            )
        } else {
            c('')
        }
    })
    
    #######################################
    ### Create a filter selection input ###
    #######################################
    output$filterselection <- renderUI({
        req(input$datasetselection != '', datasetfilters()) # wait for db and dataset to be selected
        
        # check if database is phytozome or ensembl
        # they have different filters
        if (as.character(input$dbselection == 'PHYTOZOME'))  {
            dt_filter <- "gene_name_filter"
        } else {
            dt_filter <- "ensembl_gene_id"
        }
        dt_filter
        
        shinyWidgets::pickerInput("filterselection",
                                  "Select a BioMart gene filter:",
                                  choices  = c('', datasetfilters()),
                                  multiple = FALSE,
                                  selected = dt_filter,
                                  options  = pickerOptions(actionsBox = TRUE))
    })
    
    # #########################################
    # ### Show Biomart available attributes ###
    # ### (for the selected db and dataset) ###
    # #########################################
    # output$biomartdatasets <- renderPrint({
    #     if (input$datasetselection != '') {
    #         print(
    #             head( biomaRt::listAttributes(
    #                 mart = useDataset(as.character(input$datasetselection), mart = selectedMart())), 10 )
    #         )
    #     } else {
    #         print(
    #             "Please select a BioMart dataset!"
    #         )
    #     }
    # })
    
    ##############################################
    ### Get available biomaRt attributes pages ###
    ##############################################
    attpages <- reactive({
        req(input$datasetselection != '') # wait for a dataset to be selected
        withProgress(
            message = "Getting available attribute pages in the selected dataset.",
            value   = 0.5,
            {
                retry(
                    attributePages(useDataset(as.character(input$datasetselection), mart = selectedMart()))
                ) # try connection 5 times
            }
        )
    })
    
    #################################################
    ### Create an ui to select available attPages ###
    #################################################
    output$attpageselection <- renderUI({
        req(input$datasetselection, attpages()) # wait for dataset to be selected
        shinyWidgets::pickerInput("attpageselection", "Select BioMart attributes Page:",
                                  choices  = attpages(),
                                  multiple = FALSE,
                                  options  = pickerOptions(actionsBox = TRUE))
    })
    
    ########################################
    ### Get available biomaRt attributes ###
    ########################################
    biomartatts <- reactive({
        req(input$datasetselection != '', input$attpageselection) # wait for a dataset to be selected
        withProgress(
            message = "Getting available attributes in the selected attribute page.",
            value   = 0.5,
            {
                retry(
                    as.list(
                        biomaRt::listAttributes(
                            mart = useDataset(as.character(input$datasetselection), mart = selectedMart()),
                            page = as.character(input$attpageselection)) %>%
                            unique()
                    )
                ) # try connection 5 times
            }
        )
    })
    
    ############################################
    ### Create an attributes selection input ###
    ############################################
    output$attselection <- renderUI({
        req(biomartatts()) # wait for list of available attributes
        
        # check if database is phytozome or ensembl
        # they have different attributes
        if (as.character(input$dbselection == 'PHYTOZOME'))  {
            att_list <- c("gene_name_filter","gene_chrom_start","gene_chrom_end","gene_description")
        } else {
            att_list <- c("ensembl_gene_id","start_position","end_position","description")
        }
        
        # render shinyWidgets::pickerInput
        shinyWidgets::pickerInput("attselection",
                                  "Select dataset attributes:",
                                  choices  = biomartatts(),  multiple = TRUE,
                                  selected = att_list,
                                  options  = pickerOptions(actionsBox = TRUE))
    })
    
    ########################################
    ### Select desired genes to annotate ###
    ########################################
    output$geneselection <- renderUI({
        req(genes_list()) # wait for genes
        shinyWidgets::pickerInput("geneselection",
                                  "Select genes to annotate:",
                                  choices  = genes_list(),
                                  multiple = TRUE,
                                  selected = genes_list(),
                                  options  = pickerOptions(actionsBox = TRUE))
    })
    
    ####################################
    ### Create button to trigger run ###
    ####################################
    output$triggerbutton <- renderUI({
        req(input$attselection, input$geneselection != '') # wait for complete selection
        actionButton("annotate", "Annotate selected genes!")
    })
    
    #############################
    ### BioMart output header ###
    #############################
    output$biomartresultsheader <- renderUI({
        req(input$attselection, input$geneselection != '') # wait for complete selection
        h4("BioMart results:")
    })
    
    ######################################################
    ### Create reactive data.frame to render datatable ###
    ######################################################
    annotation_biomart <- eventReactive(input$annotate, {
        
        withProgress(
            message = "Accessing biomart and annotating selected genes.",
            value   = 0.5,
            {
                retry(
                    biomaRt::getBM(
                        attributes = as.character(req(input$attselection)),
                        filters    = as.character(input$filterselection),
                        values     = as.character(req(input$geneselection)),
                        mart       = useDataset(as.character(input$datasetselection), mart = selectedMart()),
                        useCache   = FALSE
                    )
                ) # try connection 5 times
            }
        )
    })
    
    ############################
    ### BioMart output table ###
    ############################
    output$biomartresults <- renderDT(server = FALSE, {
        req(annotation_biomart()) # Input requirement
        
        # Load annotation
        anot_df <- annotation_biomart()
        
        # Render DT
        datatable(anot_df,
                  rownames = F,
                  selection = 'none',
                  filter = 'none',
                  extensions = 'Buttons',
                  options = list(pageLength = 5,
                                 lengthMenu = c(5, 10, 15, 20, 50),
                                 scrollX = "100%",
                                 autoWidth = FALSE,
                                 fixedHeader = TRUE,
                                 searching= FALSE,
                                 dom = 'Blfrtip',
                                 buttons = c('copy', 'csv', 'excel')))
    })
    
    #####################################################
    ### Grab gene universe for GO enrichment analysis ###
    #####################################################
    gene_universe <- reactive({
        req(input$universecsv) # wait for file
        as.list(read.csv(input$universecsv$datapath, stringsAsFactors = FALSE)[,1])
    })
    
    ####################################
    ### Create button to trigger run ###
    ####################################
    output$gobutton <- renderUI({
        req(gene_universe(), input$geneselection, input$filterselection) # wait for genes
        att_has_id(as.vector(genes_list()), as.vector(gene_universe()))
        actionButton("goenrich", "Run GO enrichment analysis!")
    })
    
    ##########################################################################
    ### Define right biomaRt attribute that correctly matches the gene ids ###
    ##########################################################################
    attselectionGO <- reactive({
        req(biomartatts()) # wait for attributes
        if (as.character(input$dbselection == 'PHYTOZOME'))  {
            att <- "gene_name1"
        } else {
            att <- "ensembl_gene_id"
        }
        att
    })
    
    #######################################
    ### Get GO information for universe ###
    #######################################
    goInfo <- eventReactive(input$goenrich, {
        req(gene_universe(), attselectionGO()) # wait for genes
        withProgress(
            message = "Retrieving GO information for the given gene universe.",
            value   = 0.5,
            {
                retry(
                    biomaRt::getBM(
                        attributes = c(as.character(attselectionGO()), "go_id"),
                        filters    = as.character(input$filterselection),
                        values     = as.character(req(gene_universe())),
                        mart       = useDataset(as.character(input$datasetselection), mart = selectedMart()),
                        useCache   = FALSE
                    )
                ) # try connection 5 times
            }
        )
        # %>%
        #     na_if("") %>%
        #     na.omit() %>%
        #     group_by_at(1) %>%
        #     summarise(go_id = paste(go_id, collapse=", "))
    })
    
    #########################################
    ### Get GO topGO package for universe ###
    #########################################
    topgo_universe <- reactive({
        
        req(goInfo()) # wait for data
        
        # Get GO information (together with gene ids/names)
        uni_GOs <-
            goInfo()
        
        # Write in temporary file
        write.table(
            uni_GOs,
            file = "/tmp/tmp_go_universe.txt",
            sep  = "\t",
            quote = FALSE,
            append = FALSE,
            row.names = FALSE,
            col.names = FALSE
        )
        
        # Load temporary file for topGO
        readMappings(file = "/tmp/tmp_go_universe.txt") # reads a named list
        
    })
    
    ##############################
    ### GO enrichment analysis ###
    ##############################
    enrichedGO <- reactive({
        req(topgo_universe(), input$geneselection != '') # Input requirement
        
        # check if identification is main ID
        if (as.character(input$filterselection) != as.character(attselectionGO())) {
            
            # target genes
            withProgress(
                message = "Retrieving GO information for the selected target genes.",
                value   = 0.5,
                {
                    tg_genes <- retry(
                        as.list(
                            biomaRt::getBM(
                                attributes = attselectionGO(),
                                filters    = as.character(input$filterselection),
                                values     = as.character(input$geneselection),
                                mart       = useDataset(as.character(input$datasetselection), mart = selectedMart()),
                                useCache   = FALSE
                            )[,1])
                    ) # try connection 5 times
                }
            )
            
            # universe genes
            uni_genes <- as.list(goInfo()[,1])
            
        } else {
            tg_genes  <- input$geneselection
            uni_genes <- gene_universe()
        }
        
        # Load target gene ids
        target_genes <- as.character(tg_genes)
        geneList <- factor(as.integer(uni_genes %in% target_genes))
        names(geneList) <- uni_genes
        
        # Execute topGO
        allRes <- data.frame()
        withProgress(
            message = "Performing GO enrichment analysis! This may take a while.",
            value   = 0.5,
            {
                for (GOcategory in c('MF', 'CC', 'BP')) {
                    res <- try(
                        runTopGO(GOcategory, geneList, topgo_universe()),
                        silent = TRUE
                    )
                    if (class(res) != 'try-error') {
                        allRes <- rbind(allRes, res)
                    }
                }
            }
        )
        
        # fix numeric
        allRes$KS <- as.numeric(allRes$KS)
        allRes$weightFisher <- as.numeric(allRes$weightFisher)
        
        # show table sorted
        allRes
    })
    
    ##################################
    ### Write enriched GO as table ###
    ##################################
    output$goresults <- renderDT(server = FALSE, {
        
        # Render DT
        datatable(enrichedGO() %>% arrange(!!sym(as.character(input$goplottest))),
                  rownames = F,
                  selection = 'none',
                  filter = 'none',
                  extensions = 'Buttons',
                  options = list(pageLength = 5,
                                 lengthMenu = c(5, 10, 15, 20, 50),
                                 scrollX = "100%",
                                 fixedHeader = TRUE,
                                 autoWidth = FALSE,
                                 searching= FALSE,
                                 dom = 'Blfrtip',
                                 buttons = c('copy', 'csv', 'excel')))
    })
    
    #####################################
    ### Generate plot of enriched GOs ###
    #####################################
    goplot <- reactiveValues(built = NULL)
    enrichedGOplot <- reactive({
        req(enrichedGO()) # require GO enrichment to finish
        
        # Get the number of top GOs user wants
        if (isTruthy(input$gontop)) {
            ntop <- as.numeric(input$gontop)
        } else {
            ntop <- 5
        }
        
        # get data
        ggdata <- enrichedGO()
        
        # fix numeric
        ggdata$KS <- as.numeric(ggdata$KS)
        ggdata$weightFisher <- as.numeric(ggdata$weightFisher)
        
        # Filter categories
        cat_MF <- ggdata %>% filter(Macro == "Molecular Function")
        cat_BP <- ggdata %>% filter(Macro == "Biological Process")
        cat_CC <- ggdata %>% filter(Macro == "Cellular Component")
        
        # order data by KS
        cat_MF <- cat_MF %>% arrange(!!sym(as.character(input$goplottest)))
        cat_BP <- cat_BP %>% arrange(!!sym(as.character(input$goplottest)))
        cat_CC <- cat_CC %>% arrange(!!sym(as.character(input$goplottest)))
        
        # filter ntop
        cat_MF <- cat_MF[1:ntop,]
        cat_BP <- cat_BP[1:ntop,]
        cat_CC <- cat_CC[1:ntop,]
        
        # concat
        ggdata <- na.omit(rbind(cat_MF, cat_BP, cat_CC))
        
        # fix order of terms
        ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term))
        
        # generate plot
        gg1 <- plotTopGO(ggdata, as.character(input$goplottest))
        
        # change state
        goplot$built <- TRUE
        
        # call plot
        gg1
    })
    
    #############################################
    ### Create buttons for plot customization ###
    #############################################
    # Hide customization buttons until plot is built
    output$plotbuilt <- reactive({
        return(!is.null(goplot$built))
    })
    outputOptions(output, 'plotbuilt', suspendWhenHidden= FALSE)
    
    ##############################
    ### Show enriched GOs plot ###
    ##############################
    output$goplot <- renderPlot({
        req(enrichedGOplot())
        enrichedGOplot()
    })
    
    ################################
    ### Download button for plot ###
    ################################
    output$goplot_down <- downloadHandler(
        filename = function() { paste('enrichedGOplot', as.character(input$goplotdevice), sep='.') },
        content = function(file) {
            
            withProgress(message = "Please wait, preparing the data for download.",
                         value = 0.5, {
                             
                             # Get the number of top GOs user wants
                             if (isTruthy(input$gontop)) {
                                 ntop <- as.numeric(input$gontop)
                             } else {
                                 ntop <- 5
                             }
                             
                             # get data
                             ggdata <- enrichedGO()
                             
                             # fix numeric
                             ggdata$KS <- as.numeric(ggdata$KS)
                             ggdata$weightFisher <- as.numeric(ggdata$weightFisher)
                             
                             # Filter categories
                             cat_MF <- ggdata %>% filter(Macro == "Molecular Function")
                             cat_BP <- ggdata %>% filter(Macro == "Biological Process")
                             cat_CC <- ggdata %>% filter(Macro == "Cellular Component")
                             
                             # order data by KS
                             cat_MF <- cat_MF %>% arrange(!!sym(as.character(input$goplottest)))
                             cat_BP <- cat_BP %>% arrange(!!sym(as.character(input$goplottest)))
                             cat_CC <- cat_CC %>% arrange(!!sym(as.character(input$goplottest)))
                             
                             # filter ntop
                             cat_MF <- cat_MF[1:ntop,]
                             cat_BP <- cat_BP[1:ntop,]
                             cat_CC <- cat_CC[1:ntop,]
                             
                             # concat
                             ggdata <- na.omit(rbind(cat_MF, cat_BP, cat_CC))
                             
                             # fix order of terms
                             ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term))
                             
                             # generate plot
                             gg1 <- plotTopGO(ggdata, as.character(input$goplottest))
                             
                             ggplot2::ggsave(
                                 file,
                                 plot   = gg1,
                                 device = as.character(input$goplotdevice),
                                 width  = as.numeric(input$goplotwidth),
                                 height = as.numeric(input$goplotheight),
                                 dpi    = as.numeric(input$goplotdpi),
                                 units  = 'cm'
                             )
                             
                         })
        }
        
    )
    
    
    #####################
    ### Stacked plots ###
    #####################
    
    stacked_violin_Server("stacked1")
    
    
    # output$load_rds_ui <- stacked_violin_Server("load_rds_ui")
    # output$list_of_genes_ui <- stacked_violin_Server("list_of_genes_ui")
    # output$list_of_genes <- stacked_violin_Server("stacked_violing")
    # output$stacked_violing_ui <- stacked_violin_Server("stacked_violing_ui")
    #output$load_rds_ui_stacked <- load_rds_Server("load_rds_ui_stacked")
    
    # output$TEST <- renderReactable({
    #     my_reactable(filt_features())
    # })
    
}
