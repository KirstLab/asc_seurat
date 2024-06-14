# SolusCell web app
# Version 1.0
set.seed(1407)
options(shiny.sanitize.errors = T)
options(max.print=100)

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
suppressMessages( require(shinycssloaders) )
suppressMessages( require(DT) )
suppressMessages( require(dplyr) )
suppressMessages( require(hdf5r) )
suppressMessages( require(scales) )
suppressMessages( require(utils) )
suppressMessages( require(vroom) )
suppressMessages( require(showtext) )

# Bioconductor
suppressMessages( require(ComplexHeatmap) )
suppressMessages( require(SingleCellExperiment) )
suppressMessages( require(slingshot) )
suppressMessages( require(multtest) )
suppressMessages( require(glmGamPoi) ) # Bioconductor

## New packages 
suppressMessages( require(scCustomize) )
suppressMessages( require(cowplot) )

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

font_add_google("Montserrat", "montserrat")
showtext_auto()

min_cells = 3
min_features = 1
most_var_method = "vst"

function(input, output, session) {
    
    if (dir.exists('/app/user_work')) {
        setwd('/app/user_work')
    }
    
    #####################################
    ######   Tab 1 - Clustering    ######
    #####################################
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
    
    output$select_sample_tab1_rds_ui <- renderUI({
        
        rds_list <- list.files('./RDS_files/', pattern = "*.rds")
        
        div(class = "option-group",
            shinyWidgets::pickerInput(
                inputId = "select_sample_tab1_rds",
                label = "Select the file containing the data",
                choices = sort(rds_list),
                multiple = FALSE,
                options = list(`actions-box` = TRUE)
            ))
        
    })
    
    is_load_10X_rds_clicked <- reactiveVal(FALSE)
    observeEvent( input$load_10X_rds, {
        is_load_10X_rds_clicked(TRUE)
    })
    
    is_run_clustering_clicked <- reactiveVal(FALSE)
    observeEvent( input$run_clustering, {
        is_run_clustering_clicked(TRUE)
    })
    
    single_cell_data_reac <- eventReactive( input$load_10X, {
        
        shinyFeedback::feedbackWarning("min_cells", is.na(min_cells), "Required value")
        shinyFeedback::feedbackWarning("min_features", is.na(min_features), "Required value")
        shinyFeedback::feedbackWarning("mito_regex", !shiny::isTruthy( isolate( input$mito_regex) ), "Required value")
        
        req(min_cells)
        req(min_features)
        req( isolate( input$mito_regex) )
        
        showNotification("Loading the data",
                         duration = NULL,
                         id = "p1")
        
        sing_cell_data.data <- Seurat::Read10X(data.dir = req(  isolate( input$sample_folder_tab1)) )
        
        # Initialize the Seurat object with the raw (non-normalized data).
        sing_cell_data <- Seurat::CreateSeuratObject(counts = sing_cell_data.data,
                                                     project =  isolate( input$proj_name ),
                                                     min.cells = min_cells,
                                                     min.features = min_features)
        
        # Calculate the % of mithocondrial contamination
        sing_cell_data <- Seurat::PercentageFeatureSet(sing_cell_data,
                                                       pattern =  isolate( input$mito_regex),
                                                       col.name = "percent.mt")
        
        return(sing_cell_data)
        
    })
    
    observeEvent( input$load_10X, {
        
        output$target_genes_mitho <- renderPrint({
            
            sc_object <- req(single_cell_data_reac())
            mito_regex <- req( isolate( input$mito_regex) )
            assay <- DefaultAssay(object = sc_object)
            
            mitho_genes <- grep(pattern = mito_regex,
                                x = rownames(x = sc_object[[assay]]), 
                                value = TRUE)
            
            if ( length(mitho_genes) <= 0) {
                
                print("No gene was target by the identifier!")
                
            }else {
                
                mitho_genes
                
            }
            
        })
        
    })
    
    output$VlnPlot <- renderPlot({
        
        data_set <- req( single_cell_data_reac() )
        on.exit(removeNotification(id = "p1"), add = TRUE)
        
        myVlnPlot(sc_data = data_set,
                  show_box_median = input$Vln_box)
        
    }, height = 500)
    
    observeEvent( input$run_vinplot, {
        
        data_sc <- req( single_cell_data_reac() )
        
        test_cond <- if( !is.na( isolate( input$max_count) ) && !is.na( isolate( input$min_count) ) ) {
            
            shinyFeedback::feedbackWarning("max_count",
                                           isolate( input$min_count ) >  isolate( input$max_count),
                                           "No cells will be selected by appling this parameters!")
            
            shinyFeedback::feedbackWarning("min_count",
                                           isolate( input$min_count ) >  isolate( input$max_count) ,
                                           "No cells will be selected by appling this parameters!")
            
            validate( need(  isolate( input$min_count ) <  isolate( input$max_count ),
                             "Error: No cells will be selected by appling this parameters!"))
            
        }
        
        if ( !is.na(  isolate( input$min_count) ) ) {
            
            data_sc <- base::subset(data_sc,
                                    subset = nFeature_RNA >  isolate( input$min_count) )
            
        }
        
        if ( !is.na( isolate( input$max_count) ) ) {
            
            shinyFeedback::feedbackWarning("max_count",
                                           isolate( input$max_count ) <= 0,
                                           "No cells will be selected by appling this parameters!")
            
            validate( need(  isolate( input$max_count ) > 0, "Error: No cells will be selected by appling this parameters!"))
            
            data_sc <- base::subset(data_sc,
                                    subset =  nFeature_RNA <  isolate( input$max_count) )
            
        }
        
        if ( !is.na(  isolate( input$max_mito_perc) ) ) {
            
            shinyFeedback::feedbackWarning("max_mito_perc",
                                           isolate( input$max_mito_perc ) <= 0,
                                           "No cells will be selected by appling this parameters!")
            
            validate(need(  isolate( input$max_mito_perc ) > 0, "Error: No cells will be selected by appling this parameters!"))
            
            data_sc <- base::subset(data_sc,
                                    subset =  percent.mt <  isolate( input$max_mito_perc) )
            
        }
        
        single_cell_data_filt <- data_sc
        
        output$VlnPlot_filt <- renderPlot({
            
            myVlnPlot(sc_data = single_cell_data_filt,
                      show_box_median = input$Vln_box)
            
        }, height = 500)
        
    })
    
    output$p1_down <- downloadHandler(
        
        filename = function() {
            paste("Violin_plot", ".",  isolate( input$p1_format ), sep = "")
        },
        content = function(file) {
            
            withProgress(message = "Please wait, preparing the data for download.",
                         value = 0.5, {
                             
                             data_sc <- req( single_cell_data_reac() )
                             
                             test_cond <- if( !is.na(  isolate( input$max_count) ) && !is.na(  isolate( input$min_count) ) ) {
                                 
                                 shinyFeedback::feedbackWarning("max_count",
                                                                isolate( input$min_count ) >  isolate( input$max_count) ,
                                                                "No cells will be selected by appling this parameters!")
                                 
                                 shinyFeedback::feedbackWarning("min_count",
                                                                isolate( input$min_count ) > isolate( input$max_count ),
                                                                "No cells will be selected by appling this parameters!")
                                 
                                 validate( need( isolate( input$min_count ) < isolate( input$max_count ),
                                                 "Error: No cells will be selected by appling this parameters!"))
                                 
                             }
                             
                             if ( !is.na( isolate( input$min_count) ) ) {
                                 
                                 data_sc <- base::subset(data_sc,
                                                         subset = nFeature_RNA > isolate( input$min_count) )
                                 
                             }
                             
                             if ( !is.na( isolate( input$max_count) ) ) {
                                 
                                 shinyFeedback::feedbackWarning("max_count",
                                                                isolate( input$max_count ) <= 0,
                                                                "No cells will be selected by appling this parameters!")
                                 
                                 validate( need( isolate( input$max_count ) > 0, "Error: No cells will be selected by appling this parameters!"))
                                 
                                 data_sc <- base::subset(data_sc,
                                                         subset =  nFeature_RNA < isolate( input$max_count) )
                                 
                             }
                             
                             if ( !is.na( isolate( input$max_mito_perc) ) ) {
                                 
                                 shinyFeedback::feedbackWarning("max_mito_perc",
                                                                isolate( input$max_mito_perc ) <= 0,
                                                                "No cells will be selected by appling this parameters!")
                                 
                                 validate( need( isolate( input$max_mito_perc ) > 0, "Error: No cells will be selected by appling this parameters!"))
                                 
                                 data_sc <- base::subset(data_sc,
                                                         subset =  percent.mt < isolate( input$max_mito_perc) )
                                 
                             }
                             
                             height <- as.numeric( req( isolate(  input$p1_height) ) )
                             width <- as.numeric( req( isolate(  input$p1_width) ) )
                             res <- as.numeric( req(  isolate( input$p1_res) ) )
                             
                             p <-  myVlnPlot(sc_data = single_cell_data_filt,
                                             show_box_median = input$Vln_box)
                             
                             ggplot2::ggsave(file,
                                             p,
                                             height=height,
                                             width=width,
                                             units="cm",
                                             dpi=res,
                                             bg = "#FFFFFF")
                             
                         })
            
        }
    )
    
    single_cell_data_pca <- eventReactive( c( input$run_pca, input$rerun_after_filtering ), {
        
        sc_data <- req( single_cell_data_reac() )
        
        if ( isolate( input$filter_clusters) == 1) { to_filter <- req( to_filter() )}
        
        req( isolate( input$min_count) )
        req( isolate( input$max_count) )
        req( isolate( input$max_mito_perc) )
        req( isolate( input$filter_clusters) )
        req( isolate( input$filter_clusters_opt) )
        req( isolate( most_var_method) )
        
        if ( isolate( input$normaliz_method) == 0) { #LogNormalize
            
            shinyFeedback::feedbackWarning("scale_factor", is.na(  isolate( input$scale_factor) ), "Required value")
            shinyFeedback::feedbackWarning("n_of_var_genes", is.na( isolate( input$n_of_var_genes) ), "Required value")
            
            req( isolate( input$scale_factor) )
            req( isolate( input$n_of_var_genes) )
            
            sc_data2 <- lognorm_function("Test",
                                         to_filter = to_filter,
                                         sc_data = sc_data,
                                         min_count =  isolate( input$min_count),
                                         max_count =  isolate( input$max_count),
                                         max_mito_perc =  isolate( input$max_mito_perc),
                                         scale_factor =  isolate( input$scale_factor),
                                         filter_clusters =  isolate( input$filter_clusters),
                                         filter_clusters_opt =  isolate( input$filter_clusters_opt),
                                         most_var_method =  isolate( most_var_method),
                                         n_of_var_genes =  isolate( input$n_of_var_genes) )
            
        } else if ( isolate( input$normaliz_method == 1) ) { #"SCTransform"
            
            sc_data2 <- SCTransform_function("Test",
                                             to_filter = to_filter,
                                             sc_data = sc_data,
                                             min_count =  isolate( input$min_count ),
                                             max_count =  isolate( input$max_count ),
                                             max_mito_perc =  isolate( input$max_mito_perc ),
                                             filter_clusters =  isolate( input$filter_clusters ),
                                             filter_clusters_opt =  isolate( input$filter_clusters_opt) )
            
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
    
    observeEvent( c( input$rerun_after_filtering), {
        
        output$n_of_PCAs <- renderPlot({
            
            data_sc <- req( single_cell_data_pca() )
            
            Seurat::ElbowPlot(data_sc, ndims = 50, reduction = "pca")
            
        })
        
    })
    
    # ## This remove the plots and reset the PCA value when the user select or exclude clusters. The goal is to make the user evaluate the parameters before executing
    is_remove_plots_clicked <- reactiveVal(FALSE)
    
    observeEvent( input$rerun_after_filtering, {
        is_remove_plots_clicked(TRUE)
    })
    
    observeEvent( input$rerun_after_filtering, {
        
        if ( is_remove_plots_clicked() ) {
            updateNumericInput(session, "n_of_PCs", value = NA) ## <<<<
        }
        
    })
    
    output$p2_down <- downloadHandler(
        
        filename = function() {
            paste("Elbow_Plot", ".",  isolate( input$p2_format ), sep = "")
        },
        content = function(file) {
            
            withProgress(message = "Please wait, preparing the data for download.",
                         value = 0.5, {
                             
                             shinyFeedback::feedbackWarning("p2_height", is.na( isolate( input$p2_height) ), "Required value")
                             shinyFeedback::feedbackWarning("p2_width", is.na( isolate( input$p2_width) ), "Required value")
                             
                             height <- as.numeric( req(  isolate( input$p2_height) ) )
                             width <- as.numeric( req(  isolate( input$p2_width) ) )
                             res <- as.numeric( isolate( input$p2_res) )
                             
                             p <- Seurat::ElbowPlot(req( single_cell_data_pca()), ndims = 50, reduction = "pca")
                             
                             ggplot2::ggsave(file,
                                             p,
                                             height=height,
                                             width=width,
                                             units="cm",
                                             dpi=res,
                                             bg = "#FFFFFF")
                             
                         })
        }
        
    )
    
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
                                       is.null(  isolate( input$cluster_list) ),
                                       "Required value")
        req( isolate( input$cluster_list) )
        
        if (  isolate( input$filter_clusters_opt) == 0 ) { # "select"
            
            to_filter <- base::subset(req( single_cell_data_reso_umap() ),
                                      idents = as.numeric( isolate( input$cluster_list)) )
            
            to_filter_ch <- to_filter@meta.data
            to_filter_ch <- base::rownames(to_filter_ch)
            
            to_filter_ch
            
        } else if (  isolate( input$filter_clusters_opt) == 1 ) { # exclude
            
            to_filter <- base::subset(req( single_cell_data_reso_umap() ),
                                      idents = as.numeric( isolate( input$cluster_list) ),
                                      invert = TRUE)
            
            to_filter_ch <- to_filter@meta.data
            to_filter_ch <- base::rownames(to_filter_ch)
            
            to_filter_ch
        }
        
    })
    
    # set reactive values for data, initialize, and plot type
    rv = reactiveValues(type1 = NULL, type2 = NULL, PCs = NULL)
    
    # set plot type at the time "go" was clicked and data
    observeEvent( input$run_clustering, {
        
        shinyFeedback::feedbackWarning("n_of_PCs",
                                       is.na( isolate( input$n_of_PCs) ),
                                       "Required value")
        shinyFeedback::feedbackWarning("resolution_clust", 
                                       is.na( isolate( input$resolution_clust) ),
                                       "Required value")
        
        rv$type2 = "umap"
        rv$type3 = "tsne"
        rv$PCs =  isolate( req(input$n_of_PCs ) )
    })
    
    single_cell_data_reso_umap <- eventReactive(  list( input$run_clustering, input$load_10X_rds ), {
        
        if (  isolate( input$sample_tab1_options) == 1 & is_load_10X_rds_clicked() ) { # Load file
            
            showNotification("Loading the data",
                             id = "m9",
                             duration = NULL)
            
            sc_data <- readRDS( paste0("./RDS_files/", req( isolate( input$select_sample_tab1_rds) ) ) )
            
            validate(need("seurat_clusters" %in% colnames(sc_data@meta.data),
                          "Clustering was not detected in the rds file", ""))
            
            on.exit(removeNotification(id = "m9"), add = TRUE)
            
        } else if ( isolate( input$sample_tab1_options) == 0 ) { # new analysis
            
            if (is.na( isolate( input$n_of_PCs) ) ||  isolate( input$n_of_PCs) != rv$PCs || is.null(rv$PCs)) {
                NULL
            }  else {
                
                req( isolate( input$n_of_PCs) )
                req( isolate( input$resolution_clust) )
                
                showNotification("Running the clustering step",
                                 duration = NULL,
                                 id = "m6")
                
                data_sc <- req( single_cell_data_pca() )
                
                data_sc <- Seurat::FindNeighbors(data_sc, dims = 1:isolate( isolate( input$n_of_PCs) ) )
                
                data_sc <- Seurat::FindClusters(data_sc, resolution =  isolate( input$resolution_clust) )
                
                sc_data <- Seurat::RunUMAP(data_sc, dims = 1:isolate(input$n_of_PCs))
                sc_data <- Seurat::RunTSNE(sc_data, dims = 1:isolate(input$n_of_PCs))
                
            }
        }
        sc_data
    })
    
    
    
    observeEvent( list( input$run_clustering, input$load_10X_rds ), {
        
        if ( is_load_10X_rds_clicked() | is_run_clustering_clicked() ) {
            
            output$tSNE <- renderPlot({
                isolate( input$run_clustering)
                
                if (is.na( isolate( input$n_of_PCs) ) ||  isolate( input$n_of_PCs ) != rv$PCs || is.null(rv$PCs)) {
                    NULL
                }  else {
                    
                    Seurat::DimPlot( req( single_cell_data_reso_umap() ), reduction = "tsne", label = T, pt.size = .1)
                    
                }
            })
            
            output$umap <- renderPlot({
                isolate( input$run_clustering )
                
                if (is.na( isolate( input$n_of_PCs) ) ||  isolate( input$n_of_PCs ) != rv$PCs || is.null(rv$PCs)) {
                    NULL
                }  else {
                    Seurat::DimPlot(req( single_cell_data_reso_umap()), reduction = "umap", label = T, pt.size = .1)
                }
            })
            
            output$cluster_size_plot <- renderPlot({
                
                sing_cell_data <- req( single_cell_data_reso_umap() )
                
                sc_meta <- as.data.frame(sing_cell_data[[]])
                sc_meta$cellcluster <- base::rownames(sc_meta)
                sc_meta <- sc_meta[, c("cellcluster", "seurat_clusters")]
                
                on.exit(removeNotification(id = "m6"), add = TRUE)
                
                table_n <- as.data.frame(base::table(sc_meta$seurat_clusters))
                table_n %>%
                    dplyr::rename(Cluster = "Var1",
                                  `N. of cells` = "Freq")
                
                table_n <- table_n %>%
                    dplyr::rename(Cluster = "Var1",
                                  `N. of cells` = "Freq")

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
                
            })
            
            # output$cluster_size <-  renderTable({
            #     
            #     sing_cell_data <- req( single_cell_data_reso_umap() )
            #     
            #     sc_meta <- as.data.frame(sing_cell_data[[]])
            #     sc_meta$cellcluster <- base::rownames(sc_meta)
            #     sc_meta <- sc_meta[, c("cellcluster", "seurat_clusters")]
            #     
            #     on.exit(removeNotification(id = "m6"), add = TRUE)
            #     
            #     table_n <- as.data.frame(base::table(sc_meta$seurat_clusters))
            #     table_n %>%
            #         dplyr::rename(Cluster = "Var1",
            #                       `N. of cells` = "Freq")
            #     
            # })
        }
    })
    
    output$p3_down <- downloadHandler(
        
        filename = function() {
            paste("clustering_plot_",   isolate( input$p3_down_opt ), ".",  isolate( input$p3_format ), sep = "")
        },
        content = function(file) {
            
            withProgress(message = "Please wait, preparing the data for download.",
                         value = 0.5, {
                             
                             shinyFeedback::feedbackWarning("p3_height", is.na( isolate( input$p3_height) ), "Required value")
                             shinyFeedback::feedbackWarning("p3_width", is.na( isolate( input$p3_width) ), "Required value")
                             
                             height <- as.numeric( req(  isolate( input$p3_height) ) )
                             width <- as.numeric( req(  isolate( input$p3_width) ) )
                             res <- as.numeric( isolate( input$p3_res) )
                             
                             if (  isolate( input$p3_down_opt) == "UMAP" ) {
                                 
                                 p <- Seurat::DimPlot( req( single_cell_data_reso_umap() ), reduction = "umap", label = T, pt.size = .1)
                                 
                             } else if (  isolate( input$p3_down_opt) == "t-SNE") {
                                 
                                 p <- Seurat::DimPlot(req( single_cell_data_reso_umap() ), reduction = "tsne", label = T, pt.size = .1)
                                 
                             }
                             
                             ggplot2::ggsave(file,
                                             p,
                                             height=height,
                                             width=width,
                                             units="cm",
                                             dpi=res,
                                             bg = "#FFFFFF")
                             
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
        clusters <- clusters[ !clusters ==  isolate( input$find_markers_clust_ID1_tab1 )]
        
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
                                       is.na( isolate( input$find_markers_tab1_return_thresh) ),
                                       "Required value")
        shinyFeedback::feedbackWarning("find_markers_tab1_logfc_threshold",
                                       is.na( isolate( input$find_markers_tab1_logfc_threshold) ),
                                       "Required value")
        shinyFeedback::feedbackWarning("find_markers_tab1_min_pct",
                                       is.na( isolate( input$find_markers_tab1_min_pct) ),
                                       "Required value")
        
        req( isolate( input$find_markers_tab1_return_thresh) )
        req( isolate( input$find_markers_tab1_logfc_threshold) )
        req( isolate( input$find_markers_tab1_min_pct) )
        
        sc_data <- req( single_cell_data_reso_umap() )
        
        if ( isolate( input$find_markers_tab1_opt) == 2) { #distinguishing a cluster from other
            
            ident_1 =  isolate( input$find_markers_clust_ID1_tab1 )
            ident_2 =  isolate( input$find_markers_clust_ID2_tab1 )
            
        } else if ( isolate( input$find_markers_tab1_opt) == 1 ) {
            
            ident_1 =  isolate( input$find_markers_clust_id_tab1)
            
        }
        
        showNotification("Identifing markers or D.E. genes",
                         duration = NULL,
                         id = "tab1_n1")
        
        finding_markers("finding_markers_tab1",
                        sc_data,
                        assay_choice = "RNA",
                        find_markers_tab1_opt =  isolate( input$find_markers_tab1_opt),
                        find_markers_tab1_return.thresh =  isolate( input$find_markers_tab1_return_thresh),
                        find_markers_tab1_logfc.threshold =  isolate( input$find_markers_tab1_logfc_threshold),
                        find_markers_tab1_min.pct =  isolate( input$find_markers_tab1_min_pct),
                        find_markers_tab1_test.use =  isolate( input$find_markers_tab1_test.use),
                        ident.1 = ident_1,
                        ident.2 = ident_2,
                        find_markers_tab1_filt_pos =  isolate( input$find_markers_tab1_filt_pos)
        )
        
    })
    
    output$markers_tab1_react <- renderReactable({
        
        markers_tab1 <- req( markers_tab1() )
        
        if ( isolate( input$find_markers_tab1_opt) == 1) {
            
            markers_tab1$cluster <- isolate(input$find_markers_clust_id_tab1)
            markers_tab1 <- markers_tab1[, c( 1, ncol(markers_tab1), 2 : ( ncol(markers_tab1) - 1) ) ]
            
        } else if ( isolate( input$find_markers_tab1_opt ) == 2) {
            
            markers_tab1$cluster <- paste0( isolate(input$find_markers_clust_ID1_tab1),
                                            "_vs_" ,
                                            isolate(input$find_markers_clust_ID2_tab1))
            markers_tab1 <- markers_tab1[, c( 1, ncol(markers_tab1), 2 : ( ncol(markers_tab1) - 1) ) ]
            
        } else if ( isolate( input$find_markers_tab1_opt) == 0) {
            
            markers_tab1 <- markers_tab1[, c(1, ncol(markers_tab1), 2 : ( ncol(markers_tab1) - 1) ) ]
            
        }
        
        on.exit(removeNotification(id = "tab1_n1"), add = TRUE)
        
        my_reactable(markers_tab1)
        
    })
    
    output$download_markers_tab1 <- downloadHandler(
        
        filename = function() {
            
            if (  isolate( input$find_markers_tab1_opt) == 0 ) {
                
                paste("list_of_markers_all_clusters", ".csv", sep = "")
                
            } else if ( isolate( input$find_markers_tab1_opt) == 1) {
                
                paste("list_of_markers_", "cluster",  isolate( input$find_markers_clust_id_tab1), ".csv", sep = "")
                
            } else if ( isolate( input$find_markers_tab1_opt) == 2) {
                
                paste("list_of_markers_", "cluster",  isolate( input$find_markers_clust_ID1_tab1), "vs" ,  isolate( input$find_markers_clust_ID2_tab1), ".csv", sep = "")
                
            }
            
        },
        content = function(file) {
            
            withProgress(message = "Please wait, preparing the data for download.",
                         value = 0.5, {
                             
                             markers_tab1 <- req( markers_tab1() )
                             
                             if ( isolate( input$find_markers_tab1_opt) == 1) {
                                 
                                 markers_tab1$cluster <-  isolate( input$find_markers_clust_id_tab1)
                                 markers_tab1 <- markers_tab1[, c( 1, ncol(markers_tab1), 2 : ( ncol(markers_tab1) - 1) ) ]
                                 
                             } else if ( isolate( input$find_markers_tab1_opt ) == 2) {
                                 
                                 markers_tab1$cluster <- paste0( isolate( input$find_markers_clust_ID1_tab1), "_vs_" ,  isolate( input$find_markers_clust_ID2_tab1) )
                                 markers_tab1 <- markers_tab1[, c( 1, ncol(markers_tab1), 2 : ( ncol(markers_tab1) - 1) ) ]
                                 
                             } else if ( isolate( input$find_markers_tab1_opt) == 0) {
                                 
                                 markers_tab1 <- markers_tab1[, c(1, ncol(markers_tab1), 2 : ( ncol(markers_tab1) - 1) ) ]
                                 
                             }
                             
                             write.csv(markers_tab1, file, row.names = FALSE)
                             
                         })
        }
    )
    
    features <- reactive({
        
        req( isolate( input$markers_list) )
        ext <- tools::file_ext( isolate( input$markers_list$name) )
        markers_list_file <-  isolate( input$markers_list$datapath )
        
        if (  isolate( input$markers_list_header_opt ) == "" ) {
            
            shinyFeedback::feedbackWarning("markers_list_header_opt",
                                           TRUE,
                                           "Please, inform if the file has a header.")
            
        }
        
        req( isolate( input$markers_list_header_opt) )
        
        read_file_markers("tab1_readfile",
                          markers_list_file = markers_list_file,
                          feed_ID ="markers_list",
                          ext = ext,
                          header_opt =  isolate( input$markers_list_header_opt) )
    })
    
    observeEvent(input$load_markers, {
        
        output$marker_group_selec = renderUI({
            
            pickerInput_markers_group("features_group", genes = req( features() ) )
            
        })
        
        output$marker_genes_selec = renderUI({
            
            features_group <- req( isolate( input$features_group) )
            id_choice <- req( isolate( input$genes_ids) )
            
            pickerInput_markers_genes("selected_genes",
                                      genes = req( features() ),
                                      features_group = features_group,
                                      id_choice = id_choice)
        })
        
    })
    
    filt_features <- eventReactive(input$run_heatmap, {
        
        features_f <- req( features() )
        features_group <- req( isolate( input$features_group) )
        
        features_f <- dplyr::filter(features_f,
                                    Group %in% features_group)
        
        if ( isolate( input$filter_genes_q) == 0) {
            
            shinyFeedback::feedbackWarning("selected_genes",
                                           is.null( isolate( input$selected_genes) ),
                                           "Please, select one or more genes")
            req( isolate( input$selected_genes) )
            
            if ( isolate( input$genes_ids) == "ID") {
                
                features_f <- features_f[features_f$GeneID %in%  isolate( input$selected_genes), ]
                
            } else if ( isolate( input$genes_ids) == "name") {
                
                features_f <- features_f[features_f$Name %in%  isolate( input$selected_genes ), ]
                
            }
            
        }
        
        features_f
        
    })
    
    sc_data_av_react <- eventReactive( input$run_heatmap, {
        
        sc_data_av <- Seurat::AggregateExpression( req( single_cell_data_reso_umap() ),
                                                   assays = "RNA"#,
                                                   #layer = input$slot_selection_heatmap
        )
        
        sc_data_av <- sc_data_av[[1]]
        
        on.exit(removeNotification(id = "m7"), add = TRUE)
        sc_data_av
        
    })
    
    observeEvent( input$run_heatmap, {
        
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
                sc_data_av_feat <- t(sc_data_av_feat)
                
                base::rownames(sc_data_av_feat) <- features_selec[1, 1]
                
            } else {
                
                sc_data_av_feat <- sc_data_av[base::rownames(sc_data_av) %in% features_selec[, 1], ]
                
            }
            sc_data_av_feat
        })
        
        output$heat_map <- renderPlot({
            
            req(heat_map_prep())
            heat_map_prep <- heat_map_prep()
            heat_map_prep <- as(heat_map_prep, "sparseMatrix") 
            
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
                                         genes_ids =  isolate( input$genes_ids) )
            
        })
        
        output$p4_down <- downloadHandler(
            
            filename = function() {
                paste("Heatmap", ".",  isolate( input$p4_format ), sep = "")
            },
            content = function(file) {
                
                withProgress(message = "Please wait, preparing the data for download.",
                             value = 0.5, {
                                 
                                 shinyFeedback::feedbackWarning("p4_height", is.na( isolate( input$p4_height) ), "Required value")
                                 shinyFeedback::feedbackWarning("p4_width", is.na( isolate( input$p4_width) ), "Required value")
                                 
                                 height <- as.numeric( req(  isolate( input$p4_height) ) )
                                 width <- as.numeric( req(  isolate( input$p4_width) ) )
                                 res <- as.numeric( isolate( input$p4_res) )
                                 
                                 heat_map_prep <- req( heat_map_prep() )
                                 heat_map_prep <- as(heat_map_prep, "sparseMatrix")
                                 
                                 p <- ComplexHeatmap::Heatmap(heat_map_prep, border = TRUE,
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
                                                 dpi=res,
                                                 bg = "#FFFFFF")
                                 
                             })
            }
            
        )
        
    })
    
    features_selec <- eventReactive( input$run_feature_plot, {
        
        shinyFeedback::feedbackWarning("selected_genes_for_feature_plot",
                                       is.null( isolate( input$selected_genes_for_feature_plot) ),
                                       "Please, select one or more genes.")
        req( isolate( input$selected_genes_for_feature_plot) )
        
        features <- req( filt_features() )
        
        if ( isolate( input$genes_ids) == "ID") {
            
            features_f <- dplyr::filter(features,
                                        GeneID %in% input$selected_genes_for_feature_plot)
            
        }
        else if (input$genes_ids == "name") {
            
            features_f <- dplyr::filter(features,
                                        Name %in% input$selected_genes_for_feature_plot)
            
        }
        
        as.character(unique(features_f$GeneID))
        
    })
    
    observeEvent( input$run_feature_plot, {
        
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
                                                          layer =  isolate( input$slot_selection_feature_plot),
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
                                                          layer =  isolate( input$slot_selection_feature_plot),
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
                                    layer =  isolate( input$slot_selection_feature_plot) )
                    
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
            
            filt_features <- req( isolate( input$selected_genes_for_feature_plot) )
            
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
    
    observeEvent( input$start_down_add_plots_tab1, {
        
        on.exit(removeNotification(id = "m8"), add = TRUE)
        
        shinyFeedback::feedbackWarning("select_genes_add_plot_to_down",
                                       is.null( isolate( input$select_genes_add_plot_to_down) ),
                                       "Please, select one or more genes.")
        
        genes <- req( isolate( input$select_genes_add_plot_to_down) )
        
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
                                                                       layer =  isolate( input$slot_selection_feature_plot),
                                                                       reduction = "umap") +
                                                       scale_colour_gradient2(limits=c(minimal, maximal),
                                                                              midpoint = maximal / 2,
                                                                              low = "gray80",
                                                                              mid = "gold",
                                                                              high = "red"))
                             
                             file <- paste0(path_new, "/feature_plots/", genes[i], ".",  isolate( input$add_p_tab1_feat_format) )
                             
                             ggplot2::ggsave(file,
                                             p,
                                             height= isolate( input$add_p_tab1_feat_height),
                                             width= isolate( input$add_p_tab1_feat_width),
                                             units="cm",
                                             dpi=as.numeric( isolate( input$add_p_tab1_feat_res) ),
                                             bg = "#FFFFFF")
                             
                             file <- paste0(path_new, "/violin_plots/", genes[i], ".",  isolate( input$add_p_tab1_violin_format) )
                             
                             p <- Seurat::VlnPlot(sc_data,
                                                  features = genes[i],
                                                  layer =  isolate( input$slot_selection_feature_plot) )
                             
                             ggplot2::ggsave(file,
                                             p,
                                             height= isolate( input$add_p_tab1_violin_height ),
                                             width= isolate( input$add_p_tab1_violin_width),
                                             units="cm",
                                             dpi=as.numeric( isolate( input$add_p_tab1_violin_res) ),
                                             bg = "#FFFFFF")
                             
                             file <- paste0(path_new, "/dot_plots/", genes[i], ".",  isolate( input$add_p_tab1_dot_format) )
                             
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
                                             height= isolate( input$add_p_tab1_dot_height),
                                             width= isolate( input$add_p_tab1_dot_width),
                                             units="cm",
                                             dpi=as.numeric( isolate( input$add_p_tab1_dot_res) ),
                                             bg = "#FFFFFF")
                         }
                         
                     })
        
        showNotification("All plots were downloaded!",
                         duration = 15,
                         id = "")
        
    })
    
    
    #####################################
    ### Tab 2 - Integration pipeline ####
    #####################################
    # 
    # output$load_integrated_ui <- renderUI({
    #     
    #     rds_list <- list.files('./RDS_files/', pattern = "*.rds")
    #     
    #     div(class = "option-group",
    #         shinyWidgets::pickerInput(
    #             inputId = "load_integrated",
    #             label = "Select the file containing the integrated data",
    #             choices = sort(rds_list),
    #             multiple = FALSE,
    #             options = list(`actions-box` = TRUE)
    #         ))
    #     
    # })
    # 
    # samples_list_integration <- reactive({
    #     
    #     req(input$samples_list_integration)
    #     
    #     ext <- tools::file_ext(input$samples_list_integration$name)
    #     config_input_file <- input$samples_list_integration$datapath
    #     
    #     test <- ext %in% c(
    #         'text/csv',
    #         'text/comma-separated-values',
    #         'text/tab-separated-values',
    #         'csv',
    #         'tsv')
    #     
    #     shinyFeedback::feedbackDanger("samples_list_integration",
    #                                   test == F,
    #                                   "Format is not supported! Upload a CSV or TSV file.")
    #     
    #     config_input <- switch(ext,
    #                            csv = vroom::vroom(config_input_file, delim = ","),
    #                            tsv = vroom::vroom(config_input_file, delim = "\t"),
    #                            validate("Invalid file; Please upload a .csv or .tsv file")
    #     )
    #     config_input <- as.data.frame(config_input)
    #     
    #     if ( ncol(config_input) >= 6 ) {
    #         shinyFeedback::feedbackSuccess("samples_list_integration",
    #                                        T, "")
    #     } else {
    #         shinyFeedback::feedbackDanger("samples_list_integration",
    #                                       T,
    #                                       "The config. file must have at least six columns!")
    #         
    #         validate( need( ncol(config_input) >= 6,
    #                         "The config. file must have at least six columns!",
    #                         "samples_list_integration") )
    #     }
    #     
    #     config_input
    #     
    # })
    # 
    # output$select_sample_tab2 <- renderUI({
    #     
    #     dir_list <- req( samples_list_integration() )
    #     dir_list <- unique(dir_list[, 1])
    #     
    #     div(class = "option-group",
    #         shinyWidgets::pickerInput(
    #             inputId = "sample_folder_tab2",
    #             label = "Select the samples to use",
    #             choices = sort(as.character(dir_list)),
    #             multiple = TRUE,
    #             options = list(`actions-box` = TRUE),
    #         ))
    #     
    # })
    # 
    # single_cell_data_reac_tab2 <- eventReactive( c(input$load_rds_file, input$load_rds_file3), {
    #     
    #     if ( input$integration_options == 0 ) { # new analysis
    #         
    #         shinyFeedback::feedbackWarning("sample_folder_tab2", is.null(input$sample_folder_tab2), "Required value")
    #         shinyFeedback::feedbackWarning("int_regex_mito", !shiny::isTruthy(input$int_regex_mito), "Required value")
    #         
    #         selected_samples <- req(input$sample_folder_tab2)
    #         config_csv <- req( samples_list_integration() )
    #         int_regex_mito <- req(input$int_regex_mito)
    #         
    #         # Select only the samples that the user specified.
    #         config_csv <- config_csv[config_csv[, 1] %in% selected_samples, ]
    #         
    #         project_name <- input$int_project_name
    #         
    #         if ( input$normaliz_method_tab2 == 0 ) { # LogNormalize
    #             
    #             shinyFeedback::feedbackWarning("sample_folder_tab2", is.null(input$sample_folder_tab2), "Required value")
    #             shinyFeedback::feedbackWarning("scale_factor_tab2", is.na(input$scale_factor_tab2), "Required value")
    #             shinyFeedback::feedbackWarning("n_of_var_genes_integration", is.na(input$n_of_var_genes_integration), "Required value")
    #             shinyFeedback::feedbackWarning("n_of_PCs_integration", is.na(input$n_of_PCs_integration), "Required value")
    #             
    #             scale_factor_tab2 <- req(input$scale_factor_tab2)
    #             n_of_var_genes_integration <- req(input$n_of_var_genes_integration)
    #             n_of_PCs_integration <- req(input$n_of_PCs_integration)
    #             
    #             showNotification("Integrating the data. Please wait, it can take a few minutes.",
    #                              id = "m12",
    #                              duration = NULL)
    #             
    #             sc_data <- integration_lognorm(config_csv = config_csv,
    #                                            project_name = project_name,
    #                                            int_regex_mito = int_regex_mito,
    #                                            scale_factor_tab2 = scale_factor_tab2,
    #                                            most_var_method_integration = most_var_method_tab2,
    #                                            n_of_var_genes_integration = n_of_var_genes_integration,
    #                                            n_of_PCs_integration = n_of_PCs_integration)
    #             
    #             on.exit(removeNotification(id = "m12"), add = TRUE)
    #             
    #             sc_data
    #             
    #         } else if ( input$normaliz_method_tab2 == 1 ) { # SCTransform
    #             
    #             showNotification("Integrating the data. Please wait, it can take a few minutes.",
    #                              id = "m12",
    #                              duration = NULL)
    #             
    #             sc_data <- integration_sctransform(config_csv = config_csv,
    #                                                project_name = project_name,
    #                                                int_regex_mito = int_regex_mito)
    #             
    #             on.exit(removeNotification(id = "m12"), add = TRUE)
    #             
    #             sc_data
    #             
    #         }
    #         
    #     } else if ( input$integration_options == 2 ) {
    #         
    #         ##  the user to tell what normalization was used, since we should not scale the data when using SCTransform.
    #         shinyFeedback::feedbackWarning("load_rds_int_normalization",
    #                                        input$load_rds_int_normalization == "",
    #                                        #T,
    #                                        "Please select an option")
    #         validate(need(input$load_rds_int_normalization != "",
    #                       message = "",
    #                       label = "load_rds_int_normalization"))
    #         
    #         showNotification("Loading the integrated data",
    #                          id = "m9",
    #                          duration = NULL)
    #         
    #         sc_data <- readRDS( paste0("./RDS_files/", req(input$load_integrated)) )
    #         
    #         shinyFeedback::feedbackWarning("load_integrated",
    #                                        "seurat_clusters" %in% colnames(sc_data@meta.data),
    #                                        "Clustering detected. Please use the option to load a clustered dataset")
    #         
    #         `%!in%` = Negate(`%in%`)     
    #         validate(need("seurat_clusters" %!in% colnames(sc_data@meta.data),
    #                       "Clustering detected. Please use the option to load a clustered dataset", ""))
    #         
    #         validate(need("treat" %in% colnames(sc_data@meta.data),
    #                       "Integrated data not detected. Does the rds file contain an integrated dataset?", ""))
    #         
    #         on.exit(removeNotification(id = "m9"), add = TRUE)
    #         
    #         sc_data
    #     }
    #     
    #     sc_data
    #     
    # })
    # 
    # observeEvent( input$load_rds_file, {
    #     
    #     output$target_genes_mitho_tab2 <- renderPrint({
    #         
    #         sc_object <- req(single_cell_data_reac_tab2())
    #         mito_regex <- req(input$int_regex_mito)
    #         
    #         assay <- DefaultAssay(object = sc_object)
    #         
    #         mito_genes <- grep(pattern = mito_regex,
    #                            x = rownames(x = sc_object[[assay]]), 
    #                            value = TRUE)
    #         
    #         if ( length(mito_genes) <= 0) {
    #             
    #             print("No gene target by the common identifier!")
    #             
    #         }else {
    #             
    #             mito_genes
    #             
    #         }
    #         
    #     })
    #     
    # })
    # 
    # output$download_int_data <- downloadHandler(
    #     
    #     filename = function() {
    #         paste("Integrated_datasets_without_clutering", ".rds", sep = "")
    #     },
    #     content = function(file) {
    #         
    #         withProgress(message = "Please wait, preparing the data for download.",
    #                      value = 0.5, {
    #                          
    #                          saveRDS( req( single_cell_data_reac_tab2() ), file)
    #                          
    #                      })
    #         
    #     }
    #     
    # )
    # 
    # output$VlnPlot_tab2 <- renderPlot({
    #     
    #     data_set <- req( single_cell_data_reac_tab2() )
    #     
    #     Seurat::VlnPlot(data_set,
    #                     features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    #                     ncol = 3,
    #                     split.plot = F)
    #     
    # })
    # 
    # single_cell_data_filt_tab2 <- reactive({
    #     
    #     data_sc <- req( single_cell_data_reac_tab2() )
    #     
    #     test_cond <- if( !is.na(input$max_count_tab2) && !is.na(input$min_count_tab2) ) {
    #         
    #         shinyFeedback::feedbackWarning("max_count_tab2",
    #                                        input$min_count_tab2 > input$max_count_tab2,
    #                                        "No cells will be selected by appling this parameters!")
    #         
    #         shinyFeedback::feedbackWarning("min_count_tab2",
    #                                        input$min_count_tab2 > input$max_count_tab2,
    #                                        "No cells will be selected by appling this parameters!")
    #         validate(need(input$min_count_tab2 < input$max_count_tab2, "Error: No cells will be selected by appling this parameters!"))
    #         
    #     }
    #     
    #     if ( !is.na(input$min_count_tab2) ) {
    #         
    #         data_sc <- base::subset(data_sc,
    #                                 subset = nFeature_RNA > input$min_count_tab2)
    #         
    #     }
    #     
    #     if ( !is.na(input$max_count_tab2) ) {
    #         
    #         shinyFeedback::feedbackWarning("max_count_tab2",
    #                                        input$max_count_tab2 <= 0,
    #                                        "No cells will be selected by appling this parameters!")
    #         
    #         validate(need(input$max_count_tab2 > 0, "Error: No cells will be selected by appling this parameters!"))
    #         
    #         data_sc <- base::subset(data_sc,
    #                                 subset =  nFeature_RNA < input$max_count_tab2)
    #         
    #     }
    #     
    #     if ( !is.na(input$max_mito_perc_tab2) ) {
    #         
    #         shinyFeedback::feedbackWarning("max_mito_perc_tab2", input$max_mito_perc_tab2 <= 0, "No cells will be selected by appling this parameters!")
    #         validate(need(input$max_mito_perc_tab2 > 0, "Error: No cells will be selected by appling this parameters!"))
    #         
    #         data_sc <- base::subset(data_sc,
    #                                 subset =  percent.mt < input$max_mito_perc_tab2)
    #         
    #     }
    #     
    #     data_sc
    #     
    # })
    # 
    # observeEvent(input$run_vinplot_tab2, {
    #     
    #     output$VlnPlot_filt_tab2 <- renderPlot({
    #         
    #         Seurat::VlnPlot(req( single_cell_data_filt_tab2() ),
    #                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    #                         ncol = 3,
    #                         split.plot = F)
    #         
    #     })
    #     
    #     output$p5_down <- downloadHandler(
    #         
    #         filename = function() {
    #             paste("VlnPlot", ".", input$p5_format, sep = "")
    #         },
    #         content = function(file) {
    #             
    #             withProgress(message = "Please wait, preparing the data for download.",
    #                          value = 0.5, {
    #                              
    #                              shinyFeedback::feedbackWarning("p5_height", is.na(input$p5_height), "Required value")
    #                              shinyFeedback::feedbackWarning("p5_width", is.na(input$p5_width), "Required value")
    #                              
    #                              height <- as.numeric( req( input$p5_height) )
    #                              width <- as.numeric( req( input$p5_width) )
    #                              res <- as.numeric(input$p5_res)
    #                              
    #                              p <- Seurat::VlnPlot(req( single_cell_data_filt_tab2() ),
    #                                                   features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    #                                                   ncol = 3,
    #                                                   split.plot = F)
    #                              
    #                              ggplot2::ggsave(file,
    #                                              p,
    #                                              height=height,
    #                                              width=width,
    #                                              units="cm",
    #                                              dpi=res,
    #                                              bg = "#FFFFFF")
    #                              
    #                          })
    #         }
    #         
    #     )
    #     
    # })
    # 
    # single_cell_data_scaled_tab2 <- eventReactive(
    #     c(input$run_pca_tab2_1,
    #       input$run_pca_tab2_2), {
    #           
    #           single_cell_data_filt_tab2 <- req( single_cell_data_filt_tab2() )
    #           
    #           # If loading the data and normalization is SCTransform, skip the scaling. Same if running new analysis and normalization is SCtransform
    #           
    #           if (input$integration_options != 0 && input$load_rds_int_normalization == 1) {
    #               
    #               ## Same if running new analysis and normalization is SCtransform
    #           } else if (input$integration_options == 0 && input$normaliz_method_tab2 == 1) {
    #               
    #           } else { # lognormalization
    #               
    #               showNotification("Scalling the data",
    #                                duration = NULL,
    #                                id = "tab2_m4")
    #               
    #               single_cell_data_filt_tab2 <- Seurat::ScaleData(single_cell_data_filt_tab2,
    #                                                               verbose = T)
    #               
    #               on.exit(removeNotification(id = "tab2_m4"), add = TRUE)
    #               
    #           }
    #           
    #           single_cell_data_filt_tab2
    #       })
    # 
    # single_cell_data_pca_tab2 <- eventReactive(c(input$run_pca_tab2_1,
    #                                              input$run_pca_tab2_2,
    #                                              input$rerun_after_filtering_tab2), {
    #                                                  
    #                                                  if (input$filter_clusters_tab2 == 1) {
    #                                                      
    #                                                      sc_data <- req( single_cell_data_scaled_tab2_filtered() )
    #                                                      
    #                                                  } else if ( input$filter_clusters_tab2 == 0 ) {
    #                                                      
    #                                                      sc_data <- req( single_cell_data_scaled_tab2() )
    #                                                      
    #                                                  }
    #                                                  
    #                                                  showNotification("Running PCA",
    #                                                                   duration = NULL,
    #                                                                   id = "m5")
    #                                                  
    #                                                  sc_data <-  RunPCA(sc_data,
    #                                                                     verbose = T)
    #                                                  
    #                                                  on.exit(removeNotification(id = "m5"), add = TRUE)
    #                                                  
    #                                                  sc_data
    #                                                  
    #                                              })
    # 
    # output$cluster_list_tab2_ui <- renderUI ({
    #     
    #     clusters <- req( clusters_single_cell_data_reso_umap_tab2() )
    #     
    #     shinyWidgets::pickerInput(
    #         inputId = "cluster_list_tab2",
    #         label = "Choose clusters to select or exclude",
    #         choices = sort(clusters),
    #         multiple = TRUE,
    #         options = list(`actions-box` = TRUE)
    #     )
    #     
    # })
    # 
    # to_filter_tab2 <- reactive({
    #     
    #     shinyFeedback::feedbackWarning("cluster_list_tab2", is.null(input$cluster_list_tab2), "Required value")
    #     req(input$cluster_list_tab2)
    #     
    #     sc_data <- req( single_cell_data_clustered() )
    #     
    #     if ( input$filter_clusters_opt_tab2 == 0 ) { #"select"
    #         
    #         to_filter <- base::subset(sc_data,
    #                                   idents = as.numeric(input$cluster_list_tab2)
    #         )
    #         
    #         to_filter_ch <- to_filter@meta.data
    #         to_filter_ch <- base::rownames(to_filter_ch)
    #         
    #         to_filter_ch
    #         
    #     } else if ( input$filter_clusters_opt_tab2 == 1 ) { #"exclude"
    #         
    #         to_filter <- base::subset(sc_data,
    #                                   idents = as.numeric(input$cluster_list_tab2),
    #                                   invert = TRUE)
    #         
    #         to_filter_ch <- to_filter@meta.data
    #         to_filter_ch <- base::rownames(to_filter_ch)
    #         
    #         to_filter_ch
    #     }
    #     
    # })
    # 
    # single_cell_data_scaled_tab2_filtered <- eventReactive(input$rerun_after_filtering_tab2, {
    #     
    #     sc_data <- req( single_cell_data_filt_tab2() )
    #     cells_to_filter <- req( to_filter_tab2() )
    #     
    #     if ( input$filter_clusters_tab2 == 1 ) {
    #         
    #         sc_data <- base::subset(sc_data,
    #                                 cells = cells_to_filter)
    #         
    #     }
    #     
    #     
    #     if (input$integration_options != 0 && input$load_rds_int_normalization == 1) {
    #         
    #     } else if (input$integration_options == 0 && input$normaliz_method_tab2 == 1) {
    #         
    #     } else { # lognormalization
    #         
    #         showNotification("Scalling the data",
    #                          duration = NULL,
    #                          id = "tab2_m4")
    #         
    #         sc_data <- Seurat::ScaleData(sc_data, verbose = T)
    #         
    #         on.exit(removeNotification(id = "tab2_m4"), add = TRUE)
    #         
    #         sc_data
    #     }
    #     
    #     sc_data
    # })
    # 
    # observeEvent(c(input$run_pca_tab2_1,
    #                input$run_pca_tab2_2,
    #                input$rerun_after_filtering_tab2), {
    #                    
    #                    output$n_of_PCAs_tab2 <- renderPlot({
    #                        
    #                        data_sc <- req( single_cell_data_pca_tab2() )
    #                        
    #                        Seurat::ElbowPlot(data_sc, ndims = 50, reduction = "pca")
    #                        
    #                    })
    #                    
    #                    output$p6_down <- downloadHandler(
    #                        
    #                        filename = function() {
    #                            paste("ElbowPlot", ".", input$p6_format, sep = "")
    #                        },
    #                        content = function(file) {
    #                            
    #                            withProgress(message = "Please wait, preparing the data for download.",
    #                                         value = 0.5, {
    #                                             
    #                                             shinyFeedback::feedbackWarning("p6_height", is.na(input$p6_height), "Required value")
    #                                             shinyFeedback::feedbackWarning("p6_width", is.na(input$p6_width), "Required value")
    #                                             
    #                                             height <- as.numeric( req( input$p6_height) )
    #                                             width <- as.numeric( req( input$p6_width) )
    #                                             res <- as.numeric(input$p6_res)
    #                                             
    #                                             data_sc <- req( single_cell_data_pca_tab2() )
    #                                             
    #                                             p <- Seurat::ElbowPlot(data_sc, ndims = 50, reduction = "pca")
    #                                             
    #                                             ggplot2::ggsave(file,
    #                                                             p,
    #                                                             height=height,
    #                                                             width=width,
    #                                                             units="cm",
    #                                                             dpi=res,
    #                                                             bg = "#FFFFFF")
    #                                             
    #                                         })
    #                        }
    #                        
    #                    )
    #                    
    #                })
    # 
    # single_cell_data_clustered <- eventReactive( c(input$run_clustering_tab2, input$load_rds_file2), {
    #     
    #     if ( input$integration_options == 1) { # Load file - Clustered
    #         
    #         showNotification("Loading the integrated data",
    #                          id = "m9",
    #                          duration = NULL)
    #         
    #         sc_data <- readRDS( paste0("./RDS_files/", req(input$load_integrated)) )
    #         
    #         validate(need("treat" %in% colnames(sc_data@meta.data),
    #                       "Integrated data not detected. Does the rds file contain an integrated dataset?", ""))
    #         
    #         validate(need("seurat_clusters" %in% colnames(sc_data@meta.data),
    #                       "Clustering was not detected. Please use the option for unclustered data", ""))
    #         
    #         on.exit(removeNotification(id = "m9"), add = TRUE)
    #         
    #         sc_data
    #         
    #     } else if (input$integration_options == 0 || input$integration_options == 2 ) { # new analysis or loading and integrated and unclustered data
    #         
    #         shinyFeedback::feedbackWarning("n_of_PCs_tab2", is.na(input$n_of_PCs_tab2), "Required value")
    #         shinyFeedback::feedbackWarning("resolution_clust_tab2", is.na(input$resolution_clust_tab2), "Required value")
    #         
    #         req(input$n_of_PCs_tab2)
    #         req(input$resolution_clust_tab2)
    #         
    #         showNotification("Running the clustering step",
    #                          duration = NULL,
    #                          id = "m6")
    #         
    #         sc_data <- req( single_cell_data_pca_tab2() )
    #         
    #         sc_data <- Seurat::RunUMAP(sc_data,
    #                                    reduction = "pca",
    #                                    dims = 1:input$n_of_PCs_tab2)
    #         
    #         sc_data <- Seurat::FindNeighbors(sc_data,
    #                                          reduction = "pca",
    #                                          dims = 1:input$n_of_PCs_tab2)
    #         
    #         sc_data <- Seurat::FindClusters(sc_data,
    #                                         resolution = input$resolution_clust_tab2)
    #         
    #     }
    #     
    #     sc_data
    #     
    # })
    # 
    # observeEvent(c(input$run_clustering_tab2, input$load_rds_file2), {
    #     
    #     output$umap_tab2 <- renderPlot({
    #         
    #         Seurat::DimPlot(req( single_cell_data_clustered() ), reduction = "umap", label = T, pt.size = .1)
    #         
    #     })
    #     
    #     output$umap_three_samples_comb <- renderPlot({
    #         
    #         Seurat::DimPlot(req( single_cell_data_clustered() ),
    #                         reduction = "umap",
    #                         label = F,
    #                         group.by = "treat",
    #                         pt.size = .1
    #         ) + ggtitle("Samples")
    #         
    #     })
    #     
    #     output$umap_three_samples <- renderPlot({
    #         
    #         Seurat::DimPlot(req( single_cell_data_clustered() ),
    #                         reduction = "umap",
    #                         label = T,
    #                         split.by = "treat", pt.size = .1
    #         )
    #     })
    #     
    #     output$p7_down <- downloadHandler(
    #         
    #         filename = function() {
    #             paste("Clustering_plot_", input$p7_down_opt, ".", input$p7_format, sep = "")
    #         },
    #         content = function(file) {
    #             
    #             withProgress(message = "Please wait, preparing the data for download.",
    #                          value = 0.5, {
    #                              
    #                              shinyFeedback::feedbackWarning("p7_height", is.na(input$p7_height), "Required value")
    #                              shinyFeedback::feedbackWarning("p7_width", is.na(input$p7_width), "Required value")
    #                              
    #                              height <- as.numeric( req( input$p7_height) )
    #                              width <- as.numeric( req( input$p7_width) )
    #                              res <- as.numeric(input$p7_res)
    #                              
    #                              if ( input$p7_down_opt == "UMAP" ) {
    #                                  
    #                                  p <- Seurat::DimPlot(req( single_cell_data_clustered() ),
    #                                                       reduction = "umap",
    #                                                       label = T,
    #                                                       pt.size = .1)
    #                                  
    #                              } else if ( input$p7_down_opt == "UMAP1" ) {
    #                                  
    #                                  p <- Seurat::DimPlot(req( single_cell_data_clustered() ),
    #                                                       reduction = "umap",
    #                                                       label = F,
    #                                                       group.by = "treat",
    #                                                       pt.size = .1
    #                                  ) + ggtitle("Samples")
    #                                  
    #                              } else if ( input$p7_down_opt == "UMAP2" ) {
    #                                  
    #                                  p <- Seurat::DimPlot(req( single_cell_data_clustered() ),
    #                                                       reduction = "umap",
    #                                                       label = T,
    #                                                       split.by = "treat", pt.size = .1
    #                                  )
    #                                  
    #                              }
    #                              
    #                              ggplot2::ggsave(file,
    #                                              p,
    #                                              height=height,
    #                                              width=width,
    #                                              units="cm",
    #                                              dpi=res,
    #                                              bg = "#FFFFFF")
    #                              
    #                          })
    #         }
    #         
    #     )
    #     
    #     output$cluster_size_tab2 <- renderPrint({
    #         
    #         sing_cell_data <- req( single_cell_data_clustered() )
    #         
    #         sc_meta <- as.data.frame(sing_cell_data[[]])
    #         sc_meta$cellcluster <- base::rownames(sc_meta)
    #         sc_meta <- sc_meta[, c("cellcluster", "seurat_clusters")]
    #         
    #         # returns the number of cells in each cluster
    #         cell_per_cluster <- sc_meta$seurat_clusters
    #         
    #         on.exit(removeNotification(id = "m6"), add = TRUE)
    #         
    #         table(sc_meta$seurat_clusters)
    #         
    #         
    #     })
    #     
    # } )
    # 
    # clusters_single_cell_data_reso_umap_tab2 <- reactive({
    #     
    #     sc_data <- req( single_cell_data_clustered() )
    #     as.numeric( unique( as.character(sc_data@meta.data$seurat_clusters ) ) )
    #     
    # })
    # 
    # treat_single_cell_data_reso_umap_tab2 <- reactive({
    #     
    #     sc_data <- req( single_cell_data_clustered() )
    #     unique( as.character(sc_data@meta.data$treat ) )
    #     
    # })
    # 
    # output$downloadRDS_tab2 <- downloadHandler(
    #     
    #     filename = function() {
    #         paste("Integrated_datasets_after_clutering", ".rds", sep = "")
    #     },
    #     content = function(file) {
    #         
    #         withProgress(message = "Please wait, preparing the data for download.",
    #                      value = 0.5, {
    #                          
    #                          saveRDS( req( single_cell_data_clustered() ), file)
    #                          
    #                      })
    #         
    #     }
    #     
    # )
    # 
    # output$find_markers_clust_id_tab2_ui <- renderUI ({
    #     
    #     clusters <- req( clusters_single_cell_data_reso_umap_tab2() )
    #     
    #     shinyWidgets::pickerInput(
    #         inputId = "find_markers_clust_id_tab2",
    #         label = "Select the cluster of interest",
    #         choices = sort(clusters),
    #         multiple = FALSE,
    #         options = list(`actions-box` = TRUE)
    #     )
    #     
    # })
    # 
    # output$find_markers_clust_ID1_tab2_ui <- renderUI ({
    #     
    #     clusters <- req( clusters_single_cell_data_reso_umap_tab2() )
    #     
    #     shinyWidgets::pickerInput(
    #         inputId = "find_markers_clust_ID1_tab2",
    #         label = "Select the cluster of interest",
    #         choices = sort(clusters),
    #         multiple = FALSE,
    #         options = list(`actions-box` = TRUE)
    #     )
    #     
    # })
    # 
    # output$find_markers_clust_ID2_tab2_ui <- renderUI ({
    #     
    #     clusters <- req( clusters_single_cell_data_reso_umap_tab2() )
    #     
    #     # Exclude the cluster already selected in the option 1, since the two cluster must be different
    #     clusters <- clusters[ !clusters == req(input$find_markers_clust_ID1_tab2)]
    #     
    #     shinyWidgets::pickerInput(
    #         inputId = "find_markers_clust_ID2_tab2",
    #         label = "Select the cluster(s) to compare",
    #         choices = sort(clusters),
    #         multiple = T,
    #         options = list(`actions-box` = TRUE),
    #         selected = sort(clusters)[1]
    #     )
    #     
    # })
    # 
    # output$find_markers_or_DE_tab2_cluster_ui <- renderUI ({
    #     
    #     clusters <- req( clusters_single_cell_data_reso_umap_tab2() )
    #     
    #     shinyWidgets::pickerInput(
    #         inputId = "find_markers_or_DE_tab2_cluster",
    #         label = "Select the cluster of interest",
    #         choices = sort(clusters),
    #         multiple = FALSE,
    #         options = list(`actions-box` = TRUE)
    #     )
    #     
    # })
    # 
    # output$find_markers_or_DE_tab2_treat1_ui <- renderUI({
    #     
    #     treat <- req( treat_single_cell_data_reso_umap_tab2() )
    #     
    #     shinyWidgets::pickerInput(
    #         inputId = "find_markers_or_DE_tab2_treat1",
    #         label = "Select the sample/treatment of interest",
    #         choices = treat,
    #         multiple = FALSE,
    #         options = list(`actions-box` = TRUE)
    #     )
    #     
    # })
    # 
    # output$find_markers_or_DE_tab2_treat2_ui  <- renderUI({
    #     
    #     treat <- req( treat_single_cell_data_reso_umap_tab2() )
    #     treat <- treat[!treat %in% req(input$find_markers_or_DE_tab2_treat1)]
    #     
    #     shinyWidgets::pickerInput(
    #         inputId = "find_markers_or_DE_tab2_treat2",
    #         label = "Select the sample/treatment of interest",
    #         choices = treat,
    #         multiple = T,
    #         options = list(`actions-box` = TRUE),
    #         selected = treat[1]
    #     )
    #     
    # })
    # 
    # # We change the assay here, so DE analysis and the visualization are performed in the "RNA" assay.
    # single_cell_data_clustered_to_DE_vis <- reactive({
    #     
    #     sc_data <- single_cell_data_clustered()
    #     
    #     sc_data <- NormalizeData(sc_data,
    #                              assay = "RNA",
    #                              normalization.method = "LogNormalize",
    #                              scale.factor = 10000)
    #     
    #     showNotification("Scalling the data of RNA assay",
    #                      duration = NULL,
    #                      id = "tab2_m4")
    #     
    #     all_genes <- rownames(sc_data@assays$RNA)
    #     sc_data <- Seurat::ScaleData(sc_data,
    #                                  assay = "RNA", 
    #                                  features = all_genes)
    #     
    #     on.exit(removeNotification(id = "tab2_m4"), add = TRUE)
    #     
    #     DefaultAssay(sc_data) <- "RNA"
    #     sc_data
    # })
    # 
    # # Identification of markers (tab 2)
    # markers_tab2 <- eventReactive( input$run_ident_markers_tab2, {
    #     
    #     showNotification("Identifing markers or D.E. genes",
    #                      duration = NULL,
    #                      id = "tab2_n1")
    #     
    #     sc_data <- req( single_cell_data_clustered_to_DE_vis() )
    #     
    #     if (input$find_markers_or_DE_tab2 == 0 ) {
    #         
    #         showNotification("Identifing markers or D.E. genes",
    #                          duration = NULL,
    #                          id = "tab2_n1")
    #         
    #         if ( input$find_markers_tab2_opt == 0 ) {
    #             
    #             clusters <- as.numeric(as.character(unique(sc_data@meta.data$seurat_clusters)))
    #             markers_tab2 <- data.frame()
    #             
    #             for (i in clusters) {
    #                 
    #                 markers <- FindConservedMarkers(sc_data,
    #                                                 ident.1 = i,
    #                                                 grouping.var = "treat",
    #                                                 verbose = T,
    #                                                 assay = "RNA" )
    #                 
    #                 markers$cluster <- i
    #                 markers$geneID <- base::rownames(markers)
    #                 
    #                 if ( nrow(markers) == 0) {
    #                     
    #                     markers <- rbind(rep(paste("No markers identified for cluster", i ) , ncol(markers_tab2)))
    #                     colnames(markers) <- colnames(markers_tab2)
    #                     
    #                     
    #                 }
    #                 
    #                 markers_tab2 <- dplyr::bind_rows(markers_tab2, markers)
    #                 
    #             }
    #             
    #             markers_tab2
    #             
    #         } else if ( input$find_markers_tab2_opt == 1 ) {
    #             
    #             markers_tab2 <- FindConservedMarkers(sc_data,
    #                                                  ident.1 = req(input$find_markers_clust_id_tab2),
    #                                                  grouping.var = "treat",
    #                                                  assay = "RNA" )
    #             
    #             markers_tab2$cluster <- input$find_markers_clust_id_tab2
    #             
    #             markers_tab2$geneID <- base::rownames(markers_tab2)
    #             
    #         } else if ( input$find_markers_tab2_opt == 2 ) {
    #             
    #             markers_tab2 <- FindConservedMarkers(sc_data,
    #                                                  ident.1 = req(input$find_markers_clust_ID1_tab2),
    #                                                  ident.2 = req(input$find_markers_clust_ID2_tab2),
    #                                                  grouping.var = "treat",
    #                                                  assay = "RNA" )
    #             
    #             markers_tab2$cluster <- paste0(input$find_markers_clust_ID1_tab2,
    #                                            "_vs_",
    #                                            input$find_markers_clust_ID2_tab2)
    #             
    #             markers_tab2$geneID <- base::rownames(markers_tab2)
    #         }
    #         
    #     } else if (input$find_markers_or_DE_tab2 == 1 ) {
    #         
    #         showNotification("Identifing markers or D.E. genes",
    #                          duration = NULL,
    #                          id = "tab2_n1")
    #         
    #         shinyFeedback::feedbackWarning("find_markers_or_DE_tab2_pvalue",
    #                                        is.na(input$find_markers_or_DE_tab2_pvalue),
    #                                        "Required value")
    #         
    #         req(input$find_markers_or_DE_tab2_pvalue)
    #         
    #         sc_data$celltype.treat <- paste(Idents(sc_data), sc_data$treat, sep = "_")
    #         Idents(sc_data) <- "celltype.treat"
    #         
    #         markers_tab2 <- FindMarkers(sc_data,
    #                                     ident.1 = paste(req(input$find_markers_or_DE_tab2_cluster),
    #                                                     req(input$find_markers_or_DE_tab2_treat1),
    #                                                     sep = "_"),
    #                                     ident.2 = paste(req(input$find_markers_or_DE_tab2_cluster),
    #                                                     req(input$find_markers_or_DE_tab2_treat2),
    #                                                     sep = "_"),
    #                                     test.use = input$find_markers_tab2_test.use,
    #                                     assay = "RNA",
    #                                     verbose = T,
    #         )
    #         
    #         markers_tab2 <- markers_tab2[markers_tab2$p_val_adj < input$find_markers_or_DE_tab2_pvalue, ]
    #         
    #         markers_tab2$comparison <- paste("Cluster", input$find_markers_or_DE_tab2_cluster,
    #                                          input$find_markers_or_DE_tab2_treat1,
    #                                          "vs",
    #                                          input$find_markers_or_DE_tab2_treat2,
    #                                          sep = "_")
    #         
    #         markers_tab2$geneID <- base::rownames(markers_tab2)
    #         
    #     }
    #     
    #     markers_tab2
    #     
    # })
    # 
    # output$markers_tab2_react <- renderReactable({
    #     
    #     markers_tab2 <- req( markers_tab2() )
    #     markers_tab2 <- markers_tab2[ , c( ncol(markers_tab2), (ncol(markers_tab2)-1), 1: (ncol(markers_tab2) -2) ) ]
    #     
    #     on.exit(removeNotification(id = "tab2_n1"), add = TRUE)
    #     
    #     my_reactable(markers_tab2)
    #     
    # })
    # 
    # output$download_markers_tab2 <- downloadHandler(
    #     
    #     filename = function() {
    #         
    #         if (input$find_markers_or_DE_tab2 == 0 ) {
    #             
    #             if ( input$find_markers_tab2_opt == 0 ) {
    #                 
    #                 paste("list_of_markers_int.dataset_all_clusters", ".csv", sep = "")
    #                 
    #             } else if (input$find_markers_tab2_opt == 1) {
    #                 
    #                 paste("list_of_markers_int.dataset", "cluster", input$find_markers_clust_id_tab2, ".csv", sep = "")
    #                 
    #             } else if (input$find_markers_tab2_opt == 2) {
    #                 
    #                 paste("list_of_markers_int.dataset", "cluster", input$find_markers_clust_ID1_tab2, "_vs_" , input$find_markers_clust_ID2_tab2, ".csv", sep = "")
    #                 
    #             }
    #             
    #         } else if (input$find_markers_or_DE_tab2 == 1 ) {
    #             
    #             paste("list_of_markers_int.dataset_DEGs_",
    #                   input$find_markers_or_DE_tab2_treat1,
    #                   "_vs_",
    #                   input$find_markers_or_DE_tab2_treat2,
    #                   "_cluster_",
    #                   input$find_markers_or_DE_tab2_cluster,
    #                   ".csv", sep = "")
    #             
    #         }
    #     },
    #     content = function(file) {
    #         
    #         withProgress(message = "Please wait, preparing the data for download.",
    #                      value = 0.5, {
    #                          
    #                          markers_tab2 <- req( markers_tab2() )
    #                          markers_tab2 <- markers_tab2[ , c( ncol(markers_tab2), (ncol(markers_tab2)-1), 1: (ncol(markers_tab2) -2) ) ]
    #                          
    #                          write.csv(markers_tab2, file, row.names = FALSE)
    #                          
    #                      })
    #     }
    #     
    # )
    # 
    # features_tab2 <- reactive({
    #     
    #     req(input$markers_list_tab2)
    #     req(input$markers_list_header_opt_tab2)
    #     
    #     ext <- tools::file_ext(input$markers_list_tab2$name)
    #     markers_list_file <- input$markers_list_tab2$datapath
    #     
    #     if(input$markers_list_header_opt_tab2 == "") {
    #         shinyFeedback::feedbackWarning("markers_list_header_opt_tab2",
    #                                        TRUE,
    #                                        "Please, inform if the file has a header.")
    #     }
    #     
    #     read_file_markers("tab2_readfile",
    #                       markers_list_file = markers_list_file,
    #                       feed_ID ="markers_list_tab2",
    #                       ext = ext,
    #                       header_opt = input$markers_list_header_opt_tab2)
    # })
    # 
    # observeEvent(input$load_markers_tab2, {
    #     
    #     output$marker_group_selec_tab2 = renderUI({
    #         
    #         pickerInput_markers_group("features_group_tab2", genes = req( features_tab2() ) )
    #         
    #     })
    #     
    #     output$marker_genes_selec_tab2 = renderUI({
    #         
    #         features_group <- req(input$features_group_tab2)
    #         id_choice <- req(input$genes_ids_tab2)
    #         
    #         pickerInput_markers_genes("selected_genes_tab2",
    #                                   genes = req( features_tab2() ),
    #                                   features_group = features_group,
    #                                   id_choice = id_choice)
    #     })
    #     
    # })
    # 
    # filt_features_tab2 <- eventReactive(input$run_heatmap_tab2, {
    #     
    #     features_f <- req( features_tab2() )
    #     features_group <- req(input$features_group_tab2)
    #     
    #     features_f <- dplyr::filter(features_f,
    #                                 Group %in% features_group)
    #     
    #     if (input$filter_genes_q_tab2 == 0) {
    #         
    #         shinyFeedback::feedbackWarning("selected_genes_tab2", is.null(input$selected_genes_tab2), "Please, select one or more genes")
    #         req(input$selected_genes_tab2)
    #         
    #         if (input$genes_ids_tab2 == "ID") {
    #             
    #             features_f <- features_f[features_f$GeneID %in% req(input$selected_genes_tab2), ]
    #             
    #         } else if (input$genes_ids_tab2 == "name") {
    #             
    #             features_f <- features_f[features_f$Name %in% req(input$selected_genes_tab2), ]
    #             
    #         }
    #         
    #     }
    #     
    #     features_f
    #     
    # })
    # 
    # sc_data_av_react_tab2 <- eventReactive(input$run_heatmap_tab2, {
    #     
    #     sc_data_av_tab2  <- Seurat::AggregateExpression(req( single_cell_data_clustered_to_DE_vis() ),
    #                                                     assays = "RNA"#,
    #                                                     #layer = input$slot_selection_heatmap_tab2
    #     )
    #     
    #     sc_data_av_tab2 <- sc_data_av_tab2[[1]]
    #     
    #     on.exit(removeNotification(id = "m7"), add = TRUE)
    #     sc_data_av_tab2
    #     
    # })
    # 
    # observeEvent(input$run_heatmap_tab2, {
    #     
    #     heatmap_n_genes_tab2 <- reactive({
    #         
    #         req(filt_features_tab2())
    #         
    #         heatmap_genes(features = req( filt_features_tab2() ),
    #                       sc_data_av = req( sc_data_av_react_tab2()) )
    #         
    #     })
    #     
    #     heatmap_Height_tab2 <- reactive( 150 + ( 20 * req( heatmap_n_genes_tab2() ) ) )
    #     
    #     heat_map_prep_tab2 <- reactive({
    #         
    #         features <- req(filt_features_tab2())
    #         sc_data_av <- req( sc_data_av_react_tab2() )
    #         features_selec <- as.data.frame(unique(features$GeneID))
    #         
    #         if (nrow(features) == 1) {
    #             
    #             sc_data_av_feat <- sc_data_av[base::rownames(sc_data_av) %in% features_selec[, 1], ]
    #             sc_data_av_feat <- as.data.frame(t(sc_data_av_feat))
    #             base::rownames(sc_data_av_feat) <- features_selec[1, 1]
    #             
    #         } else {
    #             
    #             sc_data_av_feat <- sc_data_av[base::rownames(sc_data_av) %in% features_selec[, 1], ]
    #             
    #         }
    #         
    #         sc_data_av_feat
    #         
    #     })
    #     
    #     output$heat_map_tab2 <- renderPlot({
    #         
    #         heat_map_prep_tab2 <- req( heat_map_prep_tab2() )
    #         
    #         heat_map_prep_tab2 <- as(heat_map_prep_tab2, "sparseMatrix")
    #         
    #         ComplexHeatmap::Heatmap(heat_map_prep_tab2, border = TRUE,
    #                                 rect_gp = gpar(col = "white", lwd = 2),
    #                                 column_title = "Clusters",
    #                                 column_title_side = "bottom",
    #                                 name = "Expression",
    #                                 show_row_dend = T)
    #         
    #     }, height = req( heatmap_Height_tab2()) )
    #     
    #     output$heat_map_ui_tab2 <- renderUI({
    #         
    #         plotOutput("heat_map_tab2", height = req( heatmap_Height_tab2()) )
    #         
    #     })
    #     
    #     output$marker_to_feature_plot_tab2 = renderUI({
    #         
    #         sc_data_av_react <- req( sc_data_av_react_tab2() )
    #         genes_names <- req( filt_features_tab2() )
    #         
    #         genes_names <- genes_names[genes_names$GeneID %in% base::rownames(sc_data_av_react), ]
    #         
    #         if (input$genes_ids_tab2 == "ID") {
    #             
    #             genes_names <- unique(genes_names$GeneID)
    #             
    #         } else if (input$genes_ids_tab2 == "name") {
    #             
    #             genes_names <- unique(genes_names$Name)
    #             
    #         }
    #         
    #         div(class = "option-group",
    #             shinyWidgets::pickerInput(
    #                 inputId = "selected_genes_for_feature_plot_tab2",
    #                 label = "Select the genes to feature plots",
    #                 choices = sort(as.character(genes_names)),
    #                 multiple = TRUE,
    #                 options = list(`actions-box` = TRUE)
    #             ))
    #         
    #     })
    #     
    #     output$p8_down <- downloadHandler(
    #         
    #         filename = function() {
    #             paste("Heatmap", ".", input$p8_format, sep = "")
    #         },
    #         content = function(file) {
    #             
    #             withProgress(message = "Please wait, preparing the data for download.",
    #                          value = 0.5, {
    #                              
    #                              shinyFeedback::feedbackWarning("p8_height", is.na(input$p8_height), "Required value")
    #                              shinyFeedback::feedbackWarning("p8_width", is.na(input$p8_width), "Required value")
    #                              
    #                              height <- as.numeric( req( input$p8_height) )
    #                              width <- as.numeric( req( input$p8_width) )
    #                              res <- as.numeric(input$p8_res)
    #                              
    #                              heat_map_prep_tab2 <- req( heat_map_prep_tab2() )
    #                              heat_map_prep_tab2 <- as(heat_map_prep_tab2, "sparseMatrix")
    #                              p <- ComplexHeatmap::Heatmap(heat_map_prep_tab2, border = TRUE,
    #                                                           rect_gp = gpar(col = "white", lwd = 2),
    #                                                           column_title = "Clusters",
    #                                                           column_title_side = "bottom",
    #                                                           name = "Expression",
    #                                                           show_row_dend = T)
    #                              
    #                              gb = grid.grabExpr(draw(p))
    #                              
    #                              ggplot2::ggsave(file,
    #                                              gb,
    #                                              height=height,
    #                                              width=width,
    #                                              units="cm",
    #                                              dpi=res,
    #                                              bg = "#FFFFFF")
    #                              
    #                          })
    #         }
    #         
    #     )
    #     
    # })
    # 
    # features_selec_tab2 <- eventReactive(input$run_feature_plot_tab2, {
    #     
    #     shinyFeedback::feedbackWarning("selected_genes_for_feature_plot_tab2",
    #                                    is.null(input$selected_genes_for_feature_plot_tab2),
    #                                    "Please, select one or more genes")
    #     req(input$selected_genes_for_feature_plot_tab2)
    #     
    #     features <- req( filt_features_tab2() )
    #     
    #     if (input$genes_ids_tab2 == "ID") {
    #         
    #         features_f <- dplyr::filter(features,
    #                                     GeneID %in% input$selected_genes_for_feature_plot_tab2)
    #         
    #     } else if (input$genes_ids_tab2 == "name") {
    #         
    #         features_f <- dplyr::filter(features,
    #                                     Name %in% input$selected_genes_for_feature_plot_tab2)
    #         
    #     }
    #     
    #     as.character(unique(features_f$GeneID))
    #     
    # })
    # 
    # observeEvent(input$run_feature_plot_tab2, {
    #     
    #     output$feature_plot_tab2 <- renderUI({
    #         
    #         feat_length <- length(req( features_selec_tab2()) )
    #         
    #         plot_output_list <- lapply(1:feat_length, function(i) {
    #             plotname <- paste("plot1_tab2", i, sep="")
    #             plotOutput(plotname, height = 300)
    #         })
    #         
    #         # Convert the list to a tagList - this is necessary for the list of items
    #         # to display properly.
    #         do.call(tagList, plot_output_list)
    #     })
    #     
    #     output$umap2_tab2 <- renderUI({
    #         
    #         feat_length <- length(req( features_selec_tab2()) )
    #         
    #         plot_output_list <- lapply(1:feat_length, function(i) {
    #             plotname <- paste("plot3_tab2", i, sep="")
    #             plotOutput(plotname, height = 300)
    #         })
    #         
    #         do.call(tagList, plot_output_list)
    #     })
    #     
    #     output$run_vln_plot_tab2 <- renderUI({
    #         
    #         feat_length <- length(req( features_selec_tab2()) )
    #         
    #         plot_output_list <- lapply(1:feat_length, function(i) {
    #             plotname <- paste("plot4_tab2", i, sep="")
    #             plotOutput(plotname, height = 300)
    #         })
    #         
    #         do.call(tagList, plot_output_list)
    #     })
    #     
    #     sc_data_tab2 <- req( single_cell_data_clustered_to_DE_vis() )
    #     assay_id <- "RNA"
    #     
    #     for (i in 1:length(req( features_selec_tab2())) ) {
    #         
    #         features <- req( features_selec_tab2() )
    #         local({
    #             
    #             my_i <- i
    #             
    #             plotname <- paste("plot1_tab2", my_i, sep="")
    #             output[[plotname]] <- renderPlot({
    #                 
    #                 p_list <- FeaturePlotSingle(sc_data_tab2,
    #                                             feature = features[my_i],
    #                                             metadata_column = "treat",
    #                                             assay_id = assay_id,
    #                                             pt.size = 0.2,
    #                                             order = TRUE,
    #                                             reduction = "umap",
    #                                             label = T)
    #                 
    #                 groups <- unique(sc_data_tab2@meta.data[, "treat"])
    #                 
    #                 if ( length(groups) > 3 ) {
    #                     
    #                     wrap_plots(p_list,
    #                                guides = 'collect',
    #                                ncol = length(groups) ) +
    #                         plot_annotation(title = paste("Gene name:",
    #                                                       features[my_i]),
    #                                         theme = theme(plot.title = element_text(size = 16,
    #                                                                                 face = "bold",
    #                                                                                 hjust = 0.5)) )
    #                 } else {
    #                     
    #                     wrap_plots(p_list,
    #                                guides = 'collect',
    #                                ncol = length(groups),
    #                                ceiling( (length(groups) / 3  + 1 ) ) ) +
    #                         plot_annotation(title = paste("Gene name:",
    #                                                       features[my_i]),
    #                                         theme = theme(plot.title = element_text(size = 16,
    #                                                                                 face = "bold",
    #                                                                                 hjust = 0.5)) )
    #                 }
    #                 
    #             })
    #             
    #             plotname <- paste("plot3_tab2", my_i, sep="")
    #             output[[plotname]] <- renderPlot({
    #                 
    #                 Seurat::DimPlot(req( single_cell_data_clustered_to_DE_vis() ), reduction = "umap", label = T, pt.size = .1)
    #                 
    #             })
    #             
    #             plotname <- paste("plot4_tab2", my_i, sep="")
    #             output[[plotname]] <- renderPlot({
    #                 
    #                 Seurat::VlnPlot(sc_data_tab2,
    #                                 features = features[my_i],
    #                                 split.by = "treat",
    #                                 assay = "RNA"
    #                 )
    #                 
    #             })
    #             
    #         })
    #         
    #     }
    #     showNotification("Generating additional plots",
    #                      duration = 30,
    #                      id = "m8")
    #     
    #     output$select_genes_add_plot_to_down_tab2_ui <- renderUI({
    #         
    #         filt_features <- req(input$selected_genes_for_feature_plot_tab2)
    #         
    #         div(class = "option-group",
    #             shinyWidgets::pickerInput(
    #                 inputId = "select_genes_add_plot_to_down_tab2",
    #                 label = "Select the genes that you want to download",
    #                 choices = sort(filt_features),
    #                 multiple = TRUE,
    #                 options = list(`actions-box` = TRUE)
    #             ))
    #         
    #     })
    # })
    # 
    # observeEvent(input$start_down_add_plots_tab2, {
    #     
    #     on.exit(removeNotification(id = "m8"), add = TRUE)
    #     
    #     shinyFeedback::feedbackWarning("select_genes_add_plot_to_down_tab2",
    #                                    is.null(input$select_genes_add_plot_to_down_tab2),
    #                                    "Please, select one or more genes")
    #     genes <- req(input$select_genes_add_plot_to_down_tab2)
    #     
    #     withProgress(message = "Please wait, preparing the data for download.",
    #                  value = 0.5, {
    #                      
    #                      if (file.exists("./images") == F ) {
    #                          dir.create('./images')
    #                      }
    #                      
    #                      # Creates new folder to keep the plots organized. Adding date and time helps to avoid overwriting the data
    #                      path_new <- paste0("./images/integrated_sample_plots", format(Sys.time(),'_%Y-%m-%d__%H%M%S'))
    #                      
    #                      dir.create(path_new)
    #                      dir.create( paste0(path_new,"/feature_plots_combined_samples") )
    #                      dir.create( paste0(path_new,"/violin_plots_combined_samples") )
    #                      
    #                      dir.create( paste0(path_new,"/feature_plots") )
    #                      dir.create( paste0(path_new,"/violin_plots") )
    #                      
    #                      sc_data <- req( single_cell_data_clustered_to_DE_vis() )
    #                      assay_id <- "RNA"
    #                      
    #                      for( i in 1:length(genes) ){
    #                          
    #                          minimal <- min( sc_data[[ assay_id ]]@data[ genes[i], ] )
    #                          maximal <- max( sc_data[[ assay_id ]]@data[ genes[i], ] )
    #                          
    #                          p <- suppressMessages(Seurat::FeaturePlot(sc_data,
    #                                                                    features = genes[i],
    #                                                                    layer = input$slot_selection_feature_plot_tab2,
    #                                                                    reduction = "umap") +
    #                                                    scale_colour_gradient2(limits=c(minimal, maximal),
    #                                                                           midpoint = maximal / 2,
    #                                                                           low = "gray80",
    #                                                                           mid = "gold",
    #                                                                           high = "red"))
    #                          
    #                          file <- paste0(path_new, "/feature_plots_combined_samples/", genes[i], ".", req(input$add_p_tab2_feat_format) )
    #                          
    #                          ggplot2::ggsave(file,
    #                                          p,
    #                                          height= req(input$add_p_tab2_feat_height),
    #                                          width=req(input$add_p_tab2_feat_width),
    #                                          units="cm",
    #                                          dpi=as.numeric(req(input$add_p_tab2_feat_res)),
    #                                          bg = "#FFFFFF")
    #                          
    #                          file <- paste0(path_new, "/violin_plots_combined_samples/", genes[i], ".", req(input$add_p_tab2_violin_format) )
    #                          
    #                          p <- Seurat::VlnPlot(sc_data,
    #                                               features = genes[i],
    #                                               layer = req(input$slot_selection_feature_plot_tab2) )
    #                          
    #                          ggplot2::ggsave(file,
    #                                          p,
    #                                          height=req(input$add_p_tab2_violin_height),
    #                                          width=req(input$add_p_tab2_violin_width),
    #                                          units="cm",
    #                                          dpi=as.numeric( req(input$add_p_tab2_violin_res) ),
    #                                          bg = "#FFFFFF")
    #                          
    #                          p_list <- FeaturePlotSingle(sc_data,
    #                                                      feature = genes[i],
    #                                                      metadata_column = "treat",
    #                                                      assay_id = assay_id,
    #                                                      pt.size = 0.2,
    #                                                      order = TRUE,
    #                                                      reduction = "umap",
    #                                                      label = T)
    #                          
    #                          groups <- unique(sc_data@meta.data[, "treat"])
    #                          
    #                          if ( length(groups) > 3 ) {
    #                              
    #                              p2 <-  wrap_plots(p_list,
    #                                                guides = 'collect',
    #                                                ncol = length(groups) ) +
    #                                  plot_annotation(title = paste("Gene name:",
    #                                                                genes[i]),
    #                                                  theme = theme(plot.title = element_text(size = 16,
    #                                                                                          face = "bold",
    #                                                                                          hjust = 0.5)) )
    #                          } else {
    #                              
    #                              p2 <- wrap_plots(p_list,
    #                                               guides = 'collect',
    #                                               ncol = length(groups),
    #                                               ceiling( (length(groups) / 3  + 1 ) ) ) +
    #                                  plot_annotation(title = paste("Gene name:",
    #                                                                genes[i]),
    #                                                  theme = theme(plot.title = element_text(size = 16,
    #                                                                                          face = "bold",
    #                                                                                          hjust = 0.5)) )
    #                          }
    #                          
    #                          file <- paste0(path_new, "/feature_plots/", genes[i], ".", req(input$add_p_tab2_feat_format) )
    #                          
    #                          ggplot2::ggsave(file,
    #                                          p2,
    #                                          height = req(input$add_p_tab2_feat_height),
    #                                          width = (req(input$add_p_tab2_feat_width) * length(groups)),
    #                                          units = "cm",
    #                                          dpi = as.numeric(req(input$add_p_tab2_feat_res)),
    #                                          bg = "#FFFFFF")
    #                          
    #                          file <- paste0(path_new, "/violin_plots/", genes[i], ".", req(input$add_p_tab2_violin_format) )
    #                          
    #                          p2 <- Seurat::VlnPlot(sc_data,
    #                                                features = genes[i],
    #                                                split.by = "treat",
    #                                                assay = "RNA",
    #                                                layer = req(input$slot_selection_feature_plot_tab2) )
    #                          
    #                          ggplot2::ggsave(file,
    #                                          p2,
    #                                          height = req(input$add_p_tab2_violin_height),
    #                                          width = ( req(input$add_p_tab2_violin_width) * length(groups) ),
    #                                          units = "cm",
    #                                          dpi = as.numeric( req(input$add_p_tab2_violin_res) ),
    #                                          bg = "#FFFFFF")
    #                          
    #                      }
    #                      
    #                  })
    #     
    #     showNotification("All plots were downloaded!",
    #                      duration = 15,
    #                      id = "")
    #})
    
    #####################
    ### Stacked plots ###
    #####################
    
    stacked_violin_Server("stacked1")
    
}
