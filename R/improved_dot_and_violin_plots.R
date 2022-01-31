stacked_violin_UI <- function(id) {
    ns <- NS(id)
    
    tagList(
        
        fluidRow (
            
            shinyFeedback::useShinyFeedback(), # include shinyFeedback
            
            column(3,
                   my_withSpinner( uiOutput( ns("load_rds_ui") )),
                   
                   actionButtonInput( ns("load_data"),
                                      "Load data"),
            ),
            conditionalPanel(
                condition = "input.load_data > 0", ns = ns,
                
                column(width = 3,
                       
                       fileInput_markers_list( ns("markers_list") ),
                       div(class = "option-group",
                           selectInput(ns("markers_list_header_opt"),
                                       "Does your file have a header?",
                                       choices = c("", "Yes", "No"),
                                       multiple = FALSE,
                                       selectize = TRUE,
                                       selected = NULL)),
                       
                       div(class = "option-group",
                           actionButtonInput( ns("load_markers"),
                                              HTML("Load markers"))),
                       
                ),
                
                column(3,
                       
                       conditionalPanel(
                           condition = "input.load_markers > 0", ns = ns,
                           
                           my_withSpinner( uiOutput( ns('marker_group_selec') ) ),
                           
                           define_if_use_all_genes_or_select( ns("filter_genes_q") ),
                           
                           conditionalPanel(
                               condition = "input.filter_genes_q == 0", ns = ns,
                               
                               define_what_id_to_use( ns("genes_ids") ),
                               
                               my_withSpinner( uiOutput( ns('marker_genes_selec2') ) )
                               
                           ),
                           
                       ),
                       
                       # ),
                       
                ),
                
                column(3,
                       actionButtonInput(ns("start_violin"),
                                         "Generate plots")
                ),
            ),
        ),
        
        fluidRow(
            conditionalPanel(
                condition = "input.start_violin > 0", ns = ns,
                
                titlePanel("Stacked Violin plot"),
                column(10,
                       my_withSpinner( uiOutput( ns("stacked_violin_ui") ) )
                ),
                column(2,
                       
                       div(class = "option-group",
                           radioButtons(ns("order_plots_xaxis"),
                                        "Set the x-axis order?",
                                        list("No" = 0,
                                             "Yes" = 1),
                                        selected = 0)),
                       
                       conditionalPanel(
                           condition = "input.order_plots_xaxis == 1", ns = ns,
                           
                           
                           my_withSpinner( uiOutput( ns("plots_xaxis_ui") ) ),
                           
                           actionButtonInput( ns("order_plot_button"),
                                              "Update plot"),
                       ),
                       
                ),
                column(2,
                       div(class = "down-group",
                           numericInput(ns("p1_height"), "Height (cm)", step=0.5, value=10),
                           numericInput(ns("p1_width"), "Width (cm)", step=0.5, value=20),
                           selectInput(ns("p1_res"),
                                       "Resolution (DPI)",
                                       choices=c("100","200","300","400","500","600"),
                                       selected=100),
                           selectInput(ns("p1_format"),
                                       "File type",
                                       choices=list("png" = "png",
                                                    "tiff" = "tiff",
                                                    "jpeg" = "jpeg",
                                                    "pdf" = "pdf",
                                                    "svg" = "svg"),
                                       selected="png",
                                       multiple=FALSE,
                                       selectize=TRUE)
                       ),
                       downloadButton(ns("p1_down"), HTML("Download Plot") )
                )
            )
        ),
        
        fluidRow(
            conditionalPanel(
                condition = "input.start_violin > 0", ns = ns,
                
                titlePanel("Multi-genes Dot plot"),
                column(10,
                       my_withSpinner( uiOutput( ns("stacked_dotpot_ui") ) )
                ),
                # download plot
                column(2,
                       div(class = "down-group",
                           numericInput(ns("p2_height"), "Height (cm)", step=0.5, value=10),
                           numericInput(ns("p2_width"), "Width (cm)", step=0.5, value=20),
                           selectInput(ns("p2_res"),
                                       "Resolution (DPI)",
                                       choices=c("100","200","300","400","500","600"),
                                       selected=100),
                           selectInput(ns("p2_format"),
                                       "File type",
                                       choices=list("png" = "png",
                                                    "tiff" = "tiff",
                                                    "jpeg" = "jpeg",
                                                    "pdf" = "pdf",
                                                    "svg" = "svg"),
                                       selected="png",
                                       multiple=FALSE,
                                       selectize=TRUE)
                       ),
                       downloadButton(ns("p2_down"), HTML("Download Plot") )
                )
            )
        ),
        
    )
    
}

stacked_violin_Server <- function(id) {
    moduleServer(id, function(input, output, session) {
        
        output$load_rds_ui <- renderUI({
            
            rds_list <- list.files('./RDS_files/', pattern = "*.rds")
            
            div(class = "option-group",
                shinyWidgets::pickerInput(
                    inputId = session$ns("load_rds"),
                    label = "Select the file containing the clustered data",
                    choices = sort(rds_list),
                    multiple = FALSE,
                    options = list(`actions-box` = TRUE)
                ))
            
        })
        
        sc_data <- eventReactive(input$start_violin, {
            
            load_rds <- req(input$load_rds)
            sc_data <- readRDS( paste0("./RDS_files/", load_rds) )
            
            validate(need("seurat_clusters" %in% colnames(sc_data@meta.data),
                          "Clustering was not detected. Was the data clustered in the other tabs?", ""))
            
            sc_data
            
        })
        
        observeEvent(input$load_data, {
            
            output$list_of_genes_ui <- renderUI({
                
                load_rds <- req(input$load_rds)
                sc_data <- readRDS( paste0("./RDS_files/", load_rds) )
                
                genes <- base::rownames(sc_data)
                
                div(class = "option-group",
                    selectInput(inputId = session$ns("list_of_genes"),
                                label = "Select the genes to shown in the stacked plot",
                                choices = sort(unique(genes)),
                                multiple = TRUE,
                                selectize = TRUE
                    )
                )
                
            })
            
        })
        
        features <- reactive({
            
            req(input$markers_list)
            ext <- tools::file_ext(input$markers_list$name)
            markers_list_file <- input$markers_list$datapath
            
            if(input$markers_list_header_opt == "") {
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
        
        expressed_genes <- reactive( {
            
            load_rds <- req(input$load_rds)
            sc_data <- readRDS( paste0("./RDS_files/", load_rds) )
            
            genes <- base::rownames(sc_data@assays$RNA)
            genes
        })
        
        observeEvent(input$load_markers, {
            
            # Painel that will apper after loading the list of markers containing the filter options
            output$marker_group_selec = renderUI({
                
                df <- features()
                groups_df <- unique(df$Group)
                
                div(class = "option-group",
                    selectInput(inputId = session$ns("features_group"),
                                label = "Select the group of markers to test",
                                choices = sort(groups_df),
                                #choices = as.list( c(0:9)),
                                multiple = TRUE,
                                selectize = TRUE,
                                selected = sort(groups_df)[1]
                    )
                )
                
            })
            
            output$marker_genes_selec2 = renderUI( {
                
                req(input$features_group)
                
                
                features_group <- input$features_group
                id_choice <- req(input$genes_ids)
                
                df <- features()
                
                genes_names <- dplyr::filter(df,
                                             Group %in% features_group)
                
                if (id_choice == "ID") {
                    
                    genes_names <- unique(genes_names$GeneID)
                    
                } else if (id_choice == "name") {
                    
                    genes_names <- unique(genes_names$Name)
                    
                }
                
                genes_names <- genes_names[!is.na(genes_names)]
                
                div(class = "option-group",
                    selectInput(inputId = session$ns("selected_genes"),
                                label = "Select the genes to show",
                                choices = sort(genes_names),
                                #choices = as.list( c(0:9)),
                                multiple = TRUE,
                                selectize = TRUE
                    )
                )
            })
            
        })
        
        # Filter accordly with the parameters selected by the user
        filt_features <- eventReactive(input$start_violin, {
            
            showNotification("Generating the plot. Please wait, it might take a few seconds.",
                             duration = NULL,
                             id = "p1")
            
            req(input$filter_genes_q)
            req(input$genes_ids)
            
            features_f <- features()
            features_group <- req(input$features_group)
            
            features_f <- dplyr::filter(features_f,
                                        Group %in% features_group)
            
            ## Check if the list of genes exist in the dataset, to avoid a crash
            expressed_genes <- req( expressed_genes() )
            features_f <- features_f[features_f$GeneID %in% expressed_genes, ]
            
            if (input$filter_genes_q == 0) {
                
                if (input$genes_ids == "ID") {
                    
                    features_f <- features_f[features_f$GeneID %in% req(input$selected_genes), ]
                    
                } else if (input$genes_ids == "name") {
                    
                    features_f <- features_f[features_f$Name %in% req(input$selected_genes), ]
                    
                }
                
            }
            features_f <- features_f[!duplicated(features_f$GeneID), ]
            
            features_f$Group <- base::factor(x = features_f$Group,
                                             levels = input$features_group)
            features_f <- features_f[order(features_f$Group), ]
            
            features_f
            
        })
        
        selected_genes <- eventReactive( input$start_violin, {
            
            gene_names <- filt_features()
            # 
            if (input$filter_genes_q == 0) {
                
                gene_names <- gene_names[order( match(input$selected_genes, gene_names$GeneID) ), ]
                gene_names <- gene_names$GeneID
                
            } else {
                
                gene_names <- gene_names$GeneID
            }
            
        })
        
        output$plots_xaxis_ui <- renderUI({
            
            sc_data <- sc_data()
            clusters <- unique(sc_data@meta.data$seurat_clusters)
            
            div(class = "option-group",
                selectInput(inputId = session$ns("selected_order"),
                            label = "Select the ordered of the clusters",
                            choices = sort(clusters),
                            #choices = as.list( c(0:9)),
                            multiple = TRUE,
                            selectize = TRUE
                )
            )
            
        })
        
        stacked_violin_Height <- reactive( 150 + ( 50 * length( selected_genes() ) ) )
        
        selected_order <- eventReactive(input$order_plot_button,{
            
            showNotification("Updating the plot. Please wait, it might take a few seconds.",
                             duration = 40,
                             id = "p2")
            
            sc_data <- sc_data()
            clusters <- unique(sc_data@meta.data$seurat_clusters)
            
            selected_order <- req(input$selected_order)
            
            shinyFeedback::feedbackWarning("selected_order",
                                           length(selected_order) != length(clusters),
                                           "Please select the order of all clusters")
            
            validate(need( length(selected_order) == length(clusters),
                           message = "Please select the order of all clusters",
                           label = "selected_order"))
            
            selected_order
            
        })
        
        output$stacked_violin <- renderPlot({
            
            sc_data <- sc_data()
            selected_genes <- selected_genes()
            
            if ( input$order_plots_xaxis == 1 ) {
                
                req(selected_order())
                sc_data@meta.data$seurat_clusters2 <- base::factor(x = sc_data@meta.data$seurat_clusters,
                                                                   levels = selected_order())
                
            } else {
                
                sc_data@meta.data$seurat_clusters2 <- sc_data@meta.data$seurat_clusters
                
            }
            
            stacked_violin(name = "stacked1",
                           rds_file = sc_data,
                           selected_genes = selected_genes,
                           pt.size = 0)
            
        }, height = function() { stacked_violin_Height() })
        
        output$stacked_violin_ui <- renderUI({
            
            plotOutput(session$ns("stacked_violin"), height = stacked_violin_Height())
            
        })
        
        stacked_dotpot_Height <- reactive( 150 + ( 30 * length( selected_genes() ) ) )
        
        output$stacked_dotpot <- renderPlot({
            
            removeNotification("p1")
            
            sc_data <- sc_data()
            selected_genes <- selected_genes()
            
            if ( input$order_plots_xaxis == 1 ) {
                
                req(selected_order())
                order_levels <- selected_order()
                
                sc_data@meta.data$seurat_clusters2 <- base::factor(x = sc_data@meta.data$seurat_clusters, levels = order_levels )
                
            } else {
                
                sc_data@meta.data$seurat_clusters2 <- sc_data@meta.data$seurat_clusters
                
            }
            
            Seurat::DotPlot(object = sc_data,
                            features = rev(selected_genes),
                            cols = c("Darkblue", "red"),
                            group.by = "seurat_clusters2",
                            assay = "RNA"
            ) +
                scale_x_discrete(position = "top") +
                xlab("")+
                ylab("") +
                ggtitle("") +
                theme(
                    #legend.position = "",
                    #  axis.title.x=element_blank(),
                    #  axis.title.y=element_blank(),
                    axis.text.y = element_text(size = rel(1)),
                ) +
                coord_flip()
            
        }, height = function() { stacked_dotpot_Height() })
        
        output$stacked_dotpot_ui <- renderUI({
            
            showNotification("Generating the plot. Please wait, it might take a few seconds.",
                             duration = 40,
                             id = "p1")
            
            plotOutput(session$ns("stacked_dotpot"), height = stacked_dotpot_Height())
            
        })
        
        ### Downloading the plots
        
        output$p1_down <- downloadHandler(
            
            filename = function() {
                paste("Stacked_violin_plot", ".", input$p1_format, sep = "")
            },
            content = function(file) {
                
                withProgress(message = "Please wait, preparing the data for download.",
                             value = 0.5, {
                                 
                                 height <- as.numeric(input$p1_height)
                                 width <- as.numeric(input$p1_width)
                                 res <- as.numeric(input$p1_res)
                                 
                                 sc_data <- sc_data()
                                 selected_genes <- selected_genes()
                                 
                                 if ( input$order_plots_xaxis == 1 ) {
                                     
                                     req(selected_order())
                                     sc_data@meta.data$seurat_clusters2 <- base::factor(x = sc_data@meta.data$seurat_clusters,
                                                                                        levels = selected_order())
                                     
                                 } else {
                                     
                                     sc_data@meta.data$seurat_clusters2 <- sc_data@meta.data$seurat_clusters
                                     
                                 }
                                 
                                 p <- stacked_violin(name = "stacked1",
                                                     rds_file = sc_data,
                                                     selected_genes = selected_genes,
                                                     pt.size = 0)
                                 
                                 ggsave(file,
                                        p,
                                        height=height,
                                        width=width,
                                        units="cm",
                                        dpi=res)
                                 
                             })
                
            }
            
        )
        
        output$p2_down <- downloadHandler(
            
            filename = function() {
                paste("Multigenes_dot_plot", ".", input$p2_format, sep = "")
            },
            content = function(file) {
                
                withProgress(message = "Please wait, preparing the data for download.",
                             value = 0.5, {
                                 
                                 sc_data <- sc_data()
                                 selected_genes <- selected_genes()
                                 
                                 if ( input$order_plots_xaxis == 1 ) {
                                     
                                     req(selected_order())
                                     sc_data@meta.data$seurat_clusters2 <- base::factor(x = sc_data@meta.data$seurat_clusters,
                                                                                        levels = selected_order())
                                     
                                 } else {
                                     
                                     sc_data@meta.data$seurat_clusters2 <- sc_data@meta.data$seurat_clusters
                                     
                                 }
                                 
                                 p <- Seurat::DotPlot(object = sc_data,
                                                      features = rev(selected_genes),
                                                      cols = c("Darkblue", "red"),
                                                      group.by = "seurat_clusters2",
                                                      assay = "RNA"
                                 ) +
                                     scale_x_discrete(position = "top") +
                                     xlab("")+
                                     ylab("") +
                                     ggtitle("") +
                                     theme(
                                         #legend.position = "",
                                         # axis.title.x=element_blank(),
                                         #axis.title.y=element_blank(),
                                         axis.text.y = element_text(size = rel(1))) +
                                     coord_flip()
                                 
                                 height <- as.numeric(input$p2_height)
                                 width <- as.numeric(input$p2_width)
                                 res <- as.numeric(input$p2_res)
                                 
                                 ggsave(file,
                                        p,
                                        height=height,
                                        width=width,
                                        units="cm",
                                        dpi=res)
                                 
                             })
                
            }
            
        )
        
    })
    
}
