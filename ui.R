# SolusCell web app
# Version 1.0
set.seed(1407)

suppressMessages( require("shiny") )
suppressMessages( require("shinyWidgets") )
suppressMessages( require("DT") )
suppressMessages( require("reactable") )
suppressMessages( require("rclipboard") )
suppressMessages( require("shinycssloaders") )
suppressMessages( require("shinyFeedback") )

source("R/ui_functions.R")
source("R/improved_dot_and_violin_plots.R")

function(request) {
    
    ## This set the app to properly work with docker
    if (dir.exists('/app/user_work')) {
        setwd('/app/user_work')
    }
    
    fluidPage(
        
        shinyFeedback::useShinyFeedback(),
        
        ####################
        #### Custom CSS ####
        ####################
        
        tags$head(
            tags$link(rel = "stylesheet", type = "text/css", href = "css/SolusCell_web_app.css")
        ),
        
        tags$style("* {font-family: 'Montserrat', sans-serif !important;}"),
        
        rclipboardSetup(),
        
        # Sets the background color of the application
        setBackgroundColor(
            color = c("#F7FBFF", "#f2f9fc", "#f2fbff", "#d9f0fa", "#ccf0ff"),
            gradient = "radial",
            direction = c("top", "left")),
        
        ###############################
        # SolusCell web app interface #
        ###############################
        
        titlePanel( column( 12, a(
            img(src = "SolusCell_logo.png",
                width = "25%",
                align = "left") ) ),
            windowTitle="MyPage"
        ),
        
        br(),
        
        tabsetPanel(
            
            ######################################
            ###           One sample           ###
            ######################################
            
            tabPanel("One sample",
                     fluidRow(
                         br(),
                         p("Choose the sample to be analyzed and the initial requirements to load the data."),
                         br(),
                         
                         column(width = 2,
                                div(class = "option-group",
                                    radioButtons("sample_tab1_options",
                                                 "Run a new analysis or read a previously saved file?",
                                                 choices = list("Run a new analysis" = 0,
                                                                "Load clustered data" = 1),
                                                 selected = c(1) )),
                                
                                conditionalPanel (
                                    condition = "input.sample_tab1_options == 0",
                                    
                                    my_withSpinner(
                                        uiOutput('select_sample_tab1'))
                                    
                                ),
                         ),
                         
                         conditionalPanel (
                             condition = "input.sample_tab1_options == 1",
                             
                             column(width = 3,
                                    my_withSpinner( uiOutput("select_sample_tab1_rds_ui") ),
                             ),
                             
                             column( width = 3,
                                     actionButtonInput("load_10X_rds",
                                                       HTML("Load clustered data"))
                             ),
                         ),
                         
                         conditionalPanel (
                             condition = "input.sample_tab1_options == 0",
                             
                             column(width = 2,
                                    input_project_name("proj_name")
                             ),
                             
                             # column(width = 2,
                             #        div(class = "option-group",
                             #            numericInput("min_cells",
                             #                         label = "Min. number of cells expressing a gene for the gene to be included",
                             #                         value = 3))
                             # ),
                             
                             # column(width = 2,
                             #        div(class = "option-group",
                             #            numericInput("min_features",
                             #                         label = "Min. number of genes a cell must express to be included",
                             #                         value = 200))
                             # ),
                             
                             column( width = 2,
                                     textInput_regex("mito_regex") ,
                                     
                             ),
                             
                             
                             column( width = 3,
                                     actionButtonInput("load_10X",
                                                       HTML("Load data of the selected sample"))
                             ),
                             
                         ),
                         
                     ), # fluid row
                     
                     # If data is already clustered, this options doesn't appear
                     conditionalPanel (
                         condition = "input.sample_tab1_options == 0",
                         
                         conditionalPanel (
                             condition = "input.load_10X != 0",
                             
                             fluidRow(
                                 titlePanel("Genes identified as mitochondrial based on the informed identifier"),
                                 column(12,
                                        my_withSpinner( verbatimTextOutput("target_genes_mitho") )
                                 )
                             ), # Fluid row
                         ),
                         
                         fluidRow(
                             titlePanel("Screening plot to define filtering parameters"),
                             br(),
                             p("Use the plot below to evaluate if defining more restrictive parameters is necessary, then select values for the boxes on the right. Next, click on \"Show/update plot of filtered data\"."),
                             br(),
                             column(width = 10,
                                    my_withSpinner( plotOutput("VlnPlot") )
                             ),
                             column(width = 2,
                                    div(class = "option-group",
                                        numericInput("min_count",
                                                     label = "Keep cells that expressed at least this number of genes",
                                                     value = 1),
                                        numericInput("max_count",
                                                     label = "Exclude cells that expressed more than this number of genes",
                                                     value = 10000)),
                                    
                                    numericInput_max_mito_perc("max_mito_perc", value = 100),
                                    actionButtonInput("run_vinplot",
                                                      HTML("Show/update plot <br/> of filtered data"))
                             )
                         ), # fluid row
                         conditionalPanel(
                             condition = "input.run_vinplot!= 0",
                             
                             fluidRow (
                                 titlePanel("Screening plot showing the remaining cells after filtering"),
                                 column(width = 10,
                                        my_withSpinner( plotOutput("VlnPlot_filt") )
                                 ),
                                 # download plot
                                 column(2,
                                        numericInput_plot_height("p1_height"),
                                        numericInput_plot_width("p1_width", value=25),
                                        selectInput_plot_res("p1_res"),
                                        selectInput_plot_format("p1_format"),
                                        downloadButton("p1_down", HTML("Download Plot") )
                                 )
                             ), # fluid row
                             
                         ), # Ends conditional
                         
                         fluidRow(
                             titlePanel("Normalization and dimension reduction analysis (PCA)"),
                             br(),
                             #p("Note that only the most variable genes are used in the dimension reduction step (PCA)."),
                             #br(),
                             column(width = 3,
                                    select_norm_methods("normaliz_method")
                             ),
                             conditionalPanel(
                                 condition = "input.normaliz_method == 0",
                                 column(width = 4,
                                        numericInput_scale_factor("scale_factor"),
                                        radioInput_most_var_method("most_var_method"),
                                        numericInput_var_genes("n_of_var_genes")
                                 ),
                             ),
                             column(width = 3,
                                    actionButtonInput("run_pca",
                                                      HTML("Run the PCA analysis"))
                             )
                             
                         ), # Fluid row,
                         
                         conditionalPanel(
                             condition = "input.run_pca!= 0",
                             
                             fluidRow(
                                 titlePanel("Elbow plot showing the variation explained by each Principal Component"),
                                 br(),
                                 p("Choose the number of components that explain the most variation."),
                                 br(),
                                 column(width = 6,
                                        my_withSpinner( plotOutput("n_of_PCAs") )
                                 ),
                                 # Plot download
                                 column(2,
                                        numericInput_plot_height("p2_height", value=12),
                                        numericInput_plot_width("p2_width", value=18),
                                        selectInput_plot_res("p2_res"),
                                        selectInput_plot_format("p2_format"),
                                        downloadButton("p2_down", HTML("Download Plot"))
                                        
                                 ),
                                 
                                 column(width = 3,
                                        numericInput_n_of_PCs("n_of_PCs")
                                 )
                                 
                             ) # Fluid row
                         ), # Ends conditional
                         
                     ),
                     
                     conditionalPanel(
                         condition = "input.load_10X_rds == 0", # <<<<<<
                         
                         fluidRow(
                             
                             conditionalPanel (
                                 condition = "input.sample_tab1_options == 0",
                                 titlePanel("Clustering"),
                                 
                                 br(),
                                 p("Be aware that this parameter is central in the cluster definition. It is recommended to try different values and define the most appropriate value according to the expectations of the cell populations present in the sample."),
                                 p("Quoting from", a(tags$a(href="https://satijalab.org/seurat/archive/v1.4/pbmc3k_tutorial.html", "Seurat's tutorial", target="_blank")), ":", em("\"We find that setting this parameter between 0.6-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets.\"")),
                                 br(),
                                 column(width = 3,
                                        numericInput_resolution_clust("resolution_clust")
                                 ),
                                 column(width = 3,
                                        actionButtonInput("run_clustering",
                                                          HTML("Run the clustering analysis"))
                                 )
                                 
                             ), # Fluid row
                         ),
                     ),
                     
                     conditionalPanel(
                         condition = "input.run_clustering!= 0 || input.load_10X_rds!= 0",
                         
                         fluidRow(
                             titlePanel("Clustering plots"),
                             
                             column(width = 6,
                                    div(id = "remove",
                                        my_withSpinner( plotOutput("umap") )
                                    )
                             ),
                             column(width = 6,
                                    div(id = "remove",
                                        my_withSpinner( plotOutput("tSNE") )
                                    )
                             )
                             
                         ), # Fluid row
                         fluidRow(
                             br(),
                             # download plot
                             column(width = 3,
                                    div(class = "down-group",
                                        radioButtons("p3_down_opt",
                                                     "Select the plot to download",
                                                     choices = list("UMAP" = "UMAP",
                                                                    "t-SNE" = "t-SNE"),
                                                     selected = c("UMAP")))
                             ),
                             column(2,
                                    numericInput_plot_height("p3_height", value=10),
                                    numericInput_plot_width("p3_width", value=18)
                             ),
                             column(2,
                                    selectInput_plot_res("p3_res"),
                                    selectInput_plot_format("p3_format")
                             ),
                             column(2,
                                    div(class = "down-group",
                                        downloadButton("p3_down", HTML("Download Plot")))
                             ),
                         ),
                         
                         fluidRow(
                             titlePanel("Number of cells per cluster"),
                             column(2,
                                    div(class = "down-group", id = "remove",
                                        my_withSpinner( tableOutput("cluster_size") )
                                    ) )
                         ), # Fluid row
                         
                     ), # Ends conditional
                     
                     conditionalPanel (
                         condition = "input.sample_tab1_options == 0",
                         
                         fluidRow(
                             titlePanel("Excluding or selecting clusters for reanalysis"),
                             br(),
                             p("Sometimes, it is helpful to exclude or select the clusters that are of interest. After excluding or selecting the cells of interest, it is recommended to repeat the clustering step using only the subset."),
                             p("After selecting the clusters, click on the blue button (Reanalyze after selection/exclusion of clusters). Asc-Seurat will run the analyses of the new subset until the PCA step.", strong("Then, users need to set the new number of components using the elbow plot (above) and click on the button \"Run the clustering analysis\" again.")),
                             br(),
                             column(3,
                                    radioButtons_filter_clusters("filter_clusters")
                             ),
                             
                             conditionalPanel(
                                 condition = "input.filter_clusters != 0",
                                 
                                 column(3,
                                        radioButtons_filter_clusters_opt("filter_clusters_opt")
                                 ),
                                 
                                 column( 3,
                                         div(class = "option-group",
                                             my_withSpinner( uiOutput("cluster_list_ui") ))
                                 ),
                                 
                                 column( 3,
                                         div(class = "option-group",
                                             actionButtonInput("rerun_after_filtering",
                                                               HTML("Reanalyze after selection/exclusion of clusters")))
                                 ),
                                 
                             ) # ends conditional
                         ),
                         
                         fluidRow(
                             br(),
                             h3("Saving the processed data for the trajectory inference analysis"),
                             p("The button below allows for saving the processed data in a file that can be used for the pseudo-time analysis."),
                             #p("Asc-Seurat will save only the most recently processed data."),
                             br(),
                             p(strong("Note:", "The processed data needs to be saved in the folder", code("RDS_files/"), "so it can be load automatically in the tab \"Trajectory inference.\"" )),
                             column(6,
                                    br(),
                                    downloadButton("downloadRDS",
                                                   "Download the processed data to use in the trajectory inference analysis",
                                                   class = "down_butt")
                             )
                         ),
                     ),
                     
                     fluidRow(
                         titlePanel("Identification of markers / Differential expression analysis"),
                         column(3,
                                radioButtons_find_markers("find_markers_tab1")
                         )
                     ),
                     
                     conditionalPanel(
                         condition = "input.find_markers_tab1 == 1",
                         fluidRow(
                             br(),
                             p("Check the Seurat's vignettes", a(tags$a(href="https://satijalab.org/seurat/v3.2/de_vignette.html", "here", target="_blank")), "and the function's manual",  a(tags$a(href="https://www.rdocumentation.org/packages/Seurat/versions/3.1.4/topics/FindAllMarkers", "here",target="_blank")),  "to see the options for each parameters of this section. "),
                             
                             br(),
                             column(3,
                                    div(class = "option-group",
                                        radioButtons("find_markers_tab1_opt",
                                                     "Select the analysis to perform?",
                                                     choices = list("Identify markers for all clusters" = 0,
                                                                    "Identify markers for one specific cluster" = 1,
                                                                    "Identify markers distinguishing a cluster from other(s) cluster(s)" = 2),
                                                     selected = c(0)))),
                             conditionalPanel (
                                 condition = "input.find_markers_tab1_opt == 1",
                                 
                                 column(3,
                                        div(class = "option-group",
                                            
                                            my_withSpinner( uiOutput("find_markers_clust_id_tab1_ui") ))),
                                 
                             ),   # ends conditional
                             
                             conditionalPanel (
                                 condition = "input.find_markers_tab1_opt == 2",
                                 
                                 column(3,
                                        div(class = "option-group",
                                            my_withSpinner( uiOutput("find_markers_clust_ID1_tab1_ui") ),
                                            my_withSpinner( uiOutput("find_markers_clust_ID2_tab1_ui") ))),
                                 
                             ),    # ends conditional
                             
                             column(3,
                                    div(class = "option-group",
                                        numericInput("find_markers_tab1_logfc_threshold",
                                                     label = "Select the (log) fold change threshold",
                                                     value = 0.25),
                                        numericInput("find_markers_tab1_min_pct",
                                                     label = "Select the minimal percentage of cells expressing a gene for this gene to be tested (0 = 0%, 1 = 100%)",
                                                     value = 0.1,
                                                     step = 0.01),
                                        pickerInput_find_markers_test("find_markers_tab1_test.use"),
                                        conditionalPanel(
                                            condition = "input.find_markers_tab1_opt != 0",
                                            radioButtons("find_markers_tab1_filt_pos",
                                                         "Filter only positive markers?",
                                                         choices = c("yes" = "TRUE",
                                                                     "no" = "FALSE"),
                                                         selected = "TRUE")),
                                        numericInput("find_markers_tab1_return_thresh",
                                                     label = "Select the (adjusted) p-value threshold",
                                                     value = "0.05",
                                                     step = 0.01),
                                    ),
                                    
                                    actionButtonInput("run_ident_markers_tab1",
                                                      HTML("Search for markers/D.E. genes"))),
                         ), # ends conditional
                     ), # ends fluidrow
                     
                     conditionalPanel(
                         condition = "input.find_markers_tab1 == 1",
                         fluidRow(
                             titlePanel("List of markers or differentially expressed genes"),
                             conditionalPanel(
                                 condition = "input.run_ident_markers_tab1 != 0",
                                 column(12,
                                        my_withSpinner( reactableOutput("markers_tab1_react")) ),
                             ),
                             br(),
                             column(4,
                                    downloadButton("download_markers_tab1",
                                                   "Download the list of markers or D.E. genes"))
                         ), # ends fluidrow
                     ), # ends conditional
                     
                     fluidRow(
                         titlePanel("Expression of markers"),
                         br(),
                         p("In this section, users can visualize the gene expression of selected genes. Start by loading a", strong("CSV"), "file with at least two columns. The first column must be the gene ID, and the second is a grouping variable (e.g., the cell type name or cluster number). A third column can be used to store the common name of the gene but is optional."),
                         br(),
                         column(width = 3,
                                fileInput_markers_list("markers_list"),
                                div(class = "option-group",
                                    selectInput("markers_list_header_opt",
                                                "Does your file have a header?",
                                                choices = c("", "Yes", "No"),
                                                multiple = FALSE,
                                                selectize = TRUE,
                                                selected = NULL)),
                                div(class = "option-group",
                                    actionButtonInput("load_markers",
                                                      HTML("Load markers")))),
                         column(width = 2,
                                conditionalPanel(
                                    condition = "input.load_markers > 0",
                                    my_withSpinner( uiOutput('marker_group_selec') ),
                                    define_if_use_all_genes_or_select("filter_genes_q"),
                                    
                                    conditionalPanel(
                                        condition = "input.filter_genes_q == 0",
                                        define_what_id_to_use("genes_ids")
                                    ))),
                         conditionalPanel(
                             condition = "input.filter_genes_q == 0",
                             
                             column(width = 3,
                                    my_withSpinner( uiOutput('marker_genes_selec') )
                             )
                         ),
                         column(width = 2,
                                conditionalPanel(
                                    condition = "input.load_markers > 0",
                                    radioButtons_slot_selection_heatmap("slot_selection_heatmap"))),
                         column(width = 2,
                                conditionalPanel(
                                    condition = "input.load_markers > 0",
                                    actionButtonInput("run_heatmap",
                                                      HTML("Show heatmap"))))
                     ), # Fluid row
                     
                     #heat map
                     conditionalPanel(
                         condition = "input.run_heatmap > 0",
                         
                         fluidRow(
                             br(),
                             titlePanel("Heatmap"),
                             br(),
                             strong("Note that the scale of colors of the heatmap is adjusted based on the expression of the selected genes."),
                             br(),
                             column(width = 10,
                                    my_withSpinner( uiOutput("heat_map_ui"))),
                             column(2,
                                    numericInput_plot_height("p4_height", value=15),
                                    numericInput_plot_width("p4_width", value=20),
                                    selectInput_plot_res("p4_res"),
                                    selectInput_plot_format("p4_format"),
                                    downloadButton("p4_down", HTML("Download Plot"))),
                             
                             column(width = 10,
                                    verbatimTextOutput("test"))
                             
                         ),
                         fluidRow(
                             titlePanel("Visualization of gene expression of each cell and additional plots"),
                             br(),
                             p("Only genes selected for the heatmap can be selected for the additional plots."),
                             br(),
                             column(width = 3,
                                    my_withSpinner( uiOutput("marker_to_feature_plot") )),
                             column(width = 2,
                                    radioButtons_slot_selection_feature_plot("slot_selection_feature_plot")),
                             column(width = 3,
                                    actionButtonInput("run_feature_plot",
                                                      HTML("Show the expression of genes at the cell level")))
                         )
                     ), # Ends conditional
                     
                     conditionalPanel (
                         condition = "input.run_feature_plot > 0",
                         fluidRow(
                             titlePanel("Feature plots"),
                             column(width = 4,
                                    my_withSpinner( uiOutput("feature_plot") )),
                             column(width = 4,
                                    my_withSpinner( uiOutput("feature_plot_dark", height = "300px") )),
                             column(width = 4,
                                    my_withSpinner( uiOutput("umap2",
                                                             height = "300px") ))
                         ), # Ends Fluid row
                         
                         fluidRow(
                             titlePanel("Violin and Dot plots"),
                             column(width = 6,
                                    my_withSpinner( uiOutput("run_vln_plot") )),
                             column(width = 6,
                                    my_withSpinner( uiOutput("run_dot_plot") ))),
                         
                     ), # Ends conditional
                     
                     conditionalPanel (
                         condition = "input.run_feature_plot > 0",
                         
                         fluidRow(
                             titlePanel("Downloading additional plots"),
                             
                             p("In this section, it is possible to download the feature, violin, and dot plots of each selected gene."),
                             p("The files will be saved in the folder:", em("images/one_sample_plots_<current date>__<current time>")),
                             br(),
                             column(width = 2,
                                    radioButtons_down_add_plots("down_add_plots_tab1")),
                             conditionalPanel (
                                 condition = "input.down_add_plots_tab1 == 1",
                                 column(width = 2,
                                        my_withSpinner( uiOutput("select_genes_add_plot_to_down_ui") ),
                                 ),
                                 column(2,
                                        div(class = "option-group",
                                            p(strong("Select the options for the feature plots")) ),
                                        numericInput_plot_height("add_p_tab1_feat_height", value=10),
                                        numericInput_plot_width("add_p_tab1_feat_width", value=14),
                                        selectInput_plot_res("add_p_tab1_feat_res"),
                                        selectInput_plot_format("add_p_tab1_feat_format")),
                                 column(2,
                                        div(class = "option-group",
                                            p(strong("Select the options for the violin plots")) ),
                                        numericInput_plot_height("add_p_tab1_violin_height", value=10),
                                        numericInput_plot_width("add_p_tab1_violin_width", value=14),
                                        selectInput_plot_res("add_p_tab1_violin_res"),
                                        selectInput_plot_format("add_p_tab1_violin_format")),
                                 column(2,
                                        div(class = "option-group",
                                            p(strong("Select the options for the dot plots")) ),
                                        numericInput_plot_height("add_p_tab1_dot_height", value=10),
                                        numericInput_plot_width("add_p_tab1_dot_width", value=14),
                                        selectInput_plot_res("add_p_tab1_dot_res"),
                                        selectInput_plot_format("add_p_tab1_dot_format")),
                                 column(width = 3,
                                        actionButtonInput("start_down_add_plots_tab1",
                                                          HTML("Download additional plots")))))
                         
                     ),
                     
                     bookmarkButton(style = "position:absolute;right:2em; background-color:#BF3EFF; color:#FFFFFF;"),
                     
                     tags$hr(),
                     p(strong("SolusCell web app, version 1.0"), "- 2024.", align = "center")
                     # Ends page
            ),
            
            ####################################
            ### Integrative analysis - Tab 2 ###
            ####################################
            # 
            # tabPanel("Integration of multiple samples",
            #          br(),
            #          p(strong("Note that at any time, you can save a bookmark (purple button at the right bottom) that will save your parameter choices. Using the saved bookmark, it is possible to re-load all selected parameters and re-execute the analysis, to reproduce the results.")),
            #          br(),
            #          p("It is necessary to choose adequate parameters before loading the data for integration to avoid bias due to poor cells interfering in the anchoring step.", strong("Therefore, we recommend exploring each sample separately in the one sample tab and defining the best parameters before the integration.")),
            #          br(),
            #          p("It is necessary to have a configuration file in the csv format to load the samples for integration. Please visit",
            #            a(tags$a(href="https://asc-seurat.readthedocs.io/en/latest/loading_data_int.html#loading-the-data-and-integration-of-multiple-samples documentation", "Asc-Seurat's documentation", target="_blank")), "for instructions on how to generate this file."),
            #          br(),
            #          p("Alternatively, it is possible to load previously integrated data, saving some time by not running the integration step."),
            #          #p("After selecting the parameters, click on the blue button to load the data."),
            #          br(),
            #          fluidRow(
            #              column(width = 3,
            #                     div(class = "option-group",
            #                         radioButtons("integration_options",
            #                                      "Run a new integration analysis or read a previously saved file?",
            #                                      choices = list("Run a new analysis" = 0,
            #                                                     "Load an integrated, unclustered, dataset" = 2,
            #                                                     "Load an integrated, clustered, dataset" = 1),
            #                                      selected = c(0) )),
            #                     conditionalPanel (
            #                         condition = "input.integration_options == 0",
            #                         
            #                         div(class = "option-group",
            #                             fileInput("samples_list_integration",
            #                                       label = "Read the configuration file containing the samples' information",
            #                                       accept = c(
            #                                           'text/csv',
            #                                           'text/comma-separated-values',
            #                                           'text/tab-separated-values',
            #                                           'csv',
            #                                           'tsv') ) ),
            #                         
            #                         my_withSpinner( uiOutput('select_sample_tab2') )
            #                         #),
            #                     ),
            #              ),
            #              
            #              
            #              conditionalPanel(
            #                  condition = "input.integration_options != 0",
            #                  column(width = 3,
            #                         my_withSpinner( uiOutput("load_integrated_ui") ),
            #                         div(class = "option-group",
            #                             selectInput(
            #                                 inputId = "load_rds_int_normalization",
            #                                 label = "Inform the normalization method used to generate the integrated dataset",
            #                                 choices = list("",
            #                                                "LogNormalization" = 0,
            #                                                "SCtransform" = 1),
            #                                 selected = "",
            #                                 multiple = F))
            #                  ),
            #              ),
            #              conditionalPanel(
            #                  condition = "input.integration_options == 1",
            #                  column(width = 3,
            #                         actionButtonInput("load_rds_file2",
            #                                           HTML("Load the integrated, clustered, data")))
            #              ),
            #              conditionalPanel(
            #                  condition = "input.integration_options == 2",
            #                  column(width = 3,
            #                         actionButtonInput("load_rds_file3",
            #                                           HTML("Load the integrated, unclustered, data")))
            #              ),
            #              
            #              conditionalPanel (
            #                  condition = "input.integration_options == 0",
            #                  fluidRow(
            #                      column(width = 2,
            #                             textInput_regex("int_regex_mito"),
            #                             numericInput_var_genes("n_of_var_genes_integration",
            #                                                    label = "N of variable genes for integration"),
            #                             numericInput_n_of_PCs("n_of_PCs_integration",
            #                                                   "N of components for the integration"),
            #                             input_project_name("int_project_name")
            #                      ),
            #                      column(width = 3,
            #                             select_norm_methods("normaliz_method_tab2"),
            #                             
            #                             conditionalPanel(
            #                                 condition = "input.normaliz_method_tab2 == 0",
            #                                 numericInput_scale_factor("scale_factor_tab2"),
            #                                 radioInput_most_var_method("most_var_method_tab2"),
            #                             ),
            #                      ),
            #                      column(width = 3,
            #                             actionButtonInput("load_rds_file",
            #                                               HTML("Execute a new integration")))
            #                  ),
            #                  
            #                  
            #              ), # ends conditional
            #              
            #          ), # ends fluidRow
            #          
            #          conditionalPanel (
            #              condition = "input.integration_options != 1",
            #              
            #              conditionalPanel (
            #                  condition = "input.load_rds_file != 0",
            #                  
            #                  fluidRow(
            #                      titlePanel("Genes targeted as mitochondrial"),
            #                      p("The genes below were targeted by the parameter \"common identifier of mitochondrial genes\", set above."),
            #                      column(12,
            #                             my_withSpinner( verbatimTextOutput("target_genes_mitho_tab2") )
            #                      )
            #                  ), # Fluid row
            #              ),
            #              
            #              conditionalPanel (
            #                  condition = "input.integration_options != 1",
            #                  
            #                  fluidRow(
            #                      tags$hr(),
            #                      p("Since the integration is a time-consuming step, it is helpful to save the integrated data before running any additional analyses. The main advantage of saving the data is that if you need to restart the analysis (e. g., to change parameters), you can skip the integration step by loading the RDS file instead."),
            #                      br(),
            #                      p( strong("The rds files must be saved in the folder"), code("RDS_files.") ),
            #                      column(6,
            #                             downloadButton("download_int_data",
            #                                            "Download RDS object containing the integrated data.")
            #                      )),
            #              ),
            #              
            #              fluidRow(
            #                  titlePanel("Screening plot to define filtering parameters"),
            #                  br(),
            #                  p("Use this plot to define more restrictive parameters and exclude cells based on their number of expressed genes and the percentage of expressed genes from the mitochondria."),
            #                  #p("The parameters can be set on the right side of the plot."),
            #                  #br(),
            #                  p("After setting the parameters, click on \"Show/update plot of filtered data\" to visualize the data after filtering."),
            #                  
            #                  column(width = 10,
            #                         my_withSpinner( plotOutput("VlnPlot_tab2") )),
            #                  column(width = 2,
            #                         div(class = "option-group",
            #                             numericInput("min_count_tab2",
            #                                          label = "Keep only cells that expressed at least this number of genes",
            #                                          value = "")),
            #                         div(class = "option-group",
            #                             numericInput("max_count_tab2",
            #                                          label = "Exclude any cell that expressed more than this number of genes (i.e., possible doublets)",
            #                                          value = "")
            #                         ),
            #                         numericInput_max_mito_perc("max_mito_perc_tab2", value = "")),
            #                  column(width = 2,
            #                         actionButtonInput("run_vinplot_tab2",
            #                                           HTML("Show/update plot of filtered data")),
            #                         conditionalPanel(
            #                             condition = "input.run_vinplot_tab2 == 0",
            #                             actionButtonInput("run_pca_tab2_1",
            #                                               HTML("Run the PCA analysis"))
            #                         ),
            #                  ),
            #              ), # ends fluidRow
            #          ),
            #          
            #          # fluid row
            #          conditionalPanel(
            #              condition = "input.run_vinplot_tab2!= 0",
            #              
            #              fluidRow(
            #                  titlePanel("Screening plot showing the remaining cells after filtering"),
            #                  column(width = 10,
            #                         my_withSpinner(  plotOutput("VlnPlot_filt_tab2") )
            #                  ),
            #                  column(2,
            #                         numericInput_plot_height("p5_height", value=12),
            #                         numericInput_plot_width("p5_width", value=20),
            #                         selectInput_plot_res("p5_res"),
            #                         selectInput_plot_format("p5_format"),
            #                         downloadButton("p5_down", HTML("Download Plot")),
            #                         actionButtonInput("run_pca_tab2_2",
            #                                           HTML("Run the PCA analysis")))
            #                  
            #              ), # fluid row
            #          ), # Ends conditional
            #          
            #          conditionalPanel(
            #              condition = "input.run_pca_tab2_1 != 0 || input.run_pca_tab2_2 != 0",
            #              
            #              fluidRow(
            #                  titlePanel("Elbow plot showing the variation explained by each component"),
            #                  br(),
            #                  p("Choose the number of components that explain the most variation."),
            #                  br(),
            #                  column(width = 6,
            #                         
            #                         my_withSpinner( plotOutput("n_of_PCAs_tab2") )
            #                         
            #                  ),
            #                  column(2,
            #                         numericInput_plot_height("p6_height", value=10),
            #                         numericInput_plot_width("p6_width", value=18),
            #                         selectInput_plot_res("p6_res"),
            #                         selectInput_plot_format("p6_format"),
            #                         downloadButton("p6_down", HTML("Download Plot"))),
            #                  column(width = 3,
            #                         numericInput_n_of_PCs("n_of_PCs_tab2"))
            #              ) # Fluid row
            #          ), # Ends conditional
            #          
            #          ## Clustering tab2
            #          
            #          fluidRow(
            #              titlePanel("Clustering"),
            #              
            #              conditionalPanel (
            #                  condition = "input.integration_options != 1",
            #                  br(),
            #                  p("Be aware that this parameter is central in the cluster definition. It is recommended to try different values and define the most appropriate according to the expectations of the cell populations present in your sample."),
            #                  p("Quoting from", a(tags$a(href="https://satijalab.org/seurat/archive/v1.4/pbmc3k_tutorial.html", "Seurat's tutorial", target="_blank")), ":", em("\"We find that setting this parameter between 0.6-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets.\"")),
            #                  br(),
            #                  column(width = 3,
            #                         numericInput_resolution_clust("resolution_clust_tab2")),
            #                  column(width = 3,
            #                         actionButtonInput("run_clustering_tab2",
            #                                           HTML("Run the clustering analysis")))
            #              ), # Fluid row
            #              
            #          ),
            #          conditionalPanel(
            #              condition = "input.run_clustering_tab2 != 0 || input.load_rds_file2 != 0",
            #              
            #              fluidRow(
            #                  titlePanel("Clustering plots"),
            #                  column(width = 6,
            #                         
            #                         my_withSpinner( plotOutput("umap_tab2") )
            #                         
            #                  ),
            #                  column(width = 6,
            #                         
            #                         my_withSpinner( plotOutput("umap_three_samples_comb") )
            #                         
            #                  )
            #              ), # Fluid row
            #              fluidRow(
            #                  titlePanel("Clustering plots (UMAP) separated by sample"),
            #                  column(width = 12,
            #                         
            #                         my_withSpinner(  plotOutput("umap_three_samples") )
            #                  )
            #              ), # Fluid row
            #              
            #              fluidRow(
            #                  br(),
            #                  # download plot
            #                  column(width = 3,
            #                         div(class = "down-group",
            #                             radioButtons("p7_down_opt",
            #                                          "Select the plot to download",
            #                                          choices = list("UMAP" = "UMAP",
            #                                                         "UMAP colored by sample" = "UMAP1",
            #                                                         "UMAP split by sample"= "UMAP2"),
            #                                          selected = "UMAP"))),
            #                  column(2,
            #                         numericInput_plot_height("p7_height", value=10),
            #                         numericInput_plot_width("p7_width", value=12)),
            #                  column(2,
            #                         selectInput_plot_res("p7_res"),
            #                         selectInput_plot_format("p7_format")),
            #                  column(2,
            #                         downloadButton("p7_down", HTML("Download Plot"))),
            #              ),# Fluid row
            #              
            #              fluidRow(
            #                  titlePanel("Number of cells per cluster"),
            #                  p("The first row shows the cluster ID. The second row shows the number of cells per cluster."),
            #                  column(12,
            #                         my_withSpinner( verbatimTextOutput("cluster_size_tab2") ))
            #              ), # Fluid row
            #              
            #          ), # Ends conditional
            #          
            #          conditionalPanel (
            #              condition = "input.integration_options != 1",
            #              
            #              fluidRow(
            #                  titlePanel("Excluding or selecting clusters for reanalysis"),
            #                  br(),
            #                  p("Sometimes, it is helpful to exclude or select the clusters that are more of interest.", "After selecting or excluding the cells of interest, it is recommended to repeat the clustering step using only the subset."),
            #                  p("After selecting the clusters, click on the blue button (Reanalyze after selection/exclusion of clusters). Asc-Seurat will run the analyses of the new subset until the PCA step.", strong("Then, you will need to set the new number of components using the elbow plot (above) and click on the button \"Run the clustering analysis\" again.")),
            #                  
            #                  column(3,
            #                         radioButtons_filter_clusters("filter_clusters_tab2")
            #                  ),
            #                  
            #                  conditionalPanel(
            #                      condition = "input.filter_clusters_tab2 != 0",
            #                      
            #                      column(3,
            #                             radioButtons_filter_clusters_opt("filter_clusters_opt_tab2")
            #                      ),
            #                      
            #                      column( 3,
            #                              div(class = "option-group",
            #                                  my_withSpinner( uiOutput("cluster_list_tab2_ui") )
            #                              )
            #                      ),
            #                      
            #                      column( 3,
            #                              div(class = "option-group",
            #                                  actionButtonInput("rerun_after_filtering_tab2",
            #                                                    HTML("Reanalyze after selection/exclusion of clusters"))
            #                              )
            #                      ),
            #                      
            #                  ) # ends conditional
            #                  
            #              ), # ends fluid row
            #              
            #              fluidRow(
            #                  
            #                  br(),
            #                  h3("Saving the processed data for the trajectory inference analysis"),
            #                  
            #                  p("The button below allows for saving the processed data in a file that can be used for the pseudo-time  (trajectory inference) analysis."),
            #                  p("Asc-Seurat will save the most recently processed data."),
            #                  br(),
            #                  p(strong("Note:", "The processed data needs to be saved in the folder", code("RDS_files/"), "so it can load automatically in the tab \"Trajectory inference\"." )),
            #                  
            #                  column(6,
            #                         br(),
            #                         downloadButton("downloadRDS_tab2",
            #                                        "Download the processed data to use in the trajectory inference analysis",
            #                                        class = "down_butt")
            #                  )
            #              ),
            #          ),# ends conditional
            #          
            #          fluidRow(
            #              titlePanel("Identification of markers / Differential expression analysis"),
            #              column(3,
            #                     radioButtons_find_markers("find_markers_tab2")
            #              )
            #          ),
            #          
            #          conditionalPanel(
            #              condition = "input.find_markers_tab2 == 1",
            #              fluidRow(
            #                  br(),
            #                  p("Check Seurat's vignettes", a(tags$a(href="https://satijalab.org/seurat/v3.2/de_vignette.html", "here", target="_blank")), "and the function's",  a(tags$a(href="https://www.rdocumentation.org/packages/Seurat/versions/3.1.4/topics/FindAllMarkers", "manual", target="_blank")),  "to see the options for each parameter of this section. "),
            #                  br(),
            #                  column(3,
            #                         div(class = "option-group",
            #                             radioButtons("find_markers_or_DE_tab2",
            #                                          "Detect conserved markers or D.E. genes",
            #                                          choices = c("Conserved markers" = 0,
            #                                                      "D.E. genes between samples" = 1),
            #                                          selected = 0
            #                             ))),
            #                  conditionalPanel(
            #                      condition = "input.find_markers_or_DE_tab2 == 0",
            #                      column(3,
            #                             div(class = "option-group",
            #                                 radioButtons("find_markers_tab2_opt",
            #                                              "Select the analysis to perform?",
            #                                              choices = list("Identify conserved markers for all clusters" = 0,
            #                                                             "Identify conserved markers for one specific cluster" = 1,
            #                                                             "Identify conserved markers for one specific cluster in comparison with another cluster(s)" = 2),
            #                                              selected = c(0)),
            #                                 
            #                             )
            #                      ),
            #                      
            #                      conditionalPanel(
            #                          condition = "input.find_markers_tab2_opt == 1",
            #                          column(3,
            #                                 div(class = "option-group",
            #                                     my_withSpinner(  uiOutput("find_markers_clust_id_tab2_ui") )
            #                                 )),
            #                      ),   # ends conditional
            #                      conditionalPanel(
            #                          condition = "input.find_markers_tab2_opt == 2",
            #                          column(3,
            #                                 div(class = "option-group",
            #                                     
            #                                     my_withSpinner( uiOutput("find_markers_clust_ID1_tab2_ui") ),
            #                                     my_withSpinner( uiOutput("find_markers_clust_ID2_tab2_ui") )
            #                                     
            #                                 )),
            #                      ),    # ends conditional
            #                  ),
            #                  conditionalPanel(
            #                      condition = "input.find_markers_or_DE_tab2 == 1",
            #                      column(3,
            #                             div(class = "option-group",
            #                                 my_withSpinner( uiOutput("find_markers_or_DE_tab2_cluster_ui") ),
            #                                 my_withSpinner( uiOutput("find_markers_or_DE_tab2_treat1_ui") ),
            #                                 my_withSpinner( uiOutput("find_markers_or_DE_tab2_treat2_ui") ),
            #                                 pickerInput_find_markers_test("find_markers_tab2_test.use"),
            #                                 
            #                                 numericInput("find_markers_or_DE_tab2_pvalue",
            #                                              label = "Select the (adjusted) p-value threshold.",
            #                                              value = 0.05),
            #                                 
            #                             ))),
            #                  
            #                  column(3,
            #                         actionButtonInput("run_ident_markers_tab2",
            #                                           HTML("Search for markers/D.E. genes"))),
            #                  
            #              ), # ends conditional
            #          ), # ends fluidrow
            #          
            #          conditionalPanel(
            #              condition = "input.find_markers_tab2 == 1",
            #              fluidRow(
            #                  titlePanel("List of markers or differentially expressed genes"),
            #                  conditionalPanel(
            #                      condition = "input.run_ident_markers_tab2 != 0",
            #                      column(12,
            #                             my_withSpinner( reactableOutput("markers_tab2_react")) ),
            #                  ),
            #                  br(),
            #                  column(4,
            #                         downloadButton("download_markers_tab2",
            #                                        "Download the list of markers or D.E. genes."))
            #              ), # ends fluidrow
            #          ), # ends conditional
            #          
            #          
            #          fluidRow(
            #              titlePanel("Expression of markers"),
            #              br(),
            #              p("In this section, users can visualize the gene expression of selected genes (e.g., tissue markers). Start by loading a", strong("csv"), "file with at least two columns. The first column must be the gene ID, and the second is a grouping variable (e.g., the tissue name). A third column can be used to store the common name of the gene but is optional."),
            #              br(),
            #              
            #              column(width = 3,
            #                     fileInput_markers_list("markers_list_tab2"),
            #                     div(class = "option-group",
            #                         selectInput("markers_list_header_opt_tab2",
            #                                     "Does your file have a header?",
            #                                     choices = c("", "Yes", "No"),
            #                                     multiple = FALSE,
            #                                     selected = NULL)),
            #                     div(class = "option-group",
            #                         actionButtonInput("load_markers_tab2",
            #                                           HTML("Load markers")))),
            #              column(width = 2,
            #                     conditionalPanel(
            #                         condition = "input.load_markers_tab2 > 0",
            #                         my_withSpinner( uiOutput('marker_group_selec_tab2') ),
            #                         define_if_use_all_genes_or_select("filter_genes_q_tab2"),
            #                         
            #                         conditionalPanel(
            #                             condition = "input.filter_genes_q_tab2 == 0",
            #                             define_what_id_to_use("genes_ids_tab2")
            #                         ))),
            #              conditionalPanel(
            #                  condition = "input.filter_genes_q_tab2 == 0",
            #                  
            #                  column(width = 3,
            #                         my_withSpinner( uiOutput('marker_genes_selec_tab2') )
            #                  )),
            #              column(width = 2,
            #                     conditionalPanel(
            #                         condition = "input.load_markers_tab2 > 0",
            #                         radioButtons_slot_selection_heatmap("slot_selection_heatmap_tab2"))),
            #              column(width = 2,
            #                     conditionalPanel(
            #                         condition = "input.load_markers_tab2 > 0",
            #                         actionButtonInput("run_heatmap_tab2",
            #                                           HTML("Show heatmap"))))
            #              
            #          ), # Fluid row
            #          
            #          #heat map
            #          conditionalPanel(
            #              
            #              condition = "input.run_heatmap_tab2 > 0",
            #              fluidRow(
            #                  br(),
            #                  
            #                  titlePanel("Heatmap"),
            #                  br(),
            #                  strong("Note that the scale of colors of the heatmap is adjusted based on the expression of the selected genes."),
            #                  br(),
            #                  p("For now, the heat map shows the average expression of all samples together. It is only helpful to identify if the markers or cell types make sense with the number of clusters."),
            #                  p("To visualize the expression in the clusters by sample, use the feature plots."),
            #                  
            #                  column(width = 10,
            #                         my_withSpinner(  uiOutput("heat_map_ui_tab2") )),
            #                  # download plot
            #                  column(2,
            #                         numericInput_plot_height("p8_height", value=15),
            #                         numericInput_plot_width("p8_width", value=20),
            #                         selectInput_plot_res("p8_res"),
            #                         selectInput_plot_format("p8_format"),
            #                         downloadButton("p8_down", HTML("Download Plot")))
            #              ),
            #              
            #              fluidRow(
            #                  titlePanel("Visualization of gene expression of each cell and additional plots"),
            #                  br(),
            #                  
            #                  p("Only genes selected for the heatmap can be selected for the additional plots."),
            #                  br(),
            #                  column(width = 3,
            #                         my_withSpinner( uiOutput("marker_to_feature_plot_tab2") )),
            #                  column(width = 2,
            #                         radioButtons_slot_selection_feature_plot("slot_selection_feature_plot_tab2")),
            #                  column(width = 3,
            #                         actionButtonInput("run_feature_plot_tab2",
            #                                           HTML("Show the expression of genes at the cell level")))
            #              )
            #          ),# Ends conditional
            #          
            #          conditionalPanel(
            #              condition = "input.run_feature_plot_tab2 > 0",
            #              fluidRow(
            #                  titlePanel("Feature plots"),
            #                  column(width = 7,
            #                         my_withSpinner( uiOutput("feature_plot_tab2") )),
            #                  column(width = 4,
            #                         my_withSpinner( uiOutput("umap2_tab2", height = "300px") ))
            #              ), # Ends Fluid row
            #              
            #              fluidRow(
            #                  titlePanel("Violin and Dot plots"),
            #                  column(width = 12,
            #                         my_withSpinner( uiOutput("run_vln_plot_tab2") ))#,
            #                  
            #              ),
            #              fluidRow(
            #                  titlePanel("Downloading additional plots"),
            #                  
            #                  p("In this section, it is possible to download the feature, violin, and dot plots of each selected gene."),
            #                  p("The files will be saved in the folder:", em("images/one_sample_plots_<current date>__<current time>")),
            #                  br(),
            #                  column(width = 2,
            #                         radioButtons_down_add_plots("down_add_plots_tab2")),
            #                  conditionalPanel (
            #                      condition = "input.down_add_plots_tab2 == 1",
            #                      column(width = 2,
            #                             my_withSpinner( uiOutput("select_genes_add_plot_to_down_tab2_ui") ),
            #                      ),
            #                      
            #                      column(2,
            #                             div(class = "option-group",
            #                                 p(strong("Select the options for the feature plots")) ),
            #                             numericInput_plot_height("add_p_tab2_feat_height", value=10),
            #                             numericInput_plot_width("add_p_tab2_feat_width", value=14),
            #                             selectInput_plot_res("add_p_tab2_feat_res"),
            #                             selectInput_plot_format("add_p_tab2_feat_format")),
            #                      column(2,
            #                             div(class = "option-group",
            #                                 p(strong("Select the options for the violin plots")) ),
            #                             numericInput_plot_height("add_p_tab2_violin_height", value=10),
            #                             numericInput_plot_width("add_p_tab2_violin_width", value=14),
            #                             selectInput_plot_res("add_p_tab2_violin_res"),
            #                             selectInput_plot_format("add_p_tab2_violin_format")),
            #                      column(width = 3,
            #                             actionButtonInput("start_down_add_plots_tab2",
            #                                               HTML("Download additional plots")))
            #                  )
            #              ), # end fluid row
            #          ), # Ends conditional
            #          bookmarkButton(style = "position:absolute;right:2em; background-color:#BF3EFF; color:#FFFFFF;"),
            #          
            #          tags$hr(),
            #          p(strong("SolusCell web app, version 1.0"), "- 2024.", align = "center")
            # ), # ends tab
            
            ##########################
            ##### Advanced plots #####
            ##########################
            tabPanel("Advanced plots",
                     
                     #   h2("Advanced Plots"),
                     
                     fluidPage(
                         fluidRow(
                             stacked_violin_UI("stacked1"),
                         )
                     )
                     
            ),
            
        ), ## Closing the ui
        width = 12)
}
