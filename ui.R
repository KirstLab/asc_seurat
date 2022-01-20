# Asc-Seurat
# Version 2.2
set.seed(1407)

suppressMessages( require("shiny") )
suppressMessages( require("shinyWidgets") )
suppressMessages( require("DT") )
suppressMessages( require("reactable") )
suppressMessages( require("rclipboard") )
suppressMessages( require("shinycssloaders") )
suppressMessages( require("shinyFeedback") )

if (dir.exists('/app/user_work')) {
    source("/app/R/ui_functions.R")
    source("/app/R/improved_dot_and_violin_plots.R")
} else {
    source("R/ui_functions.R")
    source("R/improved_dot_and_violin_plots.R")
    
}
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
            tags$link(rel = "stylesheet", type = "text/css", href = "asc_seurat.css")
        ),
        
        tags$head(
            tags$style(HTML("
            p {
                font-size:16px;
            }

            .down-group {
                        border: 1px solid #ccc;
                        border-radius: 6px;
                        padding: 0px 5px;
                        margin: 5px -10px;
                        background-color: #eadaf7;
            } "))),
        
        rclipboardSetup(),
        
        # Sets the background color of the application
        setBackgroundColor(
            color = c("#F7FBFF", "#d6c5e3"),
            gradient = "radial",
            direction = c("top", "left")),
        
        ########################
        # Asc-Seurat interface #
        ########################
        
        titlePanel(title = "Asc-Seurat - Analytical single-cell Seurat-based web application"),
        br(),
        
        tabsetPanel(
            
            ######################
            ## Introduction tab ##
            ######################
            
            tabPanel("Introduction",
                     
                     h2("Introduction to Asc-Seurat"),
                     br(),
                     p("Asc-Seurat, pronounced \"ask Seurat\", is based on the popular R package \"Seurat\", from the Satija Lab. It includes many, but not all, features of the Seurat package."),
                     p("It also takes advantage of \"Dynverse\", a collection of R packages that allows the execution of multiple trajectory inference models."),
                     p("Finally, Asc_Seurat uses BioMart, through the biomaRt R package, to promote the functional annotation of genes from many species."),
                     br(),
                     
                     p(
                         "For more information, check Seurat's manual and vignettes",
                         a(tags$a(href="https://satijalab.org/seurat/", "here,", target="_blank")),
                         "and their publications",
                         a(tags$a(href="https://satijalab.org/publications/", "here.", target="_blank")),
                         "Also, check dynverse's documentation",
                         a(tags$a(href="https://dynverse.org/", "here", target="_blank")),
                         "and its publication",
                         a(tags$a(href="https://doi.org/10.1038/s41587-019-0071-9", "here;", target="_blank")),
                         "and BioMart's documentation",
                         a(tags$a(href="https://www.ensembl.org/biomart/martview/dc4f7144c82d3d4b4c1bd8f27ca07b6c", "here,", target="_blank")),
                         "and biomaRt's vignettes",
                         a(tags$a(href="https://bioconductor.org/packages/release/bioc/html/biomaRt.html", "here.", target="_blank"))
                     ),
                     
                     tags$hr(),
                     h2("Tutorial and documentation"),
                     br(),
                     p("A step-by-step introduction of Asc-Seurat's functionalities is available at ", a(tags$a(href="https://asc-seurat.readthedocs.io/en/latest/index.html", "https://asc-seurat.readthedocs.io.", target="_blank"))),
                     p("For questions or issues related to Asc-Seurat's functionalities, please visit", a(tags$a(href="https://github.com/KirstLab/asc_seurat", "our GitHub.", target="_blank"))),
                     
                     br(),
                     
                     p("Note that both Seurat and dynverse teams are also present on GitHub. For questions or issues related to Seurat's or dynverse's functionalities, please visit", a(tags$a(href="https://github.com/satijalab/seurat/issues", "Seurat's GitHub page", target="_blank")), "or", a(tags$a(href="https://github.com/dynverse/dyno/issues", "Dynverse's GitHub page.", target="_blank"))),
                     tags$hr(),
                     h2("Funding"),
                     br(),
                     p("Asc-Seurat was developed in the context of a project supported by the US Department of Energy, Office of Science Biological and Environmental Research [DE-SC0018247]."),
                     br(),
                     tags$hr(),
                     h2("Acknowledgments"),
                     ## Adds the images
                     a(
                         img(src = 'University_of_Florida.png',
                             width = "30%",
                             align = "center",
                             style='border-right: 40px solid transparent'),
                         img(src = 'Universidade_de_Brasilia.png',
                             width = "20%",
                             align = "center",
                             style='border-right: 40px solid transparent'),
                     ),
                     tags$hr(),
                     p(strong("Asc-Seurat, version 2.1"), "- Released on May 26th, 2021.", align = "center")
            ),
            
            ######################################
            ###           One sample           ###
            ######################################
            
            tabPanel("One sample",
                     fluidRow(
                         br(),
                         p(strong("Note that at any time, users can save a bookmark (purple button at the right bottom) that will save your parameter choices. Using the saved bookmark, it is possible to re-load all selected parameters and re-execute the analysis, to reproduce the results.")),
                         br(),
                         p("Choose the sample to be analyzed and the initial requirements to load the data. Note that cells that do not match the parameters will not be load."),
                         p("These parameters are used to exclude low-quality cells and allow the data to load quickly. Users can add more restrictive parameters after visualizing the distributions in the next section."),
                         p(strong("After selecting the parameters, click on the blue button to load the data.")),
                         br(),
                         
                         # column(width = 3,
                         #        my_withSpinner(
                         #            uiOutput('select_sample_tab1'))
                         # ),
                         
                         column(width = 3,
                                div(class = "option-group",
                                    radioButtons("sample_tab1_options",
                                                 "Run a new analysis or read a previously saved file?",
                                                 choices = list("Run a new analysis" = 0,
                                                                "Load file" = 1),
                                                 selected = c(0) )),
                                
                                conditionalPanel (
                                    condition = "input.sample_tab1_options == 0",
                                    
                                    #column(width = 3,
                                    my_withSpinner(
                                        uiOutput('select_sample_tab1'))
                                    
                                    #),
                                ),
                         ),
                         
                         conditionalPanel (
                             condition = "input.sample_tab1_options == 1",
                             
                             column(width = 3,
                                    my_withSpinner( uiOutput("select_sample_tab1_rds_ui") ),
                                    # div(class = "option-group",
                                    #     selectInput(
                                    #         inputId = "select_sample_tab1_rds_normalization",
                                    #         label = "Inform the normalization method used to generate the dataset",
                                    #         choices = list("",
                                    #                        "LogNormalization" = 0,
                                    #                        "SCtransform" = 1),
                                    #         selected = "",
                                    #         multiple = F))
                                    
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
                                    # div(class = "option-group",
                                    #     textInput("proj_name",
                                    #               label = "Project name",
                                    #               value = ""))
                             ),
                             
                             column(width = 2,
                                    div(class = "option-group",
                                        numericInput("min_cells",
                                                     label = "Min. number of cells expressing a gene for the gene to be included",
                                                     value = 3))
                             ),
                             
                             column(width = 2,
                                    div(class = "option-group",
                                        numericInput("min_features",
                                                     label = "Min. number of genes a cell must express to be included",
                                                     value = 200))
                             ),
                             
                             column( width = 3,
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
                         
                         fluidRow(
                             titlePanel("Screening plot to define filtering parameters"),
                             br(),
                             p("Use this plot to define more restrictive parameters and exclude cells based on their number of expressed genes and the percentage of reads that map to the mitochondrial genome."),
                             p("The parameters can be set on the right side of the plot and must be set using a higher value than the ones above. Otherwise, they will have no effect since the cells were already excluded."),
                             br(),
                             p("After setting the parameters, click on \"Show plot of filtered data\" to visualize the data after filtering."),
                             br(),
                             br(),
                             column(width = 10,
                                    my_withSpinner( plotOutput("VlnPlot") )
                             ),
                             column(width = 2,
                                    div(class = "option-group",
                                        numericInput("min_count",
                                                     label = "Keep only cells that expressed at least this number of genes",
                                                     value = 200),
                                        numericInput("max_count",
                                                     label = "Exclude any cell that expressed more than this number of genes (i.e. possible doublets)",
                                                     value = 2200)),
                                    
                                    #  my_withSpinner( uiOutput("min_count_ui")),
                                    #  my_withSpinner( uiOutput("max_count_ui"))),
                                    
                                    #lognorm_UI("lognorm_sing_sample"),
                                    numericInput_max_mito_perc("max_mito_perc"),
                                    actionButtonInput("run_vinplot",
                                                      HTML("Show plot of filtered data"))
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
                             p("Note that only the most variable genes are used in the dimension reduction step (PCA)."),
                             br(),
                             column(width = 3,
                                    select_norm_methods("normaliz_method")
                                    
                                    # div(class = "option-group",
                                    #     radioButtons("normaliz_method",
                                    #                  "Select the normalization method",
                                    #                  choices = list("LogNormalize" = 0,
                                    #                                 "SCTransform"  = 1),
                                    #                  selected = 0))
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
                                        
                                        #   download_ElbowPlot_UI("down_elbow1")
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
                     
                     fluidRow(
                         titlePanel("Clustering"),
                         
                         conditionalPanel (
                             condition = "input.sample_tab1_options == 0",
                             
                             br(),
                             p("Be aware that this parameter is central in the cluster definition. It is recommended to try different values and define the most appropriate according to the expectations of the cell populations present in the sample."),
                             p("Quoting from", a(tags$a(href="https://satijalab.org/seurat/archive/v1.4/pbmc3k_tutorial.html", "Seurat's tutorial", target="_blank")), em("\"We find that setting this parameter between 0.6-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets.\"")),
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
                     
                     
                     conditionalPanel(
                         condition = "input.run_clustering!= 0 || input.load_10X_rds!= 0",
                         
                         fluidRow(
                             titlePanel("Clustering plots"),
                             column(width = 6,
                                    my_withSpinner( plotOutput("umap") )
                             ),
                             column(width = 6,
                                    my_withSpinner( plotOutput("tSNE") )
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
                             p("First row shows the cluster ID. The second row shows the number of cells per cluster."),
                             column(12,
                                    my_withSpinner( verbatimTextOutput("cluster_size") )
                             )
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
                             p("Asc-Seurat will save only the most recently processed data."),
                             br(),
                             p(strong("Note:", "The processed data needs to be saved in the folder", code("RDS_files/"), "so it can load automatically in the tab \"Trajectory inference.\"" )),
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
                                                     label = "Select the minimal percentage of cells expressing a gene for this gene to be tested",
                                                     value = 0.1),
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
                         p("In this section, users can visualize the gene expression of selected genes (e.g., tissue markers). Start by loading a", strong("csv"), "file with at least two columns. The first column must be the gene ID, and the second is a grouping variable (e.g., the tissue name). A third column can be used to store the common name of the gene but is optional."),
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
                     p(strong("Asc-Seurat, version 2.1"), "- Released on May 26th, 2021.", align = "center")
                     # Ends page
            ),
            
            ####################################
            ### Integrative analysis - Tab 2 ###
            ####################################
            
            tabPanel("Integration of multiple samples",
                     br(),
                     p(strong("Note that at any time, you can save a bookmark (purple button at the right bottom) that will save your parameter choices. Using the saved bookmark, it is possible to re-load all selected parameters and re-execute the analysis, to reproduce the results.")),
                     br(),
                     p("It is necessary to choose adequate parameters before loading the data for integration to avoid bias due to poor cells interfering in the anchoring step.", strong("Therefore, we recommend exploring each sample separately in the one sample tab and defining the best parameters before the integration.")),
                     br(),
                     p("To load the sample, it is necessary to have a configuration file in the csv format. In this file, each row contains the information of one sample, and each column is one parameter for filtering the sample before integration.", strong("The file must have six columns and a header, as described below:")),
                     tags$ul(
                         tags$li(strong("Folder name:"), "Name of the folder where the raw data is (must be a subfolder of data/)"),
                         tags$li(strong("Sample name:"), "If the samples are part of a timepoint series, and you want to use this information in the trajectory analysis, the names need to be numeric (0 instead of 0h, for example. If you have replicates that should be combined, use the same name for all replicates)."),
                         tags$li(strong("Minimum number of cells expressing a gene")),
                         tags$li(strong("Minimum number of genes a cell must express to be included")),
                         tags$li(strong("Maximum number of genes a cell can express and still be included")),
                         tags$li(strong("Maximum percentage of genes belonging to the mitochondrial genome"))
                     ),
                     br(),
                     p("Alternatively, it is possible to load previously integrated data, saving some time by not running the integration step. For that, save an RDS file containing the integrated data."),
                     p("After selecting the parameters, click on the blue button to load the data."),
                     br(),
                     fluidRow(
                         column(width = 3,
                                div(class = "option-group",
                                    radioButtons("integration_options",
                                                 "Run a new integration analysis or read a previously saved file?",
                                                 choices = list("Run a new analysis" = 0,
                                                                "Load file" = 1),
                                                 selected = c(0) )),
                                conditionalPanel (
                                    condition = "input.integration_options == 0",
                                    
                                    div(class = "option-group",
                                        fileInput("samples_list_integration",
                                                  label = "Read the configuration file containing the samples' information",
                                                  accept = c(
                                                      'text/csv',
                                                      'text/comma-separated-values',
                                                      'text/tab-separated-values',
                                                      'csv',
                                                      'tsv') ) ),
                                    
                                    my_withSpinner( uiOutput('select_sample_tab2') )
                                    #),
                                ),
                         ),
                         
                         
                         conditionalPanel(
                             condition = "input.integration_options == 1",
                             column(width = 3,
                                    my_withSpinner( uiOutput("load_integrated_ui") ),
                                    # div(class = "option-group",
                                    #     selectInput(
                                    #         inputId = "load_rds_int_normalization",
                                    #         label = "Inform the normalization method used to generate the integrated dataset",
                                    #         choices = list("",
                                    #                        "LogNormalization" = 0,
                                    #                        "SCtransform" = 1),
                                    #         selected = "",
                                    #         multiple = F))
                             ),
                             
                             column(width = 3,
                                    actionButtonInput("load_rds_file2",
                                                      HTML("Load the integrated data")))
                         ),
                         # conditionalPanel (
                         #     condition = "input.integration_options == 0",
                         #
                         #     column(width = 3,
                         #            my_withSpinner( uiOutput('select_sample_tab2') )
                         #     ),
                         # ),
                         
                         conditionalPanel (
                             condition = "input.integration_options == 0",
                             fluidRow(
                                 column(width = 2,
                                        textInput_regex("int_regex_mito"),
                                        numericInput_var_genes("n_of_var_genes_integration",
                                                               label = "N of variable genes for integration"),
                                        numericInput_n_of_PCs("n_of_PCs_integration",
                                                              "N of components for the integration"),
                                        input_project_name("int_project_name")
                                        # div(class = "option-group",
                                        #     textInput("int_project_name",
                                        #               label = "Type the project name",
                                        #               value = "")
                                        # )
                                 ),
                                 column(width = 3,
                                        select_norm_methods("normaliz_method_tab2"),
                                        
                                        conditionalPanel(
                                            condition = "input.normaliz_method_tab2 == 0",
                                            # column(width = 3,
                                            numericInput_scale_factor("scale_factor_tab2"),
                                            radioInput_most_var_method("most_var_method_tab2"),
                                            #numericInput_var_genes("n_of_var_genes_tab2"),
                                            
                                            #            actionButtonInput("run_pca_tab2",
                                            #                              HTML("Run the PCA analysis")))
                                        ),
                                 ),
                                 column(width = 3,
                                        actionButtonInput("load_rds_file",
                                                          HTML("Execute a new integration")))
                             ),
                             
                             
                         ), # ends conditional
                         
                         
                         
                         # column(width = 3,
                         #        actionButtonInput("load_rds_file",
                         #                          HTML("Load the integrated data <br> or <br> execute a new integration"))),
                         
                     ), # ends fluidRow
                     
                     conditionalPanel (
                         condition = "input.integration_options == 0",
                         
                         fluidRow(
                             tags$hr(),
                             p("Since the integration is a time-consuming step, it is helpful to save the integrated data before running any additional analyses. The main advantage of saving the data is that if you need to restart the analysis (e. g., to change parameters), you can skip the integration step by loading the RDS file instead."),
                             br(),
                             p(strong("The rds files must be saved in the folder"), code("RDS_files."), "Otherwise, it will not be possible to read them with the \"Load file\" option."),
                             column(6,
                                    downloadButton("download_int_data",
                                                   "Download RDS object containing the integrated data.")
                             )),
                         
                         
                         fluidRow(
                             titlePanel("Screening plot to define the filter parameters to exclude cells based on counts and % of mitochondrial"),
                             br(),
                             p("Use this plot to define more restrictive parameters and exclude cells based on their number of expressed genes and the percentage of expressed genes from the mitochondria."),
                             p("The parameters can be set on the right side of the plot."),
                             br(),
                             p("After setting the parameters, click on \"Show plot of filtered data\" to visualize the data after filtering."),
                             
                             column(width = 10,
                                    my_withSpinner( plotOutput("VlnPlot_tab2") )),
                             column(width = 2,
                                    div(class = "option-group",
                                        numericInput("min_count_tab2",
                                                     label = "Keep only cells that expressed at least this number of genes",
                                                     value = "")),
                                    div(class = "option-group",
                                        numericInput("max_count_tab2",
                                                     label = "Exclude any cell that expressed more than this number of genes (i.e., possible doublets)",
                                                     value = "")
                                    ),
                                    numericInput_max_mito_perc("max_mito_perc_tab2", value = "")),
                             column(width = 2,
                                    actionButtonInput("run_vinplot_tab2",
                                                      HTML("Show plot of filtered data")),
                                    conditionalPanel(
                                        condition = "input.run_vinplot_tab2 == 0",
                                        actionButtonInput("run_pca_tab2_1",
                                                          HTML("Run the PCA analysis"))
                                    ),
                             ),
                         ), # ends fluidRow
                     ),
                     
                     # fluid row
                     conditionalPanel(
                         condition = "input.run_vinplot_tab2!= 0",
                         
                         fluidRow(
                             titlePanel("Screening plot showing the remaining cells after filtering"),
                             column(width = 10,
                                    my_withSpinner(  plotOutput("VlnPlot_filt_tab2") )
                             ),
                             column(2,
                                    numericInput_plot_height("p5_height", value=12),
                                    numericInput_plot_width("p5_width", value=20),
                                    selectInput_plot_res("p5_res"),
                                    selectInput_plot_format("p5_format"),
                                    downloadButton("p5_down", HTML("Download Plot")),
                                    actionButtonInput("run_pca_tab2_2",
                                                      HTML("Run the PCA analysis")))
                             
                         ), # fluid row
                     ), # Ends conditional
                     
                     
                     #            select_norm_methods("normaliz_method_tab2")
                     #
                     #            # div(class = "option-group",
                     #            #     radioButtons("normaliz_method_tab2",
                     #            #                  "Select the normalization method",
                     #            #                  choices = list(#"SCTransform (recommended)" = "SCTransform",
                     #            #                      "LogNormalize" = "LogNormalize"),
                     #            #                  selected = c("LogNormalize")))
                     #
                     #     ),
                     #     conditionalPanel(
                     #         condition = "input.normaliz_method_tab2 == 0",
                     #         column(width = 4,
                     #                numericInput_scale_factor("scale_factor_tab2")),
                     #         column(width = 4,
                     #                radioInput_most_var_method("most_var_method_tab2"),
                     #                numericInput_var_genes("n_of_var_genes_tab2"),
                     #                actionButtonInput("run_pca_tab2",
                     #                                  HTML("Run the PCA analysis")))
                     #     ),
                     # ), # Fluid row,
                     
                     conditionalPanel(
                         condition = "input.run_pca_tab2_1 != 0 || input.run_pca_tab2_2 != 0",
                         
                         fluidRow(
                             titlePanel("Elbow plot showing the variation explained by each component"),
                             br(),
                             p("Choose the number of components that explain the most variation."),
                             br(),
                             column(width = 6,
                                    
                                    my_withSpinner( plotOutput("n_of_PCAs_tab2") )
                                    
                             ),
                             column(2,
                                    numericInput_plot_height("p6_height", value=10),
                                    numericInput_plot_width("p6_width", value=18),
                                    selectInput_plot_res("p6_res"),
                                    selectInput_plot_format("p6_format"),
                                    downloadButton("p6_down", HTML("Download Plot"))),
                             column(width = 3,
                                    numericInput_n_of_PCs("n_of_PCs_tab2"))
                         ) # Fluid row
                     ), # Ends conditional
                     
                     ## Clustering tab2
                     
                     fluidRow(
                         titlePanel("Clustering"),
                         
                         conditionalPanel (
                             condition = "input.integration_options == 0",
                             br(),
                             p("Be aware that this parameter is central in the cluster definition. It is recommended to try different values and define the most appropriate according to the expectations of the cell populations present in your sample."),
                             p("Quoting from", a(tags$a(href="https://satijalab.org/seurat/archive/v1.4/pbmc3k_tutorial.html", "Seurat's tutorial", target="_blank")), ":", em("\"We find that setting this parameter between 0.6-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets.\"")),
                             br(),
                             column(width = 3,
                                    numericInput_resolution_clust("resolution_clust_tab2")),
                             column(width = 3,
                                    actionButtonInput("run_clustering_tab2",
                                                      HTML("Run the clustering analysis")))
                         ), # Fluid row
                         
                     ),
                     conditionalPanel(
                         condition = "input.run_clustering_tab2 != 0 || input.load_rds_file2 != 0",
                         
                         fluidRow(
                             titlePanel("Clustering plots"),
                             column(width = 6,
                                    
                                    my_withSpinner( plotOutput("umap_tab2") )
                                    
                             ),
                             column(width = 6,
                                    
                                    my_withSpinner( plotOutput("umap_three_samples_comb") )
                                    
                             )
                         ), # Fluid row
                         fluidRow(
                             titlePanel("Clustering plots (UMAP) separated by sample"),
                             column(width = 12,
                                    
                                    my_withSpinner(  plotOutput("umap_three_samples") )
                             )
                         ), # Fluid row
                         
                         fluidRow(
                             br(),
                             # download plot
                             column(width = 3,
                                    div(class = "down-group",
                                        radioButtons("p7_down_opt",
                                                     "Select the plot to download",
                                                     choices = list("UMAP" = "UMAP",
                                                                    "UMAP colored by sample" = "UMAP1",
                                                                    "UMAP split by sample"= "UMAP2"),
                                                     selected = "UMAP"))),
                             column(2,
                                    numericInput_plot_height("p7_height", value=10),
                                    numericInput_plot_width("p7_width", value=12)),
                             column(2,
                                    selectInput_plot_res("p7_res"),
                                    selectInput_plot_format("p7_format")),
                             column(2,
                                    downloadButton("p7_down", HTML("Download Plot"))),
                         ),# Fluid row
                         
                         fluidRow(
                             titlePanel("Number of cells per cluster"),
                             p("The first row shows the cluster ID. The second row shows the number of cells per cluster."),
                             column(12,
                                    my_withSpinner( verbatimTextOutput("cluster_size_tab2") ))
                         ), # Fluid row
                         
                     ), # Ends conditional
                     
                     conditionalPanel (
                         condition = "input.integration_options == 0",
                         
                         fluidRow(
                             titlePanel("Excluding or selecting clusters for reanalysis"),
                             br(),
                             p("Sometimes, it is helpful to exclude or select the clusters that are more of interest.", "After selecting or excluding the cells of interest, it is recommended to repeat the clustering step using only the subset."),
                             p("After selecting the clusters, click on the blue button (Reanalyze after selection/exclusion of clusters). Asc-Seurat will run the analyses of the new subset until the PCA step.", strong("Then, you will need to set the new number of components using the elbow plot (above) and click on the button \"Run the clustering analysis\" again.")),
                             
                             column(3,
                                    radioButtons_filter_clusters("filter_clusters_tab2")
                             ),
                             
                             conditionalPanel(
                                 condition = "input.filter_clusters_tab2 != 0",
                                 
                                 column(3,
                                        radioButtons_filter_clusters_opt("filter_clusters_opt_tab2")
                                 ),
                                 
                                 column( 3,
                                         div(class = "option-group",
                                             my_withSpinner( uiOutput("cluster_list_tab2_ui") )
                                         )
                                 ),
                                 
                                 column( 3,
                                         div(class = "option-group",
                                             actionButtonInput("rerun_after_filtering_tab2",
                                                               HTML("Reanalyze after selection/exclusion of clusters"))
                                         )
                                 ),
                                 
                             ) # ends conditional
                             
                         ), # ends fluid row
                         
                         fluidRow(
                             
                             br(),
                             h3("Saving the processed data for the trajectory inference analysis"),
                             
                             p("The button below allows for saving the processed data in a file that can be used for the pseudo-time  (trajectory inference) analysis."),
                             p("Asc-Seurat will save the most recently processed data."),
                             br(),
                             p(strong("Note:", "The processed data needs to be saved in the folder", code("RDS_files/"), "so it can load automatically in the tab \"Trajectory inference\"." )),
                             
                             column(6,
                                    br(),
                                    downloadButton("downloadRDS_tab2",
                                                   "Download the processed data to use in the trajectory inference analysis",
                                                   class = "down_butt")
                             )
                         ),
                     ),# ends conditional
                     
                     fluidRow(
                         titlePanel("Identification of markers / Differential expression analysis"),
                         column(3,
                                radioButtons_find_markers("find_markers_tab2")
                         )
                     ),
                     
                     conditionalPanel(
                         condition = "input.find_markers_tab2 == 1",
                         fluidRow(
                             br(),
                             p("Check Seurat's vignettes", a(tags$a(href="https://satijalab.org/seurat/v3.2/de_vignette.html", "here", target="_blank")), "and the function's",  a(tags$a(href="https://www.rdocumentation.org/packages/Seurat/versions/3.1.4/topics/FindAllMarkers", "manual", target="_blank")),  "to see the options for each parameter of this section. "),
                             br(),
                             column(3,
                                    div(class = "option-group",
                                        radioButtons("find_markers_or_DE_tab2",
                                                     "Detect conserved markers or D.E. genes",
                                                     choices = c("Conserved markers" = 0,
                                                                 "D.E. genes between samples" = 1),
                                                     selected = 0
                                        ))),
                             conditionalPanel(
                                 condition = "input.find_markers_or_DE_tab2 == 0",
                                 column(3,
                                        div(class = "option-group",
                                            radioButtons("find_markers_tab2_opt",
                                                         "Select the analysis to perform?",
                                                         choices = list("Identify conserved markers for all clusters" = 0,
                                                                        "Identify conserved markers for one specific cluster" = 1,
                                                                        "Identify conserved markers for one specific cluster in comparison with another cluster(s)" = 2),
                                                         selected = c(0)),
                                            
                                        )
                                 ),
                                 
                                 conditionalPanel(
                                     condition = "input.find_markers_tab2_opt == 1",
                                     column(3,
                                            div(class = "option-group",
                                                my_withSpinner(  uiOutput("find_markers_clust_id_tab2_ui") )
                                            )),
                                 ),   # ends conditional
                                 conditionalPanel(
                                     condition = "input.find_markers_tab2_opt == 2",
                                     column(3,
                                            div(class = "option-group",
                                                
                                                my_withSpinner( uiOutput("find_markers_clust_ID1_tab2_ui") ),
                                                my_withSpinner( uiOutput("find_markers_clust_ID2_tab2_ui") )
                                                
                                            )),
                                 ),    # ends conditional
                             ),
                             conditionalPanel(
                                 condition = "input.find_markers_or_DE_tab2 == 1",
                                 column(3,
                                        div(class = "option-group",
                                            my_withSpinner( uiOutput("find_markers_or_DE_tab2_cluster_ui") ),
                                            my_withSpinner( uiOutput("find_markers_or_DE_tab2_treat1_ui") ),
                                            my_withSpinner( uiOutput("find_markers_or_DE_tab2_treat2_ui") ),
                                            pickerInput_find_markers_test("find_markers_tab2_test.use"),
                                            
                                            numericInput("find_markers_or_DE_tab2_pvalue",
                                                         label = "Select the (adjusted) p-value threshold.",
                                                         value = 0.05),
                                            
                                        ))),
                             
                             column(3,
                                    actionButtonInput("run_ident_markers_tab2",
                                                      HTML("Search for markers/D.E. genes"))),
                             
                         ), # ends conditional
                     ), # ends fluidrow
                     
                     conditionalPanel(
                         condition = "input.find_markers_tab2 == 1",
                         fluidRow(
                             titlePanel("List of markers or differentially expressed genes"),
                             conditionalPanel(
                                 condition = "input.run_ident_markers_tab2 != 0",
                                 column(12,
                                        my_withSpinner( reactableOutput("markers_tab2_react")) ),
                             ),
                             br(),
                             column(4,
                                    downloadButton("download_markers_tab2",
                                                   "Download the list of markers or D.E. genes."))
                         ), # ends fluidrow
                     ), # ends conditional
                     
                     
                     fluidRow(
                         titlePanel("Expression of markers"),
                         br(),
                         p("In this section, users can visualize the gene expression of selected genes (e.g., tissue markers). Start by loading a", strong("csv"), "file with at least two columns. The first column must be the gene ID, and the second is a grouping variable (e.g., the tissue name). A third column can be used to store the common name of the gene but is optional."),
                         br(),
                         
                         column(width = 3,
                                fileInput_markers_list("markers_list_tab2"),
                                div(class = "option-group",
                                    selectInput("markers_list_header_opt_tab2",
                                                "Does your file have a header?",
                                                choices = c("", "Yes", "No"),
                                                multiple = FALSE,
                                                selected = NULL)),
                                div(class = "option-group",
                                    actionButtonInput("load_markers_tab2",
                                                      HTML("Load markers")))),
                         column(width = 2,
                                conditionalPanel(
                                    condition = "input.load_markers_tab2 > 0",
                                    my_withSpinner( uiOutput('marker_group_selec_tab2') ),
                                    define_if_use_all_genes_or_select("filter_genes_q_tab2"),
                                    
                                    conditionalPanel(
                                        condition = "input.filter_genes_q_tab2 == 0",
                                        define_what_id_to_use("genes_ids_tab2")
                                    ))),
                         conditionalPanel(
                             condition = "input.filter_genes_q_tab2 == 0",
                             
                             column(width = 3,
                                    my_withSpinner( uiOutput('marker_genes_selec_tab2') )
                             )),
                         column(width = 2,
                                conditionalPanel(
                                    condition = "input.load_markers_tab2 > 0",
                                    radioButtons_slot_selection_heatmap("slot_selection_heatmap_tab2"))),
                         column(width = 2,
                                conditionalPanel(
                                    condition = "input.load_markers_tab2 > 0",
                                    actionButtonInput("run_heatmap_tab2",
                                                      HTML("Show heatmap"))))
                         
                         # column(width = 4,
                         #        fileInput_markers_list("markers_list_tab2"),
                         #        div(class = "option-group",
                         #            actionButtonInput("load_markers_tab2",
                         #                              HTML("Load markers")))),
                         # column(width = 3,
                         #        conditionalPanel(
                         #            condition = "input.load_markers_tab2 > 0",
                         #            my_withSpinner( uiOutput('marker_group_selec_tab2') ),
                         #            #my_withSpinner( uiOutput("marker_filter_genes_q_tab2") ),
                         #            #define_what_id_to_use("marker_genes_ids_tab2"),
                         #            define_if_use_all_genes_or_select("marker_filter_genes_q_tab2"),
                         #            conditionalPanel(
                         #                condition = "input.marker_filter_genes_q_tab2 == 0",
                         #                define_what_id_to_use("genes_ids_tab2")
                         #            )
                         #        )),
                         # column(width = 4,
                         #        conditionalPanel(
                         #            condition = "input.marker_filter_genes_q_tab2 == 0",
                         #            my_withSpinner( uiOutput('selected_genes_tab2') )
                         #        )),
                         # column(width = 2,
                         #        conditionalPanel(
                         #            condition = "input.load_markers_tab2 > 0",
                         #
                         #            radioButtons_slot_selection_heatmap("slot_selection_heatmap_tab2"))),
                         # column(width = 2,
                         #        conditionalPanel(
                         #            condition = "input.load_markers_tab2 > 0",
                         #            actionButtonInput("run_heatmap_tab2",
                         #                              HTML("Show heatmap with the average \
                         #                                     of expression per cluster"))))
                     ), # Fluid row
                     
                     #heat map
                     conditionalPanel(
                         
                         condition = "input.run_heatmap_tab2 > 0",
                         fluidRow(
                             br(),
                             
                             titlePanel("Heatmap"),
                             br(),
                             strong("Note that the scale of colors of the heatmap is adjusted based on the expression of the selected genes."),
                             br(),
                             p("For now, the heat map shows the average expression of all samples together. It is only helpful to identify if the markers or cell types make sense with the number of clusters."),
                             p("To visualize the expression in the clusters by sample, use the feature plots."),
                             
                             column(width = 10,
                                    my_withSpinner(  uiOutput("heat_map_ui_tab2") )),
                             # download plot
                             column(2,
                                    numericInput_plot_height("p8_height", value=15),
                                    numericInput_plot_width("p8_width", value=20),
                                    selectInput_plot_res("p8_res"),
                                    selectInput_plot_format("p8_format"),
                                    downloadButton("p8_down", HTML("Download Plot")))
                         ),
                         
                         fluidRow(
                             titlePanel("Visualization of gene expression of each cell and additional plots"),
                             br(),
                             
                             p("Only genes selected for the heatmap can be selected for the additional plots."),
                             br(),
                             column(width = 3,
                                    my_withSpinner( uiOutput("marker_to_feature_plot_tab2") )),
                             column(width = 2,
                                    radioButtons_slot_selection_feature_plot("slot_selection_feature_plot_tab2")),
                             column(width = 3,
                                    actionButtonInput("run_feature_plot_tab2",
                                                      HTML("Show the expression of genes at the cell level")))
                         )
                     ),# Ends conditional
                     
                     conditionalPanel(
                         condition = "input.run_feature_plot_tab2 > 0",
                         fluidRow(
                             titlePanel("Feature plots"),
                             column(width = 7,
                                    my_withSpinner( uiOutput("feature_plot_tab2") )),
                             column(width = 4,
                                    my_withSpinner( uiOutput("umap2_tab2", height = "300px") ))
                         ), # Ends Fluid row
                         
                         fluidRow(
                             titlePanel("Violin and Dot plots"),
                             column(width = 6,
                                    my_withSpinner( uiOutput("run_vln_plot_tab2") )),
                             column(width = 6,
                                    my_withSpinner( uiOutput("run_dot_plot_tab2") ))),
                         fluidRow(
                             titlePanel("Downloading additional plots"),
                             
                             p("In this section, it is possible to download the feature, violin, and dot plots of each selected gene."),
                             p("The files will be saved in the folder:", em("images/one_sample_plots_<current date>__<current time>")),
                             br(),
                             column(width = 2,
                                    radioButtons_down_add_plots("down_add_plots_tab2")),
                             conditionalPanel (
                                 condition = "input.down_add_plots_tab2 == 1",
                                 column(width = 2,
                                        my_withSpinner( uiOutput("select_genes_add_plot_to_down_tab2_ui") ),
                                 ),
                                 
                                 column(2,
                                        div(class = "option-group",
                                            p(strong("Select the options for the feature plots")) ),
                                        numericInput_plot_height("add_p_tab2_feat_height", value=10),
                                        numericInput_plot_width("add_p_tab2_feat_width", value=14),
                                        selectInput_plot_res("add_p_tab2_feat_res"),
                                        selectInput_plot_format("add_p_tab2_feat_format")),
                                 column(2,
                                        div(class = "option-group",
                                            p(strong("Select the options for the violin plots")) ),
                                        numericInput_plot_height("add_p_tab2_violin_height", value=10),
                                        numericInput_plot_width("add_p_tab2_violin_width", value=14),
                                        selectInput_plot_res("add_p_tab2_violin_res"),
                                        selectInput_plot_format("add_p_tab2_violin_format")),
                                 column(2,
                                        div(class = "option-group",
                                            p(strong("Select the options for the dot plots")) ),
                                        numericInput_plot_height("add_p_tab2_dot_height", value=10),
                                        numericInput_plot_width("add_p_tab2_dot_width", value=14),
                                        selectInput_plot_res("add_p_tab2_dot_res"),
                                        selectInput_plot_format("add_p_tab2_dot_format")),
                                 column(width = 3,
                                        actionButtonInput("start_down_add_plots_tab2",
                                                          HTML("Download additional plots")))
                             )
                         ), # end fluid row
                     ), # Ends conditional
                     bookmarkButton(style = "position:absolute;right:2em; background-color:#BF3EFF; color:#FFFFFF;"),
                     
                     tags$hr(),
                     p(strong("Asc-Seurat, version 2.1"), "- Released on May 26th, 2021.", align = "center")
            ), # ends tab
            
            ##############################
            #### trajectory inference ####
            tabPanel("Trajectory inference",
                     br(),
                     p(strong("Note that at any time, you can save a bookmark (purple button at the right bottom) that will save your parameter choices. Using the saved bookmark, it is possible to re-load all selected parameters and re-execute the analysis, to reproduce the results.")),
                     
                     fluidRow(
                         titlePanel("Trajectory inference analysis"),
                         
                         p("For the trajectory inference analyses, users need to use data containing the cluster information obtained in the other tabs of Asc-Seurat. The file must be an RDS file located in the", code("RDS_files/"),"that can be generated using the previous tabs of Asc-Seurat."),
                         br(),
                         p("For this analysis, it is possible to indicate what cluster is expected to be at the beginning and/or end of the trajectory. Depending on the selected model, some of this information might be required."),
                         br(),
                         p("To start the analysis, select the file containing the data and click on", code("Run trajectory inference model"), "button."),
                         column(width = 3,
                                
                                # div(class = "option-group",
                                #     radioButtons("rds_location_tab3",
                                #                  "How to load the rds file?",
                                #                  choices = list("From RDS_files/ folder" = 0,
                                #                                 "Load from other place" = 1),
                                #                  selected = 0)),
                                # conditionalPanel(
                                #     condition = "input.rds_location_tab3 == 0",
                                
                                my_withSpinner( uiOutput("load_integrated_ui_tab3") ),
                                # ),
                                # conditionalPanel(
                                #     condition = "input.rds_location_tab3 == 1",
                                #
                                #     div(class = "option-group",
                                #         fileInput("file_input_rds_tab3",
                                #                   label = "Input the rds file",
                                #                   accept = c(".rds")
                                #         )
                                #     )
                                # ),
                                div(class = "option-group",
                                    radioButtons("ti_select_models",
                                                 "Do you want to run slingshot or another dynverse model?",
                                                 choices = list("Slingshot" = 0,
                                                                "Dynverse model (relies on docker)" = 1
                                                 ),
                                                 selected = c(0)
                                    ))
                         ),
                         column(width = 3,
                                div(class = "option-group",
                                    radioButtons("ti_sample_number",
                                                 "Are you using more than one (integrated) sample ?",
                                                 choices = list("Yes" = 1,
                                                                "No" = 0),
                                                 selected = c(0)
                                    ),
                                )
                         ),
                         
                         column(width = 3,
                                div(class = "option-group",
                                    radioButtons("set_init_clust",
                                                 "Do you want to set an initial and/or ending cluster?",
                                                 choices = list("Yes" = 1,
                                                                "No" = 0),
                                                 selected = c(0)
                                    ))),
                         
                         conditionalPanel (
                             condition = "input.set_init_clust == 1",
                             column(width = 2,
                                    div(class = "option-group",
                                        numericInput("traj_init_clusters",
                                                     label = "Set the initial cluster",
                                                     value = "")
                                    )),
                             column(width = 2,
                                    div(class = "option-group",
                                        numericInput("traj_end_clusters",
                                                     label = "Set the ending cluster",
                                                     value = ""))),
                         ),
                         
                         conditionalPanel (
                             condition = "input.ti_select_models == 1",
                             
                             column(width = 3,
                                    my_withSpinner( uiOutput("ti_methods_list_ui") )),
                             
                             
                         ),
                         column(width = 3,
                                actionButtonInput("run_ti_model",
                                                  HTML("Run trajectory inference model"))),
                         
                     ), # end fluidRow
                     
                     conditionalPanel (
                         condition = "input.run_ti_model != 0",
                         fluidRow(
                             
                             titlePanel("Visualization of the inferred trajectory"),
                             
                             column(width = 4,
                                    my_withSpinner( plotOutput("ti_order")) ),
                             
                             column(width = 3,
                                    my_withSpinner( plotOutput("ti_traject")) ),
                             
                             column(width = 3,
                                    my_withSpinner( plotOutput("ti_graph")) ),
                             
                             conditionalPanel (
                                 condition = "input.ti_sample_number == 1",
                                 column(width = 2,
                                        div(class = "option-group",
                                            radioButtons("ti_graphs_color_choice",
                                                         "Color plots by",
                                                         choices = list("Clusters" = 1,
                                                                        "Samples (treatment)" = 0),
                                                         selected = c(1)
                                            ))),
                             ), #ends conditional
                             
                         ), # end fluidRow
                         
                         fluidRow(
                             br(),
                             # download plot
                             column(width = 3,
                                    div(class = "down-group",
                                        radioButtons("p9_down_opt",
                                                     "Select the plot to download",
                                                     choices = list("Dimension reduction" = 0,
                                                                    "Dendrogram" = 1,
                                                                    "Graph"= 2),
                                                     selected = 0))),
                             column(2,
                                    numericInput_plot_height("p9_height", value=10),
                                    numericInput_plot_width("p9_width", value=12)),
                             column(2,
                                    selectInput_plot_res("p9_res"),
                                    selectInput_plot_format("p9_format")),
                             column(2,
                                    downloadButton("p9_down", HTML("Download Plot"))),
                         ), # Fluid row
                     ), #ends conditional
                     
                     br(),
                     br(),
                     conditionalPanel (
                         condition = "input.run_ti_model != 0",
                         
                         conditionalPanel (
                             condition = "input.ti_select_models == 0",
                             
                             p("More than one lineage of cells can be present in the inferred trajectory. IIf that is the case, the order of the clusters representing the development pathway of each lineage will be shown below. Clusters that appear in multiple rows represent segments of the trajectory that are shared among lineages."),
                             
                             fluidRow(
                                 column(12,
                                        my_withSpinner( verbatimTextOutput("lineages") ))
                             ),
                         ),
                     ),
                     br(),
                     fluidRow(
                         
                         titlePanel("Gene expression in the trajectory"),
                         
                         p("In this section, it is possible to observe the expression of genes in the cells that are part of each lineage/trajectory. You can either provide a list of genes of interest or use the package dynfeature to search for genes that are \"important\" in explaining the trajectory."),
                         br(),
                         p("These genes can be identified in three different ways. 1) Global overview: genes that are important to define the whole trajectory; 2) Lineage/branch: genes that are important to define a branch interest; 3) Genes that are important to define the bifurcation points (points where the trajectory splits into different branches)"),
                         br(),
                         p(strong("Note."), "Dynfeature does not calculate a p-value for the tested genes. Instead, it defines an \"importance value\" that can be used to rank the genes by their relevance."),
                         br(),
                         column(width = 3,
                                div(class = "option-group",
                                    radioButtons("ti_expre_opt",
                                                 "Select the type of analysis:",
                                                 choices = list("Visualization of a list of genes" = 0,
                                                                "Identification of important genes using Dynfeature" = 1
                                                 ),
                                                 selected = c(0)
                                    ))),
                     ), # ends fluidRow
                     
                     ##################################################
                     ### Expression plots - Using markers as input ####
                     ##################################################
                     
                     conditionalPanel(
                         condition = "input.ti_expre_opt == 0",
                         
                         fluidRow(
                             titlePanel("Expression of markers"),
                             br(),
                             # column(width = 4,
                             #        fileInput_markers_list("markers_list_tab3"),
                             #        div(class = "option-group",
                             #            actionButtonInput("load_markers_tab3",
                             #                              HTML("Load markers")))),
                             # column(width = 3,
                             #        conditionalPanel(
                             #            condition = "input.load_markers_tab3 > 0",
                             #            my_withSpinner( uiOutput('marker_group_selec_tab3') ),
                             #            my_withSpinner( uiOutput("marker_filter_genes_q_tab3") ),
                             #            define_what_id_to_use("marker_genes_ids_tab3")
                             #            #my_withSpinner( uiOutput("marker_genes_ids_tab3") )
                             #        )),
                             # column(width = 4,
                             #        conditionalPanel(
                             #            condition = "input.filter_genes_q_tab3 == 0",
                             #            my_withSpinner( uiOutput('marker_genes_selec_tab3') )
                             #        )),
                             # column(width = 2,
                             #        conditionalPanel(
                             #            condition = "input.load_markers_tab3 > 0",
                             #            actionButtonInput("run_heatmap_tab3",
                             #                              HTML("Show heatmap demonstrating the expression \
                             #                                 within the trajectory"))))
                             column(width = 3,
                                    fileInput_markers_list("markers_list_tab3"),
                                    div(class = "option-group",
                                        selectInput("markers_list_header_opt_tab3",
                                                    "Does your file have a header?",
                                                    choices = c("", "Yes", "No"),
                                                    multiple = FALSE,
                                                    selectize = TRUE,
                                                    selected = NULL)),
                                    div(class = "option-group",
                                        actionButtonInput("load_markers_tab3",
                                                          HTML("Load markers")))),
                             column(width = 2,
                                    conditionalPanel(
                                        condition = "input.load_markers_tab3 > 0",
                                        my_withSpinner( uiOutput('marker_group_selec_tab3') ),
                                        define_if_use_all_genes_or_select("filter_genes_q_tab3"),
                                        
                                        conditionalPanel(
                                            condition = "input.filter_genes_q_tab3 == 0",
                                            define_what_id_to_use("genes_ids_tab3")
                                        ))),
                             conditionalPanel(
                                 condition = "input.filter_genes_q_tab3 == 0",
                                 
                                 column(width = 3,
                                        my_withSpinner( uiOutput('marker_genes_selec_tab3') )
                                 )),
                             column(width = 2,
                                    conditionalPanel(
                                        condition = "input.load_markers_tab3 > 0",
                                        radioButtons_slot_selection_heatmap("slot_selection_heatmap_tab3"))),
                             column(width = 2,
                                    conditionalPanel(
                                        condition = "input.load_markers_tab3 > 0",
                                        actionButtonInput("run_heatmap_tab3",
                                                          HTML("Show heatmap"))))
                         ), # Fluid row
                         
                     ), #ends conditional
                     
                     conditionalPanel(
                         
                         condition = "input.run_heatmap_tab3 > 0 & input.ti_expre_opt == 0",
                         fluidRow(
                             br(),
                             
                             titlePanel("Heatmap of the trajectory"),
                             br(),
                             column(width = 10,
                                    my_withSpinner( uiOutput("heat_map_ui_tab3") )),
                             column(2,
                                    numericInput_plot_height("p10_height", value=15),
                                    numericInput_plot_width("p10_width", value=25),
                                    selectInput_plot_res("p10_res"),
                                    selectInput_plot_format("p10_format"),
                                    downloadButton("p10_down", HTML("Download Plot"))),
                         ),
                         fluidRow(
                             titlePanel("Visualization of gene expression of each cell in the trajectory"),
                             br(),
                             column(width = 3,
                                    my_withSpinner( uiOutput("marker_to_feature_plot_tab3") )),
                             column(width = 3,
                                    actionButtonInput("run_feature_plot_tab3",
                                                      HTML("Show the expression of genes at the cell level")))
                             
                         )
                     ),# Ends conditional
                     
                     conditionalPanel (
                         condition = "input.run_feature_plot_tab3 > 0 & input.ti_expre_opt == 0",
                         fluidRow(
                             titlePanel("Feature plots"),
                             
                             column(width = 4,
                                    my_withSpinner( uiOutput("ti_order_express")) ),
                             column(width = 4,
                                    my_withSpinner( uiOutput("ti_traject_express",
                                                             height = "300px")) ),
                             column(width = 4,
                                    my_withSpinner( uiOutput("ti_graph_express",
                                                             height = "300px")) )
                         ), # Ends Fluid row
                         
                     ), # Ends conditional
                     
                     #################################################################
                     ### Expression plots - Using dynverse most relevant features ####
                     #################################################################
                     
                     conditionalPanel(
                         condition = "input.ti_expre_opt == 1",
                         
                         fluidRow (
                             titlePanel("Expression of markers"),
                             br(),
                             p("For each of the below analyses, it is possible to save the list of all expressed genes with their \"importance\" values. It is also possible to save the list of the top N genes."),
                             br(),
                             p("For the branch/linage analysis, you can provide the beginning (from) and end (to) of the branch/linage of interest. You can also select only one, either from or to, as long as the other box is blank."),
                             br(),
                             p("For the bifurcation analysis, you can provide the cluster number where a branching happened to select the genes more relevant to this point. If left in blank, the analysis will identify the \"importance\" values for each branching point."),
                             br(),
                             
                             column(width = 3,
                                    div(class = "option-group",
                                        radioButtons("dynverse_opt",
                                                     "Select the analysis to perform",
                                                     choices = list("global overview of the most predictive genes" = 0,
                                                                    "Lineage/branch markers" = 1,
                                                                    "Genes important at bifurcation points" = 2),
                                                     selected = c(0) ),
                                        numericInput("dynverse_n_genes",
                                                     label = "Set the number of genes to be used in the heatmap (top relevant genes)",
                                                     value = 50) )),
                             
                             conditionalPanel ("input.dynverse_opt == 1",
                                               column(width = 2,
                                                      
                                                      my_withSpinner( uiOutput("dynverse_branch_from_ui") ),
                                                      my_withSpinner( uiOutput("dynverse_branch_to_ui") )),
                                               
                             ), # ends conditional
                             
                             conditionalPanel ("input.dynverse_opt == 2",
                                               column(width = 2,
                                                      div(class = "option-group",
                                                          numericInput("branching_milestone",
                                                                       label = "Select the branching point",
                                                                       value = NULL))),
                             ), # ends conditional
                             column(width = 4,
                                    actionButtonInput("dynverse_def_imp_genes",
                                                      HTML("Defines the most important genes")),
                                    
                                    conditionalPanel("input.dynverse_def_imp_genes > 0",
                                                     #column(width = 4,
                                                     downloadButton("download_dynverse_genes",
                                                                    HTML("Download the list of all important <br> genes in all branches <br> and their \"importance\" score."))
                                    )
                             ),
                         ), # Fluid row
                         
                     ), #ends conditional
                     
                     #heat map
                     conditionalPanel(
                         
                         condition = "input.dynverse_def_imp_genes > 0 & input.ti_expre_opt == 1",
                         fluidRow(
                             br(),
                             
                             titlePanel("Heatmap of the trajectory"),
                             
                             br(),
                             column(width = 10,
                                    my_withSpinner( uiOutput("heat_map_ui_tab3_dynv") )),
                             column(2,
                                    actionButtonInput("run_heatmap_tab3_dynverse",
                                                      HTML("Update heatmap")),
                                    numericInput_plot_height("p11_height", value=15),
                                    numericInput_plot_width("p11_width", value=25),
                                    selectInput_plot_res("p11_res"),
                                    selectInput_plot_format("p11_format"),
                                    downloadButton("p11_down", HTML("Download Plot")),
                                    downloadButton("download_dynverse_genes_filt",
                                                   HTML("Download the fittered list <br> of genes and their <br> \"importance\" score."))
                             ),
                         ), # end row
                         fluidRow(
                             titlePanel("Additional plots"),
                             br(),
                             
                             column(width = 3,
                                    my_withSpinner( uiOutput("marker_to_feature_plot_tab3_dynv") )),
                             column(width = 3,
                                    actionButtonInput("run_feature_plot_tab3_dynv",
                                                      HTML("Show the expression of genes at the cell level")))
                         )
                     ),# Ends conditional
                     
                     conditionalPanel (
                         condition = "input.run_feature_plot_tab3_dynv > 0 & input.ti_expre_opt == 1",
                         fluidRow(
                             titlePanel("Gene expression within the trajectory"),
                             
                             column(width = 4,
                                    my_withSpinner( uiOutput("ti_order_express_dynv")) ),
                             column(width = 4,
                                    my_withSpinner( uiOutput("ti_traject_express_dynv",
                                                             height = "300px")) ),
                             column(width = 4,
                                    my_withSpinner( uiOutput("ti_graph_express_dynv",
                                                             height = "300px")) )
                         ), # Ends Fluid row
                         
                     ), # Ends conditional
                     
                     
                     conditionalPanel (
                         condition = "input.run_feature_plot_tab3 > 0",
                         
                         fluidRow(
                             titlePanel("Downloading plots of gene expression within the trajectory"),
                             
                             p("The files will be saved in the folder:", em("images/one_sample_plots_<current date>__<current time>")),
                             br(),
                             column(width = 2,
                                    radioButtons_down_add_plots("down_add_plots_tab3")),
                             conditionalPanel (
                                 condition = "input.down_add_plots_tab3 == 1",
                                 column(width = 2,
                                        my_withSpinner( uiOutput("select_genes_add_plot_to_down_tab3_ui") ),
                                 ),
                                 
                                 column(2,
                                        p(strong("Select the options for the dimension reduction plots")),
                                        numericInput_plot_height("add_p_tab3_feat_height", value=10),
                                        numericInput_plot_width("add_p_tab3_feat_width", value=14),
                                        selectInput_plot_res("add_p_tab3_feat_res"),
                                        selectInput_plot_format("add_p_tab3_feat_format")),
                                 column(2,
                                        p(strong("Select the options for the dendrogram plots")),
                                        numericInput_plot_height("add_p_tab3_violin_height", value=10),
                                        numericInput_plot_width("add_p_tab3_violin_width", value=18),
                                        selectInput_plot_res("add_p_tab3_violin_res"),
                                        selectInput_plot_format("add_p_tab3_violin_format")),
                                 
                                 column(2,
                                        p(strong("Select the options for the graph plots")),
                                        numericInput_plot_height("add_p_tab3_dot_height", value=10),
                                        numericInput_plot_width("add_p_tab3_dot_width", value=14),
                                        selectInput_plot_res("add_p_tab3_dot_res"),
                                        selectInput_plot_format("add_p_tab3_dot_format")),
                                 column(width = 3,
                                        actionButtonInput("start_down_add_plots_tab3",
                                                          HTML("Download additional plots")))
                                 
                             ) # ends conditional
                         ), # ends fluid row
                     ), # ends conditional
                     
                     conditionalPanel (
                         condition = "input.run_feature_plot_tab3_dynv > 0",
                         
                         fluidRow(
                             
                             titlePanel("Downloading plots showing the gene expression within the trajectory"),
                             
                             p("The files will be saved in the folder:", em("images/one_sample_plots_<current date>__<current time>")),
                             br(),
                             column(width = 2,
                                    radioButtons_down_add_plots("down_add_plots_tab3_dynv")),
                             conditionalPanel (
                                 condition = "input.down_add_plots_tab3_dynv == 1",
                                 column(width = 2,
                                        my_withSpinner( uiOutput("select_genes_add_plot_to_down_tab3_ui_dynv") ),
                                 ),
                                 
                                 column(2,
                                        p(strong("Select the options for the dimension reduction plots")),
                                        numericInput_plot_height("add_p_tab3_feat_height_dynv", value=10),
                                        numericInput_plot_width("add_p_tab3_feat_width_dynv", value=14),
                                        selectInput_plot_res("add_p_tab3_feat_res_dynv"),
                                        selectInput_plot_format("add_p_tab3_feat_format_dynv")),
                                 column(2,
                                        p(strong("Select the options for the dendrogram plots")),
                                        numericInput_plot_height("add_p_tab3_violin_height_dynv", value=10),
                                        numericInput_plot_width("add_p_tab3_violin_width_dynv", value=14),
                                        selectInput_plot_res("add_p_tab3_violin_res_dynv"),
                                        selectInput_plot_format("add_p_tab3_violin_format_dynv")),
                                 column(2,
                                        p(strong("Select the options for the graph plots")),
                                        numericInput_plot_height("add_p_tab3_dot_height_dynv", value=10),
                                        numericInput_plot_width("add_p_tab3_dot_width_dynv", value=14),
                                        selectInput_plot_res("add_p_tab3_dot_res_dynv"),
                                        selectInput_plot_format("add_p_tab3_dot_format_dynv")),
                                 column(width = 3,
                                        actionButtonInput("start_down_add_plots_tab3_dynv",
                                                          HTML("Download additional plots")))
                                 
                             ) # ends conditional
                         ), # ends fluid row
                     ), # ends conditional
                     
                     bookmarkButton(style = "position:absolute;right:2em; background-color:#BF3EFF; color:#FFFFFF;"),
                     
                     tags$hr(),
                     p(strong("Asc-Seurat, version 2.1"), "- Released on May 26th, 2021.", align = "center")
            ),
            
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
            
            ################################
            ##### Annotation - BioMart #####
            ################################
            
            tabPanel("BioMart annotation",
                     
                     sidebarLayout(
                         sidebarPanel(
                             width = 3,
                             h5("BioMart parameters"),
                             my_withSpinner(uiOutput("dbselection"), hide_ui = T),
                             my_withSpinner(uiOutput("datasetselection"), hide_ui = T),
                             my_withSpinner(uiOutput("filterselection"), hide_ui = T),
                             my_withSpinner(uiOutput("attpageselection"), hide_ui = T),
                             my_withSpinner(uiOutput("attselection"), hide_ui = T),
                             fileInput("genescsv", "Input genes for annotation (csv):",
                                       multiple = FALSE, placeholder = "Please select your csv file containing the gene ids"),
                             uiOutput("geneselection"),
                             radioButtons("go_enrichment_analysis",
                                          "Do you want to perform a GO enrichment analysis?",
                                          choices = list("Yes" = 1,
                                                         "No" = 0),
                                          selected = 0),
                             conditionalPanel(
                                 condition = "input.go_enrichment_analysis == 1",
                                 h5("GO enrichment parameters"),
                                 fileInput("universecsv", "Input gene universe file in csv:",
                                           multiple = FALSE, placeholder = "Please select your csv file containing the gene universe ids")
                             )
                         ),
                         
                         ###########################################
                         ### Main panel set up (usually outputs) ###
                         ###########################################
                         mainPanel(
                             width = 9,
                             h3("About"),
                             markdown("This annotation module is based on the [biomaRt R package](https://www.bioconductor.org/packages/devel/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html) that is tailored to access the data available in Ensembl, or Phytozome, rapidly. The package needs an active internet connection. Users can face delays while loading the complete list of parameters in the sidebar and sometimes not get annotation results due to connection problems to the biomaRt server."),
                             markdown("A basic understanding of how BioMart queries are built is required to select the filters and attributes needed. Visit biomaRt’s [vignettes](https://www.bioconductor.org/packages/devel/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html) for an overview."),
                             markdown("Please make sure that your query species is available in the database and that its dataset is appropriately selected in the dataset selection menu that is available in the sidebar.."),
                             my_withSpinner(uiOutput("biomartresultsheader"), hide_ui = T),
                             my_withSpinner(uiOutput("triggerbutton"), hide_ui = T),
                             HTML(paste("<br/>")),
                             my_withSpinner(DTOutput("biomartresults"), hide_ui = T),
                             conditionalPanel(
                                 condition = "input.go_enrichment_analysis == 1",
                                 h4("GO enrichment analysis:"),
                                 my_withSpinner(uiOutput("gobutton"), hide_ui = T),
                                 HTML(paste("<br/>")),
                                 my_withSpinner(DTOutput("goresults"), hide_ui = T),
                                 HTML(paste("<br/>")),
                                 my_withSpinner(plotOutput("goplot"), hide_ui = T),
                                 HTML(paste("<br/>")),
                                 conditionalPanel(
                                     condition = "output.plotbuilt", # Appears after plot is built
                                     fluidRow(
                                         column(3,
                                                div(class = "down-group",
                                                    selectInput("goplottest",
                                                                "Statistical test to plot",
                                                                choices=list("Kolmogorov-Smirnov"  = "KS",
                                                                             "Fisher"              = "weightFisher"),
                                                                selected="Kolmogorov-Smirnov",
                                                                multiple=FALSE,
                                                                selectize=TRUE),
                                                    numericInput("gontop", "Plot (N) top GOs per category", 5, min = 1, step = 1)
                                                )),
                                         column(3,
                                                numericInput_plot_height("goplotheight", value=10),
                                                numericInput_plot_width("goplotwidth", value=10)),
                                         column(3,
                                                selectInput_plot_res("goplotdpi"),
                                                selectInput_plot_format("goplotdevice")),
                                         downloadButton("goplot_down", HTML("Download Plot"))
                                     )
                                 ),
                                 HTML(paste("<br/>"))
                                 
                             ),
                             
                             bookmarkButton(style = "position:absolute;right:2em; background-color:#BF3EFF; color:#FFFFFF;"),
                             
                             tags$hr(),
                             p(strong("Asc-Seurat, version 2.1"), "- Released on May 26th, 2021.", align = "center")
                         )
                     )
            )
            
        ), ## Closing the ui
        width = 12)
}
