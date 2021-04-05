# Asc-Seurat
# Version 1.0
set.seed(1407)

suppressMessages(require("shiny"))
suppressMessages(require("shinyWidgets") )
suppressMessages(require("DT") )
suppressMessages(require("reactable") )
suppressMessages(require("rclipboard") )
suppressMessages(require("shinycssloaders") )
suppressMessages(require("shinyFeedback") )

# help function
my_withSpinner <- function(ui_object_name = ui_object_name, hide_ui = FALSE) {
    
    shinycssloaders::withSpinner(
        ui_object_name,
        type = 4,
        color = "#6504b0",
        size = 0.75,
        hide.ui = hide_ui
    )
}

function(request) {
    
    shinyFeedback::useShinyFeedback()
    
    ##############################################
    #### Properly load user data using docker ####
    ##############################################
    if (dir.exists('/app/user_work')) {
        setwd('/app/user_work')
    }
    fluidPage(
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
        
        # Application title
        titlePanel(title = "Asc-Seurat - Analytical single-cell Seurat-based web application"),
        br(),
        # titlePanel(title = div( img(src = 'NITFIX_logo.png',
        #                             height = "10%",
        #                             width = "10%",
        #                             style='border-right: 10px solid transparent'),
        #                         "- Single-cell RNA-seq") ),
        
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
                             # height = "%",
                             width = "30%",
                             align = "center",
                             style='border-right: 40px solid transparent'),
                         img(src = 'Universidade_de_Brasilia.png',
                             # height = "15%",
                             width = "20%",
                             align = "center",
                             style='border-right: 40px solid transparent'),
                     ),
                     tags$hr(),
                     p(strong("Asc-Seurat, version 1.0"), "- Released on March 19th, 2021.", align = "center")
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
                         #br(),
                         p(strong("After selecting the parameters, click on the blue button to load the data.")),
                         br(),
                         
                         column(width = 3,
                                my_withSpinner(
                                    uiOutput('select_sample_tab1'))
                         ),
                         column(width = 2,
                                div(class = "option-group",
                                    textInput("proj_name",
                                              label = "Project name",
                                              value = ""))),
                         column(width = 2,
                                div(class = "option-group",
                                    numericInput("min_cells",
                                                 label = "Min. number of cells expressing a gene for the gene to be included",
                                                 value = 3))),
                         column(width = 2,
                                div(class = "option-group",
                                    numericInput("min_features",
                                                 label = "Min. number of genes a cell must express to be included",
                                                 value = 200))),
                         column(width = 3,
                                div(class = "option-group",
                                    textInput("mito_regex",
                                              label = "Common identifier of mitochondrial genes",
                                              value = "MT")),
                                actionButton("load_10X",
                                             HTML("Load data of the selected sample"),
                                             style="background-color: #d6fffe"))
                         
                     ), # fluid row
                     
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
                                    my_withSpinner( uiOutput("min_count_ui"))
                                ),
                                div(class = "option-group",
                                    my_withSpinner( uiOutput("max_count_ui"))
                                ),
                                div(class = "option-group",
                                    numericInput("max_mito_perc",
                                                 label = "Exclude any cell with more than this percentage of transcripts belonging to the mitochondrial genes",
                                                 value = 2)),
                                #column(width = 2,
                                actionButton("run_vinplot",
                                             HTML("Show plot of filtered data"),
                                             style="background-color: #d6fffe"))
                         #)
                         
                     ), # fluid row
                     conditionalPanel(
                         condition = "input.run_vinplot!= 0",
                         
                         fluidRow(
                             titlePanel("Screening plot showing the remaining cells after filtering"),
                             column(width = 10,
                                    my_withSpinner( plotOutput("VlnPlot_filt"))
                             ),
                             # download plot
                             column(2,
                                    div(class = "down-group",
                                        numericInput("p1_height","Height (cm)",step=0.5,value=15),
                                        numericInput("p1_width","Width (cm)",step=0.5,value=25),
                                        selectInput("p1_res","Resolution (DPI)",choices=c("100","200","300","400","500","600"),
                                                    selected="100"),
                                        selectInput("p1_format",
                                                    "File type",
                                                    choices=list("png" = "png",
                                                                 "tiff" = "tiff",
                                                                 "jpeg" = "jpeg",
                                                                 "pdf" = "pdf",
                                                                 "svg" = "svg"),
                                                    selected="png",
                                                    multiple=FALSE,
                                                    selectize=TRUE),
                                        downloadButton("p1_down", HTML("Download Plot")))),
                             
                         ), # fluid row
                         
                     ), # Ends conditional
                     
                     fluidRow(
                         titlePanel("Normalizing, centering, and dimension reduction analysis (PCA)"),
                         br(),
                         p("Note that only the most variable genes are used in the dimension reduction step (PCA)."),
                         br(),
                         column(width = 4,
                                div(class = "option-group",
                                    radioButtons("normaliz_method",
                                                 "Select the normalization method",
                                                 choices = list(#"SCTransform (recommended)" = "SCTransform",
                                                     "LogNormalize" = "LogNormalize"),
                                                 selected = c("LogNormalize")))),
                         column(width = 4,
                                div(class = "option-group",
                                    numericInput("scale_factor",
                                                 label = "Scale factor",
                                                 value = 10000))),
                         column(width = 4,
                                div(class = "option-group",
                                    radioButtons("most_var_method",
                                                 "Select the method to detect the most variable genes",
                                                 choices = list("vst" = "vst",
                                                                "mean.var.plot (mvp)" = "mvp",
                                                                "dispersion (disp)" = "disp"),
                                                 selected = c("vst")),
                                    numericInput("n_of_var_genes",
                                                 label = "Number of variable genes for PCA",
                                                 value = 3000)),
                                actionButton("run_pca",
                                             HTML("Run the PCA analysis"),
                                             style="background-color: #d6fffe"))
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
                                    div(class = "down-group",
                                        numericInput("p2_height","Height (cm)",step=0.5,value=12),
                                        numericInput("p2_width","Width (cm)",step=0.5,value=18),
                                        selectInput("p2_res","Resolution (DPI)",choices=c("100","200","300","400","500","600"),
                                                    selected="100"),
                                        selectInput("p2_format",
                                                    "File type",
                                                    choices=list("png" = "png",
                                                                 "tiff" = "tiff",
                                                                 "jpeg" = "jpeg",
                                                                 "pdf" = "pdf",
                                                                 "svg" = "svg"),
                                                    selected="png",
                                                    multiple=FALSE,
                                                    selectize=TRUE),
                                        downloadButton("p2_down", HTML("Download Plot")))),
                             
                             column(width = 3,
                                    div(class = "option-group",
                                        numericInput("n_of_PCs",
                                                     label = "Select the number of components",
                                                     value = 30)))
                             
                         ) # Fluid row
                     ), # Ends conditional
                     
                     fluidRow(
                         titlePanel("Clustering"),
                         br(),
                         p("Be aware that this parameter is central in the cluster definition. It is recommended to try different values and define the most appropriate according to the expectations of the cell populations present in the sample."),
                         p("Quoting from", a(tags$a(href="https://satijalab.org/seurat/archive/v1.4/pbmc3k_tutorial.html", "Seurat's tutorial", target="_blank")), em("\"We find that setting this parameter between 0.6-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets.\"")),
                         br(),
                         column(width = 3,
                                div(class = "option-group",
                                    numericInput("resolution_clust",
                                                 label =
                                                     "Select the resolution for clustering",
                                                 value = 1.2,
                                                 step = 0.1))),
                         column(width = 3,
                                actionButton("run_clustering",
                                             HTML("Run the clustering analysis"),
                                             style="background-color: #d6fffe"))
                         
                     ), # Fluid row
                     
                     conditionalPanel(
                         condition = "input.run_clustering!= 0",
                         
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
                                                     selected = c("UMAP")))),
                             column(2,
                                    div(class = "down-group",
                                        numericInput("p3_height","Height (cm)",step=0.5,value=10),
                                        numericInput("p3_width","Width (cm)",step=0.5,value=12))),
                             column(2,
                                    div(class = "down-group",
                                        selectInput("p3_res","Resolution (DPI)",choices=c("100","200","300","400","500","600"),
                                                    selected="100"),
                                        selectInput("p3_format",
                                                    "File type",
                                                    choices=list("png" = "png",
                                                                 "tiff" = "tiff",
                                                                 "jpeg" = "jpeg",
                                                                 "pdf" = "pdf",
                                                                 "svg" = "svg"),
                                                    selected="png",
                                                    multiple=FALSE,
                                                    selectize=TRUE))),
                             column(2,
                                    div(class = "down-group",
                                        downloadButton("p3_down", HTML("Download Plot")))),
                         ),
                         
                         fluidRow(
                             titlePanel("Number of cells per cluster"),
                             p("First row shows the cluster ID. The second row shows the number of cells per cluster."),
                             column(12,
                                    my_withSpinner( verbatimTextOutput("cluster_size") ))
                         ), # Fluid row
                         
                     ), # Ends conditional
                     
                     fluidRow(
                         titlePanel("Excluding or selecting clusters for reanalysis"),
                         br(),
                         p("Sometimes, it is helpful to exclude or select the clusters that are of interest. After excluding or selecting the cells of interest, it is recommended to repeat the clustering step using only the subset."),
                         p("After selecting the clusters, click on the blue button (Reanalyze after selection/exclusion of clusters). Asc-Seurat will run the analyses of the new subset until the PCA step.", strong("Then, users need to set the new number of components using the elbow plot (above) and click on the button \"Run the clustering analysis\" again.")),
                         br(),
                         column(3,
                                div(class = "option-group",
                                    radioButtons("filter_clusters",
                                                 "Do you want to select or exclude clusters of cells and reanalyze the data?",
                                                 choices = list("Yes" = 1,
                                                                "No" = 0),
                                                 selected = 0),
                                    
                                    
                                )
                         ),
                         
                         conditionalPanel(
                             condition = "input.filter_clusters != 0",
                             
                             column(3,
                                    div(class = "option-group",
                                        radioButtons("filter_clusters_opt",
                                                     "Do you want to select or exclude the clusters?",
                                                     choices = list("Select" = "select",
                                                                    "Exclude" = "exclude"),
                                                     selected = c("select")),
                                        
                                        
                                    )
                             ),
                             
                             column( 3,
                                     div(class = "option-group",
                                         my_withSpinner( uiOutput("cluster_list_ui") )
                                     )
                             ),
                             
                             column( 3,
                                     div(class = "option-group",
                                         actionButton("rerun_after_filtering",
                                                      HTML("Reanalyze after selection/exclusion of clusters"),
                                                      style="background-color: #d6fffe")
                                     )
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
                         )),
                     
                     fluidRow(
                         titlePanel("Identification of markers / Differential expression analysis"),
                         column(3,
                                div(class = "option-group",
                                    radioButtons("find_markers_tab1",
                                                 "Do you want to identify markers or D.E. genes?",
                                                 choices = list("Yes" = 1,
                                                                "No" = 0),
                                                 selected = c(0)),
                                    
                                )
                         )),
                     
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
                                        numericInput("find_markers_tab1_logfc.threshold",
                                                     label = "Select the (log) fold change threshold",
                                                     value = 0.25),
                                        numericInput("find_markers_tab1_min.pct",
                                                     label = "Select the minimal percentage of cells expressing a gene for this gene to be tested",
                                                     value = 0.1),
                                        pickerInput(
                                            inputId = "find_markers_tab1_test.use",
                                            label = "Select the statistical test",
                                            choices = c("wilcox", "bimod", "roc", "t"),
                                            selected = "wilcox",
                                            multiple = FALSE,
                                            options = list(`actions-box` = TRUE)),
                                        conditionalPanel(
                                            condition = "input.find_markers_tab1_opt != 0",
                                            radioButtons("find_markers_tab1_filt_pos",
                                                         "Filter only positive markers?",
                                                         choices = c("yes" = "TRUE",
                                                                     "no" = "FALSE"),
                                                         selected = "TRUE")),
                                        numericInput("find_markers_tab1_return.thresh",
                                                     label = "Select the (adjusted) p-value threshold. Leave it blank if no filtering is required",
                                                     value = "",
                                                     step = 0.01),
                                    ),
                                    actionButton("run_ident_markers_tab1",
                                                 HTML("Search for markers/D.E. genes"),
                                                 style="background-color: #d6fffe")),
                             
                         ), # ends conditional
                     ), # ends fluidrow
                     
                     conditionalPanel(
                         condition = "input.find_markers_tab1 == 1",
                         fluidRow(
                             titlePanel("List of markers or differentially expressed genes"),
                             column(12,
                                    my_withSpinner( reactableOutput("markers_tab1_react")) ),
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
                         column(width = 4,
                                div(class = "option-group",
                                    fileInput("markers_list",
                                              label = "Input the list of markers in the format: 'GeneID,Group,Name' (no header)",
                                              accept = c("text/csv",
                                                         "text/comma-separated-values,text/plain",
                                                         ".csv"))
                                ),
                                div(class = "option-group",
                                    actionButton("load_markers",
                                                 HTML("Load markers"),
                                                 style="background-color: #d6fffe"))),
                         column(width = 3,
                                conditionalPanel(
                                    condition = "input.load_markers > 0",
                                    my_withSpinner( uiOutput('marker_group_selec') ),
                                    my_withSpinner( uiOutput("marker_filter_genes_q") ),
                                    my_withSpinner( uiOutput("marker_genes_ids"))
                                )),
                         column(width = 4,
                                conditionalPanel(
                                    condition = "input.filter_genes_q == 0",
                                    my_withSpinner( uiOutput('marker_genes_selec') )
                                )),
                         column(width = 2,
                                conditionalPanel(
                                    condition = "input.load_markers > 0",
                                    div(class = "option-group",
                                        radioButtons("slot_selection_heatmap",
                                                     "Select the expression values to show",
                                                     choices = list("counts" = "counts",
                                                                    "normalized" = "data",
                                                                    "normalized and scaled" = "scale.data"),
                                                     selected = c("scale.data"))))),
                         column(width = 2,
                                conditionalPanel(
                                    condition = "input.load_markers > 0",
                                    actionButton("run_heatmap",
                                                 HTML("Show heatmap with the average \
                                                             of expression per cluster"),
                                                 style="background-color: #d6fffe")))
                         
                     ), # Fluid row
                     
                     #heat map
                     conditionalPanel(
                         
                         condition = "input.run_heatmap > 0",
                         fluidRow(
                             br(),
                             #       p(strong("obs."), "The heatmap will crash if only one (or zero) gene from your list is found in the dataset. A fix is in progress for the next update. For now, use the feature plots if you want to check a single gene."),
                             #       br(),
                             titlePanel("Heatmap"),
                             br(),
                             strong("Note that the scale of colors of the heatmap is adjusted based on the expression of the selected genes."),
                             br(),
                             column(width = 10,
                                    my_withSpinner( uiOutput("heat_map_ui"))),
                             # download plot
                             column(2,
                                    div(class = "down-group",
                                        numericInput("p4_height","Height (cm)",step=0.5,value=15),
                                        numericInput("p4_width","Width (cm)",step=0.5,value=20),
                                        selectInput("p4_res","Resolution (DPI)",choices=c("100","200","300","400","500","600"),
                                                    selected="100"),
                                        selectInput("p4_format",
                                                    "File type",
                                                    choices=list("png" = "png",
                                                                 "tiff" = "tiff",
                                                                 "jpeg" = "jpeg",
                                                                 "pdf" = "pdf",
                                                                 "svg" = "svg"),
                                                    selected="png",
                                                    multiple=FALSE,
                                                    selectize=TRUE),
                                        downloadButton("p4_down", HTML("Download Plot")))),
                             
                             column(width = 10,
                                    verbatimTextOutput("test"))
                             
                         ),
                         fluidRow(
                             titlePanel("Visualization of gene expression of each cell and additional plots"),
                             br(),
                             # p("Note that, depending of the capacity of your computer, the webpage can crash if you select many genes at once. We recommend selecting no more than 20 genes, so the page scrolls smoothly."),
                             # br(),
                             p("Only genes selected for the heatmap can be selected for the additional plots."),
                             br(),
                             column(width = 3,
                                    my_withSpinner( uiOutput("marker_to_feature_plot") )),
                             column(width = 2,
                                    div(class = "option-group",
                                        radioButtons("slot_selection_feature_plot",
                                                     "Select the expression values to show",
                                                     choices = list("counts" = "counts",
                                                                    "normalized (recommended)" = "data"),
                                                     selected = c("data")))),
                             column(width = 3,
                                    actionButton("run_feature_plot",
                                                 HTML("Show the expression of genes at the cell level"),
                                                 style="background-color: #d6fffe"))
                             
                         )
                     ),# Ends conditional
                     
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
                                    div(class = "option-group",
                                        radioButtons("down_add_plots_tab1",
                                                     "Do you want to download the plots?",
                                                     choices = list("Yes" = 1,
                                                                    "No" = 0),
                                                     selected = 0))),
                             conditionalPanel (
                                 condition = "input.down_add_plots_tab1 == 1",
                                 column(width = 2,
                                        my_withSpinner( uiOutput("select_genes_add_plot_to_down_ui") ),
                                        #uiOutput("select_genes_add_plot_to_down_ui")
                                 ),
                                 column(2,
                                        div(class = "down-group",
                                            p(strong("Select the options for the feature plots")),
                                            numericInput("add_p_tab1_feat_height","Height (cm)",step=0.5,value=10),
                                            numericInput("add_p_tab1_feat_width","Width (cm)",step=0.5,value=14),
                                            selectInput("add_p_tab1_feat_res","Resolution (DPI)",choices=c("100","200","300","400","500","600"),
                                                        selected="100"),
                                            selectInput("add_p_tab1_feat_format",
                                                        "File type",
                                                        choices=list("png" = "png",
                                                                     "tiff" = "tiff",
                                                                     "jpeg" = "jpeg",
                                                                     "pdf" = "pdf",
                                                                     "svg" = "svg"),
                                                        selected="png",
                                                        multiple=FALSE,
                                                        selectize=TRUE))),
                                 column(2,
                                        div(class = "down-group",
                                            p(strong("Select the options for the violin plots")),
                                            numericInput("add_p_tab1_violin_height","Height (cm)",step=0.5,value=10),
                                            numericInput("add_p_tab1_violin_width","Width (cm)",step=0.5,value=18),
                                            selectInput("add_p_tab1_violin_res","Resolution (DPI)",choices=c("100","200","300","400","500","600"),
                                                        selected="100"),
                                            selectInput("add_p_tab1_violin_format",
                                                        "File type",
                                                        choices=list("png" = "png",
                                                                     "tiff" = "tiff",
                                                                     "jpeg" = "jpeg",
                                                                     "pdf" = "pdf",
                                                                     "svg" = "svg"),
                                                        selected="png",
                                                        multiple=FALSE,
                                                        selectize=TRUE))),
                                 column(2,
                                        div(class = "down-group",
                                            p(strong("Select the options for the dot plots")),
                                            numericInput("add_p_tab1_dot_height","Height (cm)",step=0.5,value=10),
                                            numericInput("add_p_tab1_dot_width","Width (cm)",step=0.5,value=14),
                                            selectInput("add_p_tab1_dot_res","Resolution (DPI)",choices=c("100","200","300","400","500","600"),
                                                        selected="100"),
                                            selectInput("add_p_tab1_dot_format",
                                                        "File type",
                                                        choices=list("png" = "png",
                                                                     "tiff" = "tiff",
                                                                     "jpeg" = "jpeg",
                                                                     "pdf" = "pdf",
                                                                     "svg" = "svg"),
                                                        selected="png",
                                                        multiple=FALSE,
                                                        selectize=TRUE))),
                                 column(width = 3,
                                        actionButton("start_down_add_plots_tab1",
                                                     HTML("Download additional plots"),
                                                     style="background-color: #d6fffe"))
                                 
                             )),
                     ),
                     
                     bookmarkButton(style = "position:absolute;right:2em; background-color:#BF3EFF; color:#FFFFFF;"),
                     
                     
                     tags$hr(),
                     p(strong("Asc-Seurat, version 1.0"), "- Released on March 19th, 2021.", align = "center")
                     # Ends page
                     
                     
            ),
            
            # Set the second tav
            
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
                     p("Alternatively, it is possible to load previously integrated data, saving some time by not running the integration step. For that, save an RDS file containing the integrated data."), #strong("Note that it must be saved before the execution of PCA and clustering steps. Otherwise, some of the functions will not work.")
                     p("After selecting the parameters, click on the blue button to load the data."),
                     br(),
                     fluidRow(
                         column(width = 3,
                                div(class = "option-group",
                                    radioButtons("integration_options",
                                                 "Run a new integration analysis or read a previously saved file?",
                                                 choices = list("Run a new analysis" = 0,
                                                                "Load file" = 1),
                                                 selected = c(0) ))),
                         
                         conditionalPanel(
                             condition = "input.integration_options == 1",
                             column(width = 3,
                                    my_withSpinner( uiOutput("load_integrated_ui") )
                             )),
                         
                         conditionalPanel (
                             condition = "input.integration_options == 0",
                             
                             column(width = 3,
                                    div(class = "option-group",
                                        fileInput("samples_list_integration",
                                                  label = "Read the configuration file containing the samples' information",
                                                  accept = c("text/csv",
                                                             "text/comma-separated-values,text/plain",
                                                             ".csv")))),
                             column(width = 3,
                                    my_withSpinner(  uiOutput('select_sample_tab2') )),
                             column(width = 3,
                                    div(class = "option-group",
                                        radioButtons("most_var_method_integration",
                                                     "Select the method to detect the most variable genes",
                                                     choices = list("vst" = "vst",
                                                                    "mean.var.plot (mvp)" = "mvp",
                                                                    "dispersion (disp)" = "disp"),
                                                     selected = c("vst")))),
                             br(),
                             p("To choose the number of components for the integration, use the largest number of components used for the PCA of a single sample."),
                             br(),
                             column(width = 2,
                                    div(class = "option-group",
                                        textInput("int_project_name",
                                                  label = "Type the project name",
                                                  value = "sc-RNA-seq"))),
                             column(width = 2,
                                    div(class = "option-group",
                                        textInput("int_regex_mito",
                                                  label = "Type regex to detect mitochondrial genes",
                                                  value = "MT"))),
                             column(width = 2,
                                    div(class = "option-group",
                                        numericInput("n_of_var_genes_integration",
                                                     label = "N of variable genes for integration",
                                                     value = 3000))),
                             column(width = 2,
                                    div(class = "option-group",
                                        numericInput("n_of_PCs_integration",
                                                     label = "N of components for the integration",
                                                     value = 30))),
                             
                         ), # ends conditional
                         
                         column(width = 3,
                                actionButton("load_rds_file",
                                             HTML("Load the integrated data <br> or <br> execute a new integration"),
                                             style="background-color: #d6fffe")),
                         
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
                     ),
                     
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
                                                 value = "")),
                                div(class = "option-group",
                                    numericInput("max_mito_perc_tab2",
                                                 label = "Exclude any cell with more than this percentage of transcripts belonging to the mitochondrial genes",
                                                 value = ""))),
                         column(width = 2,
                                actionButton("run_vinplot_tab2",
                                             HTML("Show plot of filtered data"),
                                             style="background-color: #d6fffe"))
                         
                     ), # ends fluidRow
                     
                     # fluid row
                     conditionalPanel(
                         condition = "input.run_vinplot_tab2!= 0",
                         
                         fluidRow(
                             titlePanel("Screening plot showing the remaining cells after filtering"),
                             column(width = 10,
                                    my_withSpinner(  plotOutput("VlnPlot_filt_tab2") )
                             ),
                             column(2,
                                    div(class = "down-group",
                                        numericInput("p5_height","Height (cm)",step=0.5,value=8),
                                        numericInput("p5_width","Width (cm)",step=0.5,value=15),
                                        selectInput("p5_res","Resolution (DPI)",choices=c("100","200","300","400","500","600"),
                                                    selected="100"),
                                        selectInput("p5_format",
                                                    "File type",
                                                    choices=list("png" = "png",
                                                                 "tiff" = "tiff",
                                                                 "jpeg" = "jpeg",
                                                                 "pdf" = "pdf",
                                                                 "svg" = "svg"),
                                                    selected="png",
                                                    multiple=FALSE,
                                                    selectize=TRUE),
                                        downloadButton("p5_down", HTML("Download Plot"))))
                             
                         ), # fluid row
                     ), # Ends conditional
                     
                     fluidRow(
                         titlePanel("Normalizing, centering, and dimension reduction analysis (PCA)"),
                         br(),
                         p("Note that only the most variable genes are used in the dimension reduction step (PCA)."),
                         br(),
                         column(width = 4,
                                div(class = "option-group",
                                    radioButtons("normaliz_method_tab2",
                                                 "Select the normalization method",
                                                 choices = list(#"SCTransform (recommended)" = "SCTransform",
                                                     "LogNormalize" = "LogNormalize"),
                                                 selected = c("LogNormalize")))),
                         column(width = 4,
                                div(class = "option-group",
                                    numericInput("scale_factor_tab2",
                                                 label = "Scale factor",
                                                 value = 10000))),
                         column(width = 4,
                                div(class = "option-group",
                                    radioButtons("most_var_method_tab2",
                                                 "Select the method to detect the most variable genes",
                                                 choices = list("vst" = "vst",
                                                                "mean.var.plot (mvp)" = "mvp",
                                                                "dispersion (disp)" = "disp"),
                                                 selected = c("vst")),
                                    numericInput("n_of_var_genes_tab2",
                                                 label = "Number of variable genes for PCA",
                                                 value = 3000)),
                                actionButton("run_pca_tab2",
                                             HTML("Run the PCA analysis"),
                                             style="background-color: #d6fffe"))
                     ), # Fluid row,
                     
                     conditionalPanel(
                         condition = "input.run_pca_tab2!= 0",
                         
                         fluidRow(
                             titlePanel("Elbow plot showing the variation explained by each component"),
                             br(),
                             p("Choose the number of components that explain the most variation."),
                             br(),
                             column(width = 6,
                                    
                                    my_withSpinner( plotOutput("n_of_PCAs_tab2") )
                                    
                             ),
                             column(2,
                                    div(class = "down-group",
                                        numericInput("p6_height","Height (cm)",step=0.5,value=6),
                                        numericInput("p6_width","Width (cm)",step=0.5,value=10),
                                        selectInput("p6_res","Resolution (DPI)",choices=c("100","200","300","400","500","600"),
                                                    selected="100"),
                                        selectInput("p6_format",
                                                    "File type",
                                                    choices=list("png" = "png",
                                                                 "tiff" = "tiff",
                                                                 "jpeg" = "jpeg",
                                                                 "pdf" = "pdf",
                                                                 "svg" = "svg"),
                                                    selected="png",
                                                    multiple=FALSE,
                                                    selectize=TRUE),
                                        downloadButton("p6_down", HTML("Download Plot")))),
                             column(width = 3,
                                    div(class = "option-group",
                                        numericInput("n_of_PCs_tab2",
                                                     label = "Select the number of PCs",
                                                     value = 30)))
                             
                         ) # Fluid row
                     ), # Ends conditional
                     
                     ## Clustering tab2
                     
                     fluidRow(
                         titlePanel("Clustering"),
                         br(),
                         p("Be aware that this parameter is central in the cluster definition. It is recommended to try different values and define the most appropriate according to the expectations of the cell populations present in your sample."),
                         p("Quoting from", a(tags$a(href="https://satijalab.org/seurat/archive/v1.4/pbmc3k_tutorial.html", "Seurat's tutorial", target="_blank")), ":", em("\"We find that setting this parameter between 0.6-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets.\"")),
                         br(),
                         column(width = 3,
                                div(class = "option-group",
                                    numericInput("resolution_clust_tab2",
                                                 label =
                                                     "Select the resolution for clustering",
                                                 value = 1.2,
                                                 step = 0.1))),
                         column(width = 3,
                                actionButton("run_clustering_tab2",
                                             HTML("Run the clustering analysis"),
                                             style="background-color: #d6fffe"))
                     ), # Fluid row
                     
                     conditionalPanel(
                         condition = "input.run_clustering_tab2!= 0",
                         
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
                                    div(class = "down-group",
                                        numericInput("p7_height","Height (cm)",step=0.5,value=10),
                                        numericInput("p7_width","Width (cm)",step=0.5,value=12))),
                             column(2,
                                    div(class = "down-group",
                                        selectInput("p7_res","Resolution (DPI)", choices=c("100","200","300","400","500","600"),
                                                    selected="100"),
                                        selectInput("p7_format",
                                                    "File type",
                                                    choices=list("png" = "png",
                                                                 "tiff" = "tiff",
                                                                 "jpeg" = "jpeg",
                                                                 "pdf" = "pdf",
                                                                 "svg" = "svg"),
                                                    selected="png",
                                                    multiple=FALSE,
                                                    selectize=TRUE))),
                             column(2,
                                    div(class = "down-group",
                                        downloadButton("p7_down", HTML("Download Plot")))),
                         ),# Fluid row
                         
                         fluidRow(
                             titlePanel("Number of cells per cluster"),
                             p("The first row shows the cluster ID. The second row shows the number of cells per cluster."),
                             column(12,
                                    my_withSpinner( verbatimTextOutput("cluster_size_tab2") ))
                         ), # Fluid row
                         
                     ), # Ends conditional
                     fluidRow(
                         titlePanel("Excluding or selecting clusters for reanalysis"),
                         br(),
                         p("Sometimes, it is helpful to exclude or select the clusters that are more of interest.", "After selecting or excluding the cells of interest, it is recommended to repeat the clustering step using only the subset."),
                         p("After selecting the clusters, click on the blue button (Reanalyze after selection/exclusion of clusters). Asc-Seurat will run the analyses of the new subset until the PCA step.", strong("Then, you will need to set the new number of components using the elbow plot (above) and click on the button \"Run the clustering analysis\" again.")),
                         
                         column(3,
                                div(class = "option-group",
                                    radioButtons("filter_clusters_tab2",
                                                 "Do you want to select or exclude clusters of cells and reanalyze the data?",
                                                 choices = list("Yes" = 1,
                                                                "No" = 0),
                                                 selected = 0),
                                )
                         ),
                         
                         conditionalPanel(
                             condition = "input.filter_clusters_tab2 != 0",
                             
                             column(3,
                                    div(class = "option-group",
                                        radioButtons("filter_clusters_opt_tab2",
                                                     "Do you want to select or exclude the clusters?",
                                                     choices = list("Select" = "select",
                                                                    "Exclude" = "exclude"),
                                                     selected = c("select")),
                                    )
                             ),
                             
                             column( 3,
                                     div(class = "option-group",
                                         my_withSpinner( uiOutput("cluster_list_tab2_ui") )
                                     )
                             ),
                             
                             column( 3,
                                     div(class = "option-group",
                                         actionButton("rerun_after_filtering_tab2",
                                                      HTML("Reanalyze after selection/exclusion of clusters"),
                                                      style="background-color: #d6fffe")
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
                         )),
                     
                     
                     fluidRow(
                         titlePanel("Identification of markers / Differential expression analysis"),
                         column(3,
                                div(class = "option-group",
                                    radioButtons("find_markers_tab2",
                                                 "Do you want to identify markers and/or D.E. genes?",
                                                 choices = list("Yes" = 1,
                                                                "No" = 0),
                                                 selected = c(0)),
                                    
                                )
                         )),
                     
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
                                                                 "D.E. genes between treat vs. control" = 1),
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
                                            pickerInput(
                                                inputId = "find_markers_tab2_test.use",
                                                label = "Select the statistical test",
                                                choices = c("wilcox", "bimod", "roc", "t"), #"negbinom", "poisson", "LR", "MAST" were excluded since they need to set latent.vars and we do not provide options to the user to set it.
                                                selected = "wilcox",
                                                multiple = FALSE,
                                                options = list(`actions-box` = TRUE)),
                                            numericInput("find_markers_or_DE_tab2_pvalue",
                                                         label = "Select the (adjusted) p-value threshold. Leave it blank if no filtering is required",
                                                         value = 0.05),
                                            
                                        ))),
                             
                             column(3,
                                    actionButton("run_ident_markers_tab2",
                                                 HTML("Search for markers/D.E. genes"),
                                                 style="background-color: #d6fffe")),
                             
                         ), # ends conditional
                     ), # ends fluidrow
                     
                     conditionalPanel(
                         condition = "input.find_markers_tab2 == 1",
                         fluidRow(
                             titlePanel("List of markers of differentially expressed genes"),
                             column(12,
                                    my_withSpinner( reactableOutput("markers_tab2_react")) ),
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
                         column(width = 4,
                                div(class = "option-group",
                                    fileInput("markers_list_tab2",
                                              label = "Input the list of markers in the format: 'GeneID,Group,Name' (no header)",
                                              accept = c("text/csv",
                                                         "text/comma-separated-values,text/plain",
                                                         ".csv"))
                                    
                                    
                                ),
                                div(class = "option-group",
                                    actionButton("load_markers_tab2",
                                                 HTML("Load markers"),
                                                 style="background-color: #d6fffe"))),
                         column(width = 3,
                                conditionalPanel(
                                    condition = "input.load_markers_tab2 > 0",
                                    my_withSpinner( uiOutput('marker_group_selec_tab2') ),
                                    my_withSpinner( uiOutput("marker_filter_genes_q_tab2") ),
                                    my_withSpinner( uiOutput("marker_genes_ids_tab2") )
                                )),
                         column(width = 4,
                                conditionalPanel(
                                    condition = "input.filter_genes_q_tab2 == 0",
                                    my_withSpinner( uiOutput('marker_genes_selec_tab2') )
                                )),
                         column(width = 2,
                                conditionalPanel(
                                    condition = "input.load_markers_tab2 > 0",
                                    div(class = "option-group",
                                        radioButtons("slot_selection_heatmap_tab2",
                                                     "Select the expression values to show",
                                                     choices = list("counts" = "counts",
                                                                    "normalized" = "data",
                                                                    "normalized and scaled" = "scale.data"),
                                                     selected = c("scale.data"))))),
                         column(width = 2,
                                conditionalPanel(
                                    condition = "input.load_markers_tab2 > 0",
                                    actionButton("run_heatmap_tab2",
                                                 HTML("Show heatmap with the average \
                                                             of expression per cluster"),
                                                 style="background-color: #d6fffe")))
                         
                     ), # Fluid row
                     
                     #heat map
                     conditionalPanel(
                         
                         condition = "input.run_heatmap_tab2 > 0",
                         fluidRow(
                             br(),
                             #          p(strong("Obs."), "The heatmap will crash if only one (or zero) gene of your list is found in the dataset. I am going to fix that soon. For now, use the feature plots if you want to check a single gene."),
                             #          br(),
                             titlePanel("Heatmap"),
                             br(),
                             strong("Note that the scale of colors of the heatmap is adjusted based on the expression of the selected genes."),
                             br(),
                             p("For now, the heat map shows the average expression of all samples together. It is only helpful to identify if the markers or cell types make sense with the number of clusters."),
                             p("To visualize the expression in the clusters by sample, use the feature plots."),
                             #        p(strong("obs."), "The heatmap will crash if only one (or zero) gene of your list is found in the dataset. I am going to fix that soon. For now, use the feature plots if you want to check a single gene."),
                             column(width = 10,
                                    my_withSpinner(  uiOutput("heat_map_ui_tab2") )),
                             # download plot
                             column(2,
                                    div(class = "down-group",
                                        numericInput("p8_height","Height (cm)",step=0.5,value=15),
                                        numericInput("p8_width","Width (cm)",step=0.5,value=20),
                                        selectInput("p8_res","Resolution (DPI)",choices=c("100","200","300","400","500","600"),
                                                    selected="100"),
                                        selectInput("p8_format",
                                                    "File type",
                                                    choices=list("png" = "png",
                                                                 "tiff" = "tiff",
                                                                 "jpeg" = "jpeg",
                                                                 "pdf" = "pdf",
                                                                 "svg" = "svg"),
                                                    selected="png",
                                                    multiple=FALSE,
                                                    selectize=TRUE),
                                        downloadButton("p8_down", HTML("Download Plot")))),
                             
                         ),
                         
                         
                         
                         
                         fluidRow(
                             titlePanel("Visualization of gene expression of each cell and additional plots"),
                             br(),
                             # p("Note that, depending of the capacity of your computer, the webpage can crash if you select many genes at once. It is recommended to select no more than 20 genes, so the page scrolls smoothly."),
                             # br(),
                             p("Only genes selected for the heatmap can be selected for the additional plots."),
                             br(),
                             column(width = 3,
                                    my_withSpinner( uiOutput("marker_to_feature_plot_tab2") )),
                             column(width = 2,
                                    div(class = "option-group",
                                        radioButtons("slot_selection_feature_plot_tab2",
                                                     "Select the expression values to show",
                                                     choices = list("counts" = "counts",
                                                                    "normalized (recommended)" = "data"),
                                                     selected = c("data")))),
                             column(width = 3,
                                    actionButton("run_feature_plot_tab2",
                                                 HTML("Show the expression of genes at the cell level"),
                                                 style="background-color: #d6fffe"))
                             
                         )
                     ),# Ends conditional
                     
                     conditionalPanel(
                         condition = "input.run_feature_plot_tab2 > 0",
                         fluidRow(
                             titlePanel("Feature plots"),
                             column(width = 8,
                                    my_withSpinner( uiOutput("feature_plot_tab2") )),
                             column(width = 4,
                                    my_withSpinner( uiOutput("umap2_tab2", height = "300px") ))
                         ), # Ends Fluid row
                         
                         fluidRow(
                             titlePanel("Violin and Dot plots"),
                             column(width = 8,
                                    my_withSpinner( uiOutput("run_vln_plot_tab2") )),
                             column(width = 4,
                                    my_withSpinner( uiOutput("run_dot_plot_tab2") ))),
                         fluidRow(
                             titlePanel("Downloading additional plots"),
                             
                             p("In this section, it is possible to download the feature, violin, and dot plots of each selected gene."),
                             p("The files will be saved in the folder:", em("images/one_sample_plots_<current date>__<current time>")),
                             br(),
                             column(width = 2,
                                    div(class = "option-group",
                                        radioButtons("down_add_plots_tab2",
                                                     "Do you want to download the plots?",
                                                     choices = list("Yes" = 1,
                                                                    "No" = 0),
                                                     selected = 0))),
                             conditionalPanel (
                                 condition = "input.down_add_plots_tab2 == 1",
                                 column(width = 2,
                                        my_withSpinner( uiOutput("select_genes_add_plot_to_down_tab2_ui") ),
                                        #uiOutput("select_genes_add_plot_to_down_ui")
                                 ),
                                 column(2,
                                        div(class = "down-group",
                                            p(strong("Select the options for the feature plots")),
                                            numericInput("add_p_tab2_feat_height","Height (cm)",step=0.5,value=10),
                                            numericInput("add_p_tab2_feat_width","Width (cm)",step=0.5,value=14),
                                            selectInput("add_p_tab2_feat_res","Resolution (DPI)",choices=c("100","200","300","400","500","600"),
                                                        selected="100"),
                                            selectInput("add_p_tab2_feat_format",
                                                        "File type",
                                                        choices=list("png" = "png",
                                                                     "tiff" = "tiff",
                                                                     "jpeg" = "jpeg",
                                                                     "pdf" = "pdf",
                                                                     "svg" = "svg"),
                                                        selected="png",
                                                        multiple=FALSE,
                                                        selectize=TRUE))),
                                 column(2,
                                        div(class = "down-group",
                                            p(strong("Select the options for the violin plots")),
                                            numericInput("add_p_tab2_violin_height","Height (cm)",step=0.5,value=10),
                                            numericInput("add_p_tab2_violin_width","Width (cm)",step=0.5,value=18),
                                            selectInput("add_p_tab2_violin_res","Resolution (DPI)",choices=c("100","200","300","400","500","600"),
                                                        selected="100"),
                                            selectInput("add_p_tab2_violin_format",
                                                        "File type",
                                                        choices=list("png" = "png",
                                                                     "tiff" = "tiff",
                                                                     "jpeg" = "jpeg",
                                                                     "pdf" = "pdf",
                                                                     "svg" = "svg"),
                                                        selected="png",
                                                        multiple=FALSE,
                                                        selectize=TRUE))),
                                 column(2,
                                        div(class = "down-group",
                                            p(strong("Select the options for the dot plots")),
                                            numericInput("add_p_tab2_dot_height","Height (cm)",step=0.5,value=10),
                                            numericInput("add_p_tab2_dot_width","Width (cm)",step=0.5,value=14),
                                            selectInput("add_p_tab2_dot_res","Resolution (DPI)",choices=c("100","200","300","400","500","600"),
                                                        selected="100"),
                                            selectInput("add_p_tab2_dot_format",
                                                        "File type",
                                                        choices=list("png" = "png",
                                                                     "tiff" = "tiff",
                                                                     "jpeg" = "jpeg",
                                                                     "pdf" = "pdf",
                                                                     "svg" = "svg"),
                                                        selected="png",
                                                        multiple=FALSE,
                                                        selectize=TRUE))),
                                 column(width = 3,
                                        actionButton("start_down_add_plots_tab2",
                                                     HTML("Download additional plots"),
                                                     style="background-color: #d6fffe"))
                                 
                             )
                         ), # end fluid row
                     ), # Ends conditional
                     bookmarkButton(style = "position:absolute;right:2em; background-color:#BF3EFF; color:#FFFFFF;"),
                     
                     tags$hr(),
                     p(strong("Asc-Seurat, version 1.0"), "- Released on March 19th, 2021.", align = "center")
            ), # ends tab
            
            ##############################
            #### trajectory inference ####
            tabPanel("Trajectory inference",
                     br(),
                     p(strong("Note that at any time, you can save a bookmark (purple button at the right bottom) that will save your parameter choices. Using the saved bookmark, it is possible to re-load all selected parameters and re-execute the analysis, to reproduce the results.")),
                     
                     fluidRow(
                         titlePanel("Trajectory inference analysis"),
                         # br(),
                         # p("In this tab it is possible to execute the trajectory inference (pseudo-time) analysis using slingshot. For more information about the package, click here", a(tags$a(href="https://www.bioconductor.org/packages/release/bioc/html/sli", "here."))),
                         # br(),
                         p("For the trajectory inference analyses, users need to use data containing the cluster information obtained in the other tabs of Asc-Seurat. The file must be an RDS file located in the", code("RDS_files/"),"that can be generated using the previous tabs of Asc-Seurat."),
                         br(),
                         p("For this analysis, it is possible to indicate what cluster is expected to be at the beginning and/or end of the trajectory. Depending on the selected model, some of this information might be required."),
                         br(),
                         p("To start the analysis, select the file containing the data and click on", code("Run trajectory inference model"), "button."),
                         column(width = 3,
                                my_withSpinner( uiOutput("load_integrated_ui_tab3") ),
                                div(class = "option-group",
                                    radioButtons("ti_select_models",
                                                 "Do you want to run slingshot or another dynverse model?",
                                                 choices = list("Slingshot" = 0,
                                                                "Dynverse model (relies on docker)" = 1
                                                 ),
                                                 selected = c(0)
                                    ))),
                         column(width = 3,
                                div(class = "option-group",
                                    radioButtons("ti_sample_number",
                                                 "Are you using more than one (integrated) sample ?",
                                                 choices = list("Yes" = 1,
                                                                "No" = 0),
                                                 selected = c(0)
                                    ),
                                    # conditionalPanel (
                                    #     condition = "input.ti_sample_number == 1",
                                    #     radioButtons("ti_timepoint",
                                    #                  "If using multiple sample, are they timepoints?",
                                    #                  choices = list("Yes, they are timepoints" = 1,
                                    #                                 "No, they are not timepoints" = 0),
                                    #                  selected = c(1)
                                    #     ),
                                    #     conditionalPanel (
                                    #         condition = "input.ti_timepoint == 1",
                                    #         radioButtons("ti_timepoint_cont_disct",
                                    #                      "Continuous or discrete?",
                                    #                      choices = list("Continuous" = 1,
                                    #                                     "Discrete" = 0),
                                    #                      selected = c(1)
                                    #         ),
                                    #     ),
                                    # ),
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
                                actionButton("run_ti_model",
                                             HTML("Run trajectory inference model"),
                                             style="background-color: #d6fffe")),
                         
                     ), # end fluidRow
                     
                     conditionalPanel (
                         condition = "input.run_ti_model != 0",
                         fluidRow(
                             
                             titlePanel("Visualization of the inferred trajectory"),
                             
                             #p("Below, it is possible to visualize the trajectory in three different representation of the data."),
                             
                             
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
                                    div(class = "down-group",
                                        numericInput("p9_height","Height (cm)",step=0.5,value=10),
                                        numericInput("p9_width","Width (cm)",step=0.5,value=12))),
                             column(2,
                                    div(class = "down-group",
                                        selectInput("p9_res","Resolution (DPI)", choices=c("100","200","300","400","500","600"),
                                                    selected="100"),
                                        selectInput("p9_format",
                                                    "File type",
                                                    choices=list("png" = "png",
                                                                 "tiff" = "tiff",
                                                                 "jpeg" = "jpeg",
                                                                 "pdf" = "pdf",
                                                                 "svg" = "svg"),
                                                    selected="png",
                                                    multiple=FALSE,
                                                    selectize=TRUE))),
                             column(2,
                                    div(class = "down-group",
                                        downloadButton("p9_down", HTML("Download Plot")))),
                         ), # Fluid row
                     ), #ends conditional
                     
                     br(),
                     br(),
                     conditionalPanel (
                         condition = "input.run_ti_model != 0",
                         
                         p("More than one lineage of cells can be present in the inferred trajectory. IIf that is the case, the order of the clusters representing the development pathway of each lineage will be shown below. Clusters that appear in multiple rows represent segments of the trajectory that are shared among lineages."),
                         
                         fluidRow(
                             column(12,
                                    my_withSpinner( verbatimTextOutput("lineages") )
                             )
                         ),
                     ),
                     br(),
                     fluidRow(
                         
                         titlePanel("Gene expression in the trajectory"),
                         
                         p("In this section, it is possible to observe the expression of genes in the cells that are part of each lineage/trajectory. You can either provide a list of genes of interest or use the package dynfeature to search for genes that are \"important\" in explaining the trajectory."),
                         br(),
                         p("These genes can be identified in three different ways. 1) Global overview: genes that are important to define the whole trajectory; 2) Lineage/branch: genes that are important to define a branch interest; 3) Genes that are important to define the bifurcation points (points where the trajectory splits into different branches)"),
                         br(),
                         # p("D.E. genes == differentialy expressed genes"),
                         p(strong("Note."), "Dynfeature does not calculate a p-value for the tested genes. Instead, it defines an \"importance value\" that can be used to rank the genes by their relevance."),
                         br(),
                         
                         #p("tradSEq does provide a p-value estimative"),
                         
                         column(width = 3,
                                div(class = "option-group",
                                    radioButtons("ti_expre_opt",
                                                 "Select the type of analysis:",
                                                 choices = list("Visualization of a list of genes" = 0,
                                                                "Identification of important genes using Dynfeature" = 1
                                                                #,"D.E. genes based on tradSEQ" = 2
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
                             column(width = 4,
                                    div(class = "option-group",
                                        fileInput("markers_list_tab3",
                                                  label = "Input the list of markers in the format: 'GeneID,Group,Name' (no header)",
                                                  accept = c("text/csv",
                                                             "text/comma-separated-values,text/plain",
                                                             ".csv"))
                                        
                                        
                                    ),
                                    div(class = "option-group",
                                        actionButton("load_markers_tab3",
                                                     HTML("Load markers"),
                                                     style="background-color: #d6fffe"))),
                             column(width = 3,
                                    conditionalPanel(
                                        condition = "input.load_markers_tab3 > 0",
                                        my_withSpinner( uiOutput('marker_group_selec_tab3') ),
                                        my_withSpinner( uiOutput("marker_filter_genes_q_tab3") ),
                                        my_withSpinner( uiOutput("marker_genes_ids_tab3") )
                                    )),
                             column(width = 4,
                                    conditionalPanel(
                                        condition = "input.filter_genes_q_tab3 == 0",
                                        my_withSpinner( uiOutput('marker_genes_selec_tab3') )
                                    )),
                             column(width = 2,
                                    conditionalPanel(
                                        condition = "input.load_markers_tab3 > 0",
                                        actionButton("run_heatmap_tab3",
                                                     HTML("Show heatmap demonstrating the expression \
                                                             within the trajectory"),
                                                     style="background-color: #d6fffe")))
                             
                         ), # Fluid row
                         
                     ), #ends conditional
                     
                     conditionalPanel(
                         
                         condition = "input.run_heatmap_tab3 > 0 & input.ti_expre_opt == 0",
                         fluidRow(
                             br(),
                             #    p(strong("obs."), "The heatmap will crash if only one (or zero) gene from your list is found in the dataset. A fix is in progress for the next update. For now, use the feature plots if you want to check a single gene."),
                             #    br(),
                             titlePanel("Heatmap of the trajectory"),
                             br(),
                             column(width = 10,
                                    my_withSpinner( uiOutput("heat_map_ui_tab3") )),
                             column(2,
                                    div(class = "down-group",
                                        numericInput("p10_height","Height (cm)",step=0.5,value=15),
                                        numericInput("p10_width","Width (cm)",step=0.5,value=25),
                                        selectInput("p10_res","Resolution (DPI)",choices=c("100","200","300","400","500","600"),
                                                    selected="100"),
                                        selectInput("p10_format",
                                                    "File type",
                                                    choices=list("png" = "png",
                                                                 "tiff" = "tiff",
                                                                 "jpeg" = "jpeg",
                                                                 "pdf" = "pdf",
                                                                 "svg" = "svg"),
                                                    selected="png",
                                                    multiple=FALSE,
                                                    selectize=TRUE),
                                        downloadButton("p10_down", HTML("Download Plot")))),
                         ),
                         fluidRow(
                             titlePanel("Visualization of gene expression of each cell in the trajectory"),
                             br(),
                             # p("Note that, depending of the capacity of your computer, the webpage can crash if you select many genes at once. We recommend selecting no more than 20 genes, so the page scrolls smoothly."),
                             # br(),
                             column(width = 3,
                                    my_withSpinner( uiOutput("marker_to_feature_plot_tab3") )),
                             column(width = 3,
                                    actionButton("run_feature_plot_tab3",
                                                 HTML("Show the expression of genes at the cell level"),
                                                 style="background-color: #d6fffe"))
                             
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
                                    actionButton("dynverse_def_imp_genes",
                                                 HTML("Defines the most important genes"),
                                                 style="background-color: #d6fffe"),
                                    
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
                             #       p(strong("obs."), "The heatmap will crash if only one (or zero) gene from your list is found in the dataset. A fix is in progress for the next update. For now, use the feature plots if you want to check a single gene."),
                             #       br(),
                             #p("Be careful with the name of the genes."),
                             titlePanel("Heatmap of the trajectory"),
                             
                             br(),
                             column(width = 10,
                                    my_withSpinner( uiOutput("heat_map_ui_tab3_dynv") )),
                             column(2,
                                    actionButton("run_heatmap_tab3_dynverse",
                                                 HTML("Update heatmap"),
                                                 style="background-color: #d6fffe"),
                                    
                                    div(class = "down-group",
                                        numericInput("p11_height","Height (cm)",step=0.5,value=15),
                                        numericInput("p11_width","Width (cm)",step=0.5,value=25),
                                        selectInput("p11_res","Resolution (DPI)",choices=c("100","200","300","400","500","600"),
                                                    selected="100"),
                                        selectInput("p11_format",
                                                    "File type",
                                                    choices=list("png" = "png",
                                                                 "tiff" = "tiff",
                                                                 "jpeg" = "jpeg",
                                                                 "pdf" = "pdf",
                                                                 "svg" = "svg"),
                                                    selected="png",
                                                    multiple=FALSE,
                                                    selectize=TRUE),
                                        downloadButton("p11_down", HTML("Download Plot"))),
                                    downloadButton("download_dynverse_genes_filt",
                                                   HTML("Download the fittered list <br> of genes and their <br> \"importance\" score."))
                             ),
                         ), # end row
                         fluidRow(
                             titlePanel("Additional plots"),
                             br(),
                             # p("Note that, depending of the capacity of your computer, the webpage can crash if you select many genes at once. We recommend selecting no more than 20 genes, so the page scrolls smoothly."),
                             # br(),
                             
                             column(width = 3,
                                    my_withSpinner( uiOutput("marker_to_feature_plot_tab3_dynv") )),
                             column(width = 3,
                                    actionButton("run_feature_plot_tab3_dynv",
                                                 HTML("Show the expression of genes at the cell level"),
                                                 style="background-color: #d6fffe"))
                             
                         )
                     ),# Ends conditional
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     conditionalPanel (
                         condition = "input.run_feature_plot_tab3_dynv > 0 & input.ti_expre_opt == 1",
                         fluidRow(
                             #<<<<<<< HEAD
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
                                    div(class = "option-group",
                                        radioButtons("down_add_plots_tab3",
                                                     "Do you want to download the plots?",
                                                     choices = list("Yes" = 1,
                                                                    "No" = 0),
                                                     selected = 0))),
                             conditionalPanel (
                                 condition = "input.down_add_plots_tab3 == 1",
                                 column(width = 2,
                                        my_withSpinner( uiOutput("select_genes_add_plot_to_down_tab3_ui") ),
                                        #uiOutput("select_genes_add_plot_to_down_ui")
                                 ),
                                 column(2,
                                        div(class = "down-group",
                                            p(strong("Select the options for the dimension reduction plots")),
                                            numericInput("add_p_tab3_feat_height","Height (cm)",step=0.5,value=10),
                                            numericInput("add_p_tab3_feat_width","Width (cm)",step=0.5,value=14),
                                            selectInput("add_p_tab3_feat_res","Resolution (DPI)",choices=c("100","200","300","400","500","600"),
                                                        selected="100"),
                                            selectInput("add_p_tab3_feat_format",
                                                        "File type",
                                                        choices=list("png" = "png",
                                                                     "tiff" = "tiff",
                                                                     "jpeg" = "jpeg",
                                                                     "pdf" = "pdf",
                                                                     "svg" = "svg"),
                                                        selected="png",
                                                        multiple=FALSE,
                                                        selectize=TRUE))),
                                 column(2,
                                        div(class = "down-group",
                                            p(strong("Select the options for the dendrogram plots")),
                                            numericInput("add_p_tab3_violin_height","Height (cm)",step=0.5,value=10),
                                            numericInput("add_p_tab3_violin_width","Width (cm)",step=0.5,value=18),
                                            selectInput("add_p_tab3_violin_res","Resolution (DPI)",choices=c("100","200","300","400","500","600"),
                                                        selected="100"),
                                            selectInput("add_p_tab3_violin_format",
                                                        "File type",
                                                        choices=list("png" = "png",
                                                                     "tiff" = "tiff",
                                                                     "jpeg" = "jpeg",
                                                                     "pdf" = "pdf",
                                                                     "svg" = "svg"),
                                                        selected="png",
                                                        multiple=FALSE,
                                                        selectize=TRUE))),
                                 column(2,
                                        div(class = "down-group",
                                            p(strong("Select the options for the graph plots")),
                                            numericInput("add_p_tab3_dot_height","Height (cm)",step=0.5,value=10),
                                            numericInput("add_p_tab3_dot_width","Width (cm)",step=0.5,value=14),
                                            selectInput("add_p_tab3_dot_res","Resolution (DPI)",choices=c("100","200","300","400","500","600"),
                                                        selected="100"),
                                            selectInput("add_p_tab3_dot_format",
                                                        "File type",
                                                        choices=list("png" = "png",
                                                                     "tiff" = "tiff",
                                                                     "jpeg" = "jpeg",
                                                                     "pdf" = "pdf",
                                                                     "svg" = "svg"),
                                                        selected="png",
                                                        multiple=FALSE,
                                                        selectize=TRUE))),
                                 column(width = 3,
                                        actionButton("start_down_add_plots_tab3",
                                                     HTML("Download additional plots"),
                                                     style="background-color: #d6fffe"))
                                 
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
                                    div(class = "option-group",
                                        radioButtons("down_add_plots_tab3_dynv",
                                                     "Do you want to download the plots?",
                                                     choices = list("Yes" = 1,
                                                                    "No" = 0),
                                                     selected = 0))),
                             conditionalPanel (
                                 condition = "input.down_add_plots_tab3_dynv == 1",
                                 column(width = 2,
                                        my_withSpinner( uiOutput("select_genes_add_plot_to_down_tab3_ui_dynv") ),
                                        #uiOutput("select_genes_add_plot_to_down_ui")
                                 ),
                                 column(2,
                                        div(class = "down-group",
                                            p(strong("Select the options for the dimension reduction plots")),
                                            numericInput("add_p_tab3_feat_height_dynv","Height (cm)",step=0.5,value=10),
                                            numericInput("add_p_tab3_feat_width_dynv","Width (cm)",step=0.5,value=14),
                                            selectInput("add_p_tab3_feat_res_dynv","Resolution (DPI)",choices=c("100","200","300","400","500","600"),
                                                        selected="100"),
                                            selectInput("add_p_tab3_feat_format_dynv",
                                                        "File type",
                                                        choices=list("png" = "png",
                                                                     "tiff" = "tiff",
                                                                     "jpeg" = "jpeg",
                                                                     "pdf" = "pdf",
                                                                     "svg" = "svg"),
                                                        selected="png",
                                                        multiple=FALSE,
                                                        selectize=TRUE))),
                                 column(2,
                                        div(class = "down-group",
                                            p(strong("Select the options for the dendrogram plots")),
                                            numericInput("add_p_tab3_violin_height_dynv","Height (cm)",step=0.5,value=10),
                                            numericInput("add_p_tab3_violin_width_dynv","Width (cm)",step=0.5,value=18),
                                            selectInput("add_p_tab3_violin_res_dynv","Resolution (DPI)",choices=c("100","200","300","400","500","600"),
                                                        selected="100"),
                                            selectInput("add_p_tab3_violin_format_dynv",
                                                        "File type",
                                                        choices=list("png" = "png",
                                                                     "tiff" = "tiff",
                                                                     "jpeg" = "jpeg",
                                                                     "pdf" = "pdf",
                                                                     "svg" = "svg"),
                                                        selected="png",
                                                        multiple=FALSE,
                                                        selectize=TRUE))),
                                 column(2,
                                        div(class = "down-group",
                                            p(strong("Select the options for the graph plots")),
                                            numericInput("add_p_tab3_dot_height_dynv","Height (cm)",step=0.5,value=10),
                                            numericInput("add_p_tab3_dot_width_dynv","Width (cm)",step=0.5,value=14),
                                            selectInput("add_p_tab3_dot_res_dynv","Resolution (DPI)",choices=c("100","200","300","400","500","600"),
                                                        selected="100"),
                                            selectInput("add_p_tab3_dot_format_dynv",
                                                        "File type",
                                                        choices=list("png" = "png",
                                                                     "tiff" = "tiff",
                                                                     "jpeg" = "jpeg",
                                                                     "pdf" = "pdf",
                                                                     "svg" = "svg"),
                                                        selected="png",
                                                        multiple=FALSE,
                                                        selectize=TRUE))),
                                 column(width = 3,
                                        actionButton("start_down_add_plots_tab3_dynv",
                                                     HTML("Download additional plots"),
                                                     style="background-color: #d6fffe"))
                                 
                             ) # ends conditional
                         ), # ends fluid row
                     ), # ends conditional
                     
                     bookmarkButton(style = "position:absolute;right:2em; background-color:#BF3EFF; color:#FFFFFF;"),
                     
                     tags$hr(),
                     p(strong("Asc-Seurat, version 1.0"), "- Released on March 19th, 2021.", align = "center")
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
                                                div(class = "down-group",
                                                    numericInput("goplotwidth", "Width (cm)", 10, min = 1, step = 1),
                                                    numericInput("goplotheight", "Height (cm)", 10, min = 1, step = 1)
                                                )),
                                         column(3,
                                                div(class = "down-group",
                                                    selectInput("goplotdpi","Resolution (DPI)",
                                                                choices=c("100","200","300","400","500","600", "700", "800", "900", "1000"),
                                                                selected="300"),
                                                    selectInput("goplotdevice",
                                                                "File type",
                                                                choices=list("png"  = "png",
                                                                             "tiff" = "tiff",
                                                                             "jpeg" = "jpeg",
                                                                             "pdf"  = "pdf",
                                                                             "svg"  = "svg"),
                                                                selected="svg",
                                                                multiple=FALSE,
                                                                selectize=TRUE)
                                                )),
                                         downloadButton("goplot_down", HTML("Download Plot"))
                                     )
                                 ),
                                 HTML(paste("<br/>"))
                                 
                             ),
                             
                             bookmarkButton(style = "position:absolute;right:2em; background-color:#BF3EFF; color:#FFFFFF;"),
                             
                             tags$hr(),
                             p(strong("Asc-Seurat, version 1.0"), "- Released on March 19th, 2021.", align = "center")
                         )
                     )
                     
            )
            
        ), ## Closing the ui
        width = 12)
}
