my_withSpinner <- function(ui_object_name = ui_object_name, hide_ui = FALSE) {
    
    shinycssloaders::withSpinner(
        ui_object_name,
        type = 4,
        color = "#6504b0",
        size = 0.75,
        hide.ui = hide_ui
    )
    
}

textInput_regex <- function(id, value = " ") {
    
    div(class = "option-group",
        textInput(id,
                  label = "Common identifier of mitochondrial genes",
                  value = value)
    )
    
}

numericInput_var_genes <- function(id,
                                   label = "Number of variable genes for PCA") {
    div(class = "option-group",
        numericInput(id,
                     label = label,
                     value = 3000)
    )
    
}

numericInput_n_of_PCs <- function(id,
                                  label = "Select the number of components",
                                  value = 30) {
    
    div(class = "option-group",
        numericInput(id,
                     label = label,
                     value = value)
    )
    
}

radioInput_most_var_method <- function(id) {
    
    div(class = "option-group",
        radioButtons(id,
                     "Select the method to detect the most variable genes",
                     choices = list("vst" = "vst",
                                    "mean.var.plot (mvp)" = "mvp",
                                    "dispersion (disp)" = "disp"),
                     selected = c("vst"))
        
    )
    
}

numericInput_scale_factor <- function(id) {
    
    div(class = "option-group",
        numericInput(id,
                     label = "Scale factor",
                     value = 10000)
    )
    
}

numericInput_max_mito_perc <- function(id, value = 2) {
    
    div(class = "option-group",
        numericInput(id,
                     label = "Exclude any cell with more than this percentage of transcripts belonging to the mitochondrial genes",
                     value = value)
    )
    
}

numericInput_plot_height <- function(id, value= 15) {
    
    div(class = "down-group",
        numericInput(id, "Height (cm)", step=0.5, value=value)
    )
    
}

numericInput_plot_width  <- function(id, value = 15) {
    
    div(class = "down-group",
        numericInput(id, "Width (cm)", step=0.5, value=value)
    )
    
}

selectInput_plot_res <- function(id) {
    
    div(class = "down-group",
        selectInput(id,
                    "Resolution (DPI)",
                    choices=c("100","200","300","400","500","600"),
                    selected="100")
    )
    
}

selectInput_plot_format <- function(id) {
    
    div(class = "down-group",
        selectInput(id,
                    "File type",
                    choices=list("png" = "png",
                                 "tiff" = "tiff",
                                 "jpeg" = "jpeg",
                                 "pdf" = "pdf",
                                 "svg" = "svg"),
                    selected="png",
                    multiple=FALSE,
                    selectize=TRUE)
    )
    
}

actionButtonInput <- function(id, label = label) {
    
    actionButton(id,
                 label = label,
                 style="background-color: #d6fffe",
                 class = "btn-block")
    
}

numericInput_resolution_clust <- function(id) {
    
    div(class = "option-group",
        numericInput(id,
                     label = "Select the resolution for clustering",
                     value = 0.6,
                     step = 0.1)
    )
}

radioButtons_filter_clusters  <- function(id) {
    
    div(class = "option-group",
        radioButtons(id,
                     "Do you want to select or exclude clusters of cells and reanalyze the data?",
                     choices = list("Yes" = 1,
                                    "No" = 0),
                     selected = 0)
    )
    
}

radioButtons_filter_clusters_opt <- function(id) {
    
    div(class = "option-group",
        radioButtons(id,
                     "Do you want to select or exclude the clusters?",
                     choices = list("Select" = 0,
                                    "Exclude" = 1),
                     selected = 0)
    )
    
}

radioButtons_find_markers <- function(id) {
    
    div(class = "option-group",
        radioButtons(id,
                     "Do you want to identify markers or D.E. genes?",
                     choices = list("Yes" = 1,
                                    "No" = 0),
                     selected = c(0))
        
    )
    
}

pickerInput_find_markers_test <- function(id) {
    
    shinyWidgets::pickerInput(
        inputId = id,
        label = "Select the statistical test",
        choices = c("wilcox", "bimod", "roc", "t"),
        selected = "wilcox",
        multiple = FALSE,
        options = list(`actions-box` = TRUE)
    )
    
}

fileInput_markers_list <- function(id) {
    
    div(class = "option-group",
        fileInput(id,
                  label = "Input the list of markers",
                  accept = c("text/csv",
                             "text/comma-separated-values,",
                             ".csv",
                             "text/tsv",
                             "text/tab-separated-values,",
                             ".tsv"
                  ))
    )
    
}

radioButtons_slot_selection_heatmap <- function(id) {
    
    div(class = "option-group",
        radioButtons(id,
                     "Select the expression values to show",
                     choices = list("counts" = "counts",
                                    "normalized" = "data",
                                    "normalized and scaled" = "scale.data"),
                     selected = "scale.data"))
    
}

radioButtons_slot_selection_feature_plot <- function(id) {
    
    div(class = "option-group",
        radioButtons(id,
                     "Select the expression values to show",
                     choices = list("counts" = "counts",
                                    "normalized (recommended)" = "data"),
                     selected = c("data")))
    
}

radioButtons_down_add_plots <- function(id) {
    
    div(class = "option-group",
        radioButtons(id,
                     "Do you want to download the plots?",
                     choices = list("Yes" = 1,
                                    "No" = 0),
                     selected = 0))
    
}

define_what_id_to_use <- function(id) {
    
    div(class = "option-group",
        radioButtons(id,
                     "ID options",
                     choices = list("Use IDs" = "ID",
                                    "Use Name" = "name"),
                     selected = "ID"))
    
}

define_if_use_all_genes_or_select <- function(id) {
    
    div(class = "option-group",
        radioButtons(id,
                     "Visualization options",
                     choices = list("Show all genes" = 1,
                                    "Select genes to show" = 0),
                     selected = 1))
    
}

select_norm_methods <- function(id) {
    
    div(class = "option-group",
        radioButtons(id,
                     "Select the normalization method",
                     choices = list("LogNormalize" = 0,
                                    "SCTransform"  = 1),
                     selected = 0))
    
}

input_project_name <-  function(id){
    
    div(class = "option-group",
        textInput(id,
                  label = "Project name",
                  value = ""))
    
}