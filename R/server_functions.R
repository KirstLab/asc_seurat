my_reactable <- function(my_input) {

reactable(my_input,
          defaultColDef = colDef(align = "center"),
          rownames = F,
          filterable = TRUE,
          bordered = TRUE,
          highlight = TRUE,
          searchable = TRUE,
          showPageSizeOptions = TRUE,
          resizable = T,
          defaultPageSize = 10,
          pageSizeOptions = c(5, 10, 25, 50, 100, 200),
          showPagination = TRUE)

}

heatmap_genes <- function(id,
                          features = features,
                          sc_data_av = sc_data_av){

  features_selec <- as.data.frame(unique(features))

  if (nrow(features) == 1) {

    sc_data_av_feat <- sc_data_av[base::rownames(sc_data_av) %in% features_selec[, 1], ]
    sc_data_av_feat <- as.data.frame(t(sc_data_av_feat))

    size_hetmap <- as.numeric(nrow(sc_data_av_feat))

  } else {

    sc_data_av_feat <- sc_data_av[base::rownames(sc_data_av) %in% features_selec[, 1], ]

    size_hetmap <- as.numeric(nrow(sc_data_av_feat))

  }
}

feature_plots_gene_selection <- function(id,
                                         label = "Select the genes to feature plots",
                                         genes_names = genes_names,
                                         sc_data_av_react = sc_data_av_react,
                                         genes_ids = genes_ids) {

  genes_names <- genes_names[ genes_names$GeneID %in% base::rownames(sc_data_av_react), ]

  if (genes_ids == "ID") {

    genes_names <- unique(genes_names$GeneID)

  } else if (genes_ids == "name") {

    genes_names <- unique(genes_names$Name)

  }

  div(class = "option-group",
      shinyWidgets::pickerInput(
        inputId = id,
        label = label,
        choices = sort(as.character(genes_names)),
        multiple = TRUE,
        options = list(`actions-box` = TRUE)
      ))

}