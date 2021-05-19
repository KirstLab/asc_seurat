read_file_markers <- function(id,
                              markers_list_file = markers_list_file,
                              feed_ID = feed_ID,
                              ext = ext,
                              header_opt = header_opt) {

    validate(need(header_opt != "", "Please, inform if the file has a header."))
    
    test <- ext %in% c(
        'text/csv',
        'text/comma-separated-values',
        'text/tab-separated-values',
        'csv',
        'tsv')

    shinyFeedback::feedbackDanger(feed_ID, test == F, "Format is not supported! Upload a CSV or TSV file.")

    if ( header_opt == "Yes" ) {
    
    features <- switch(ext,
           csv = vroom::vroom(markers_list_file, delim = ","),
           tsv = vroom::vroom(markers_list_file, delim = "\t"),
           validate("Invalid file; Please upload a .csv or .tsv file")
           )
    
    } else if ( header_opt == "No" ) {
        
        features <- switch(ext,
                           csv = vroom::vroom(markers_list_file, delim = ",", col_names = F),
                           tsv = vroom::vroom(markers_list_file, delim = "\t", col_names = F),
                           validate("Invalid file; Please upload a .csv or .tsv file")
        )  
        
    }
    

    if ( ncol(features) >= 2 ) {
        shinyFeedback::feedbackSuccess(feed_ID, ncol(features) >= 2, "")
    } else {
        shinyFeedback::feedbackDanger(feed_ID, ncol(features) < 2, "The input file must have at least two columns!")
        validate( need( ncol(features) >= 2, "The input file must have at least two columns", feed_ID) )
    }

    if ( ncol(features) == 2 ) {

        colnames(features)[1:2] <- c("GeneID", "Group")

    } else if (ncol(features) > 2) {

        colnames(features)[1:3] <- c("GeneID", "Group", "Name")

    }

    features

}

pickerInput_markers_group <- function(id,
                                      genes = genes) {

    groups_names <- unique(genes$Group)

    div(class = "option-group",
        shinyWidgets::pickerInput(
            inputId = id,
            label = "Select the group of markers to test",
            choices = sort( as.character(groups_names) ),
            multiple = TRUE,
            options = list(`actions-box` = TRUE),
            selected = sort( as.character(groups_names) )[1]
        ))

}

pickerInput_markers_genes <- function(id,
                                      genes = genes,
                                      features_group = features_group,
                                      id_choice = id_choice) {

    groups_names <- genes
    genes_names <- dplyr::filter(groups_names,
                                 Group %in% features_group)

    if (id_choice == "ID") {

        genes_names <- unique(genes_names$GeneID)

    } else if (id_choice == "name") {

        genes_names <- unique(genes_names$Name)

    }

    genes_names <- genes_names[!is.na(genes_names)]

    div(class = "option-group",
        shinyWidgets::pickerInput(
            inputId = id,
            label = "Select the genes to show",
            choices = sort(as.character(genes_names)),
            multiple = TRUE,
            options = list(`actions-box` = TRUE)
        ))

}