# Asc-Seurat
# Version 1.3
set.seed(1407)

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

# Bioconductor
suppressMessages( require(ComplexHeatmap) )
suppressMessages( require(tradeSeq) )
suppressMessages( require(SingleCellExperiment) )
suppressMessages( require(slingshot) )
suppressMessages( require(multtest) )
suppressMessages( require(biomaRt) )
suppressMessages( require(topGO) )

# dynverse packages
suppressMessages( require(dynplot) )
suppressMessages( require(dynwrap) )
suppressMessages( require(dynfeature) )

## Help function to improve Seurat's feature plots
FeaturePlotSingle <- function(obj, feature, metadata_column, ...) {

    all_cells <- colnames(obj)
    groups <- unique(obj@meta.data[, metadata_column])

    # the minimal and maximal of the value to make the legend scale the same.
    minimal <- min(obj[['RNA']]@data[feature, ])
    maximal <- max(obj[['RNA']]@data[feature, ])

    ps <- list()

    for (group in groups) {

        subset_indx <- obj@meta.data[, metadata_column] == group
        subset_cells <- all_cells[subset_indx]

        p <- suppressMessages(FeaturePlot(obj, features = feature, cells = subset_cells, ...) +
                                  scale_colour_gradient2(limits=c(minimal, maximal),
                                                         midpoint = maximal / 2,
                                                         low = "gray80",
                                                         mid = "gold",
                                                         high = "red")) +

            ggtitle(group) +
            theme(plot.title = element_text(size = 14, face = "bold"),
            )

        ps[[group]]<- p
    }


    return(ps)
}

## Help function to run slingshot without the docker images
run_fun_slingshot <- function(expression, priors, verbose, seed) {

    start_id <- priors$start_id
    end_id <- priors$end_id
    dimred <- priors$dimred
    groups_id <- priors$groups_id

    #####################################
    ###        INFER TRAJECTORY       ###
    #####################################

    #   ____________________________________________________________________________
    #   Preprocessing                                                           ####

    start_cell <- if (!is.null(start_id)) { sample(start_id, 1) } else { NULL }

    # TIMING: done with preproc
    checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

    #   ____________________________________________________________________________
    #   Dimensionality reduction                                                ####
    # only do dimred if it is not yet given by prior information
    if (is.null(dimred)) {
        ndim <- parameters$ndim
        if (ncol(expression) <= ndim) {
            message(paste0(
                "ndim is ", ndim, " but number of dimensions is ", ncol(expression),
                ". Won't do dimensionality reduction."
            ))
            rd <- as.matrix(expression)
        } else {
            pca <- irlba::prcomp_irlba(expression, n = ndim)

            # select optimal number of dimensions if ndim is large enough
            if (ndim > 3) {
                # this code is adapted from the expermclust() function in TSCAN
                # the only difference is in how PCA is performed
                # (they specify scale. = TRUE and we leave it as FALSE)
                x <- 1:ndim
                optpoint1 <- which.min(sapply(2:10, function(i) {
                    x2 <- pmax(0, x - i)
                    sum(lm(pca$sdev[1:ndim] ~ x + x2)$residuals^2 * rep(1:2,each = 10))
                }))

                # this is a simple method for finding the "elbow" of a curve, from
                # https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve
                x <- cbind(1:ndim, pca$sdev[1:ndim])
                line <- x[c(1, nrow(x)),]
                proj <- princurve::project_to_curve(x, line)
                optpoint2 <- which.max(proj$dist_ind)-1

                # we will take more than 3 PCs only if both methods recommend it
                optpoint <- max(c(min(c(optpoint1, optpoint2)), 3))
            } else {
                optpoint <- ndim
            }

            rd <- pca$x[, seq_len(optpoint)]
            rownames(rd) <- rownames(expression)
        }
    } else {
        message("Using given dimred")
        rd <- dimred
    }


    #   ____________________________________________________________________________
    #   Clustering                                                              ####
    # only do clustering if it is not yet given by prior information
    if (is.null(groups_id)) {
        # max clusters equal to number of cells
        max_clusters <- min(nrow(expression)-1, 10)

        # select clustering
        if (parameters$cluster_method == "pam") {
            if (nrow(rd) > 10000) {
                warning("PAM (the default clustering method) does not scale well to a lot of cells. You might encounter memory issues. This can be resolved by using the CLARA clustering method, i.e. cluster_method = 'clara'.")
            }
            clusterings <- lapply(3:max_clusters, function(K){
                cluster::pam(rd, K) # we generally prefer PAM as a more robust alternative to k-means
            })
        } else if (parameters$cluster_method == "clara") {
            clusterings <- lapply(3:max_clusters, function(K){
                cluster::clara(rd, K) # we generally prefer PAM as a more robust alternative to k-means
            })
        }

        # take one more than the optimal number of clusters based on average silhouette width
        # (max of 10; the extra cluster improves flexibility when learning the topology,
        # silhouette width tends to pick too few clusters, otherwise)
        wh.cl <- which.max(sapply(clusterings, function(x){ x$silinfo$avg.width })) + 1
        labels <- clusterings[[min(c(wh.cl, 8))]]$clustering
    } else {
        message("Using given groups/clustering")
        labels <- groups_id %>% deframe()
    }

    start.clus <-
        if(!is.null(start_cell)) {
            labels[[start_cell]]
        } else {
            NULL
        }
    end.clus <-
        if(!is.null(end_id)) {
            unique(labels[end_id])
        } else {
            NULL
        }

    #   ____________________________________________________________________________
    #   Infer trajectory                                                        ####
    sds <- slingshot::slingshot(
        rd,
        clusterLabels = labels,
        start.clus = start.clus,
        end.clus = end.clus,
        # shrink = parameters$shrink,
        # reweight = parameters$reweight,
        # reassign = parameters$reassign,
        # thresh = parameters$thresh,
        # maxit = parameters$maxit,
        # stretch = parameters$stretch,
        # smoother = parameters$smoother,
        # shrink.method = parameters$shrink.method
        reducedDim = 'PCA',
        #clusterLabels = 'seurat_clusters',
        shrink = 1L,
        reweight = TRUE,
        reassign = TRUE,
        maxit = 10L,
        smoother = "smooth.spline",
    )

    start_cell <- apply(slingshot::slingPseudotime(sds), 1, min) %>% sort() %>% head(1) %>% names()
    start.clus <- labels[[start_cell]]

    # TIMING: done with method
    checkpoints$method_aftermethod <- as.numeric(Sys.time())

    #   ____________________________________________________________________________
    #   Create output                                                           ####

    # satisfy r cmd check
    from <- to <- NULL

    # collect milestone network
    lineages <- slingLineages(sds)
    lineage_ctrl <- slingParams(sds)
    cluster_network <- lineages %>%
        map_df(~ tibble(from = .[-length(.)], to = .[-1])) %>%
        unique() %>%
        mutate(
            length = lineage_ctrl$dist[cbind(from, to)],
            directed = TRUE
        )

    # collect dimred
    dimred <- reducedDim(sds)

    # collect clusters
    cluster <- slingClusterLabels(sds)

    # collect progressions
    adj <- slingAdjacency(sds)
    lin_assign <- apply(slingCurveWeights(sds), 1, which.max)

    progressions <- map_df(seq_along(lineages), function(l) {
        ind <- lin_assign == l
        lin <- lineages[[l]]
        pst.full <- slingPseudotime(sds, na = FALSE)[,l]
        pst <- pst.full[ind]
        means <- sapply(lin, function(clID){
            stats::weighted.mean(pst.full, cluster[,clID])
        })
        non_ends <- means[-c(1,length(means))]
        edgeID.l <- as.numeric(cut(pst, breaks = c(-Inf, non_ends, Inf)))
        from.l <- lineages[[l]][edgeID.l]
        to.l <- lineages[[l]][edgeID.l + 1]
        m.from <- means[from.l]
        m.to <- means[to.l]

        pct <- (pst - m.from) / (m.to - m.from)
        pct[pct < 0] <- 0
        pct[pct > 1] <- 1

        tibble(cell_id = names(which(ind)), from = from.l, to = to.l, percentage = pct)
    })

    #   ____________________________________________________________________________
    #   Save output                                                             ####
    output <-
        dynwrap::wrap_data(
            cell_ids = rownames(expression)
        ) %>%
        dynwrap::add_trajectory(
            milestone_network = cluster_network,
            progressions = progressions,
            lineages = lineages
        ) %>%
        dynwrap::add_dimred(
            dimred = dimred
        ) %>%
        dynwrap::add_timings(checkpoints)
}

#######################################
### Function to retry an expression ###
#######################################
retry <- function(expr, isError=function(x) "try-error" %in% class(x), maxErrors=5, sleep=0) {
    attempts = 0
    retval = try(eval(expr))
    while (isError(retval)) {
        attempts = attempts + 1
        if (attempts >= maxErrors) {
            msg = sprintf("retry: too many retries [[%s]]", capture.output(str(retval)))
            #flog.fatal(msg)
            stop(msg)
        } else {
            msg = sprintf("retry: error in attempt %i/%i [[%s]]", attempts, maxErrors, 
                          capture.output(str(retval)))
            #flog.error(msg)
            warning(msg)
        }
        if (sleep > 0) Sys.sleep(sleep)
        retval = try(eval(expr))
    }
    return(retval)
}

###########################
### Load phytozome mart ###
###########################
phytozome_mart <- new("Mart",
                      biomart = "phytozome_mart",
                      vschema = "zome_mart",
                      host    = "https://phytozome.jgi.doe.gov:443/biomart/martservice")

######################################
### Function to TEST the attribute ###
### for IDs in the GO enrichment   ###
######################################
att_has_id <- function(selection, universe) {
    if (any(selection %in% universe)) {
        NULL
    } else {
        paste0("PROBLEM! -- The gene identification (first column of csv) of gene universe does not match",
               " those of the input target genes. Please, re-format your data. Both CSVs must have the",
               " same gene identification!")
    }
}

#############################
### Function to run topGO ###
#############################
runTopGO <- function(gocat, allgenes, topgouniverse) {
    GOdata <- new(
        "topGOdata",
        ontology = as.character(gocat),
        allGenes = allgenes,
        annot = annFUN.gene2GO,
        gene2GO = topgouniverse
    )
    
    # Calculate
    #sigGenes(GOdata) # show significant genes
    resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    resultFis
    
    # test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
    # resultKS <- getSigGroups(GOdata, test.stat)
    resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
    resultKS
    
    # test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
    # resultWeight <- getSigGroups(GOdata, test.stat)
    resultWeight <- runTest(GOdata, statistic = "fisher")
    resultWeight
    
    # With different statistics
    allRes <- GenTable(
        GOdata,
        classicFisher = resultFis,
        KS            = resultKS,
        weightFisher  = resultWeight,
        orderBy       = "weightFisher",
        ranksOf       = "classicFisher"
    )
    
    # add macro GO
    if (!is.null(dim(allRes)[1])) {
        if(as.character(gocat) == 'BP') {
            macro <- 'Biological Process'
        }
        if(as.character(gocat) == 'CC') {
            macro <- 'Cellular Component'
        }
        if(as.character(gocat) == 'MF') {
            macro <- 'Molecular Function'
        }
        allRes$Macro <- macro
    }
    
    return(allRes)
}

##############################
### Function to plot topGO ###
##############################
plotTopGO <- function(ggdata, test_selected) {
    
    if (test_selected == "KS") {
        
        gg_out <- ggplot(ggdata,
                      aes(x = Term,
                          y = -log10(KS),
                          size = 3
                      ))
        
    } else {
        
        gg_out <- ggplot(ggdata,
                      aes(x = Term,
                          y = -log10(weightFisher),
                          size = 3
                      ))
        
    }
    
    gg_out <- gg_out + 
        expand_limits(y = 1) +
        geom_point(shape = 21, fill = "grey50") +
        scale_size(range = c(2.5, 12.5)) +
        xlab('') + ylab( paste0('Enrichment score -log10(', test_selected, ')') ) +
        labs(
            title = 'GO Enrichment Analysis',
            caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001') +
        facet_grid(rows = vars(Macro), drop = T, scales = 'free') +
        geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
                   linetype = rep(c("dotted", "longdash", "solid"), 3),
                   colour = rep(c("black", "black", "black"), 3),
                   size = rep(c(0.5, 1.5, 3), 3)) +
        
        theme_bw(base_size = 24) +
        theme(
            legend.position = "none",
            plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
            plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
            plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1, hjust = 0.5),
            strip.text = element_text(angle = 0, size = 10, face = 'bold', vjust = 1),
            
            axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
            axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
            axis.title = element_text(size = 12, face = 'bold'),
            axis.title.x = element_text(size = 12, face = 'bold'),
            axis.title.y = element_text(size = 12, face = 'bold'),
            axis.line = element_line(colour = 'black')) +
        coord_flip()
    
    return(gg_out)
}

function(input, output, session) {

    ##############################################
    #### Properly load user data using docker ####
    ##############################################
    if (dir.exists('/app/user_work')) {
        setwd('/app/user_work')
    }

    #####################################
    ######   Tab 1 - Clustering    ######
    #####################################

    output$min_count_ui <- renderUI({

        numericInput("min_count",
                     label = "Keep only cells that expressed at least this number of genes",
                     value = input$min_features)
    })

    output$max_count_ui <- renderUI({

        max_value <- input$min_features + 2000

        numericInput("max_count",
                     label = "Exclude any cell that expressed more than this number of genes (i.e. possible doublets)",
                     value = max_value)
    })

    output$select_sample_tab1 = renderUI({

        dir_list <- list.dirs('./data', recursive=FALSE)

        div(class = "option-group",
            pickerInput(
                inputId = "sample_folder_tab1",
                label = "Select the sample to use",
                choices = sort(dir_list),
                multiple = FALSE,
                options = list(`actions-box` = TRUE)
            ))
    })

    sing_cell_data.data <- eventReactive(input$load_10X, {

        showNotification("Loading the data",
                         duration = NULL,
                         id = "p1")

        Read10X(data.dir = req(input$sample_folder_tab1))

    })

    # Creates the seurat object considering the criteria set by the user
    single_cell_data_reac <- reactive({

        # Initialize the Seurat object with the raw (non-normalized data).
        sing_cell_data <- CreateSeuratObject(counts = sing_cell_data.data(),
                                             project = input$proj_name,
                                             min.cells = input$min_cells,
                                             min.features = input$min_features)

        # Calculate the % of mithocondrial contamination
        sing_cell_data <- PercentageFeatureSet(sing_cell_data,
                                               pattern = input$mito_regex,
                                               col.name = "percent.mt")

        return(sing_cell_data)

    })

    # Visualize QC metrics as a violin plot
    output$VlnPlot <- renderPlot({

        data_set <- single_cell_data_reac()
        removeNotification(id = "p1")

        VlnPlot(data_set,
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                ncol = 3,
                split.plot = F)

    })

    # Filtering and removing cells based on counts and % of mitochondrial
    single_cell_data_filt <- reactive({
        data_sc <- single_cell_data_reac()

        # Filtering features and cells based on the counts and % of mito contamination.
        subset(data_sc,
               subset = nFeature_RNA > input$min_count &
                   nFeature_RNA < input$max_count &
                   percent.mt < input$max_mito_perc)

    })

    observeEvent(input$run_vinplot, {

        output$VlnPlot_filt <- renderPlot({

            VlnPlot(single_cell_data_filt(),
                    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    ncol = 3,
                    split.plot = F)

        })

    })

    output$p1_down <- downloadHandler(

        filename = function() {
            paste(input$dataset, ".", input$p1_format, sep = "")
        },
        content = function(file) {

            withProgress(message = "Please wait, preparing the data for download.",
                         value = 0.5, {

                             height <- as.numeric(input$p1_height)
                             width <- as.numeric(input$p1_width)
                             res <- as.numeric(input$p1_res)


                             p <- VlnPlot(single_cell_data_filt(),
                                          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                                          ncol = 3,
                                          split.plot = F)

                             ggsave(file,
                                    p,
                                    height=height,
                                    width=width,
                                    units="cm",
                                    dpi=res)

                         })

        }

    )

    # # If user select SCTransform
    # single_cell_data_SCT <- eventReactive(input$run_pca, {
    #
    #     showNotification("Normalizing the data (SCTransform)",
    #                      duration = NULL,
    #                      id = "m1")
    #
    #     sc_data <- single_cell_data_filt()
    #
    #     ret_data <-  suppressWarnings(SCTransform(sc_data,
    #                                               vars.to.regress = "percent.mt",
    #                                               verbose = FALSE))
    #
    #     removeNotification(id = "m1")
    #     ret_data
    # })
    #
    # single_cell_data_pca_SCT <- reactive({
    #     showNotification("Running PCA",
    #                      duration = NULL,
    #                      id = "m2")
    #
    #     ret_data <- RunPCA(single_cell_data_SCT(), verbose = FALSE)
    #
    #     removeNotification(id = "m2")
    #     ret_data
    # })


    # If user select LogNormalize
    single_cell_data_norm <- eventReactive(input$run_pca, {
        showNotification("Normalizing the data (LogNormalize)",
                         duration = NULL,
                         id = "m3")

        data_sc <- single_cell_data_filt()
        ret_data <- NormalizeData(data_sc,
                                  normalization.method = "LogNormalize",
                                  scale.factor = input$scale_factor)

        removeNotification(id = "m3")
        ret_data
    })

    clusters_single_cell_data_reso_umap <- reactive({

        sc_data <- single_cell_data_reso_umap()
        as.numeric( unique( as.character(sc_data@meta.data$seurat_clusters ) ) )


    })

    output$cluster_list_ui <- renderUI ({

        clusters <- clusters_single_cell_data_reso_umap()

        pickerInput(
            inputId = "cluster_list",
            label = "Choose clusters to select or exclude",
            choices = sort(clusters),
            multiple = TRUE,
            options = list(`actions-box` = TRUE)
        )

    })

    to_filter <- reactive({

        if ( input$filter_clusters_opt == "select" ) {

            to_filter <- subset(single_cell_data_reso_umap(),
                                idents = as.numeric(req(input$cluster_list))
            )

            to_filter_ch <- to_filter@meta.data
            to_filter_ch <- rownames(to_filter_ch)

            to_filter_ch

        } else if ( input$filter_clusters_opt == "exclude" ) {

            to_filter <- subset(single_cell_data_reso_umap(),
                                idents = as.numeric(req(input$cluster_list)),
                                invert = TRUE)

            to_filter_ch <- to_filter@meta.data
            to_filter_ch <- rownames(to_filter_ch)

            to_filter_ch
        }



    })

    output$range <- renderPrint({ to_filter() })

    single_cell_data_scaled <- eventReactive(input$run_pca, {

        showNotification("Scalling the data",
                         duration = NULL,
                         id = "m4")

        data_sc <- FindVariableFeatures(single_cell_data_norm(),
                                        selection.method = input$most_var_method,
                                        nfeatures = input$n_of_var_genes)

        allgenes <- rownames(data_sc)
        ret_data <- ScaleData(data_sc, features = allgenes)

        removeNotification(id = "m4")
        ret_data

    })

    single_cell_data_scaled_filtered <- eventReactive(c(input$run_pca, input$rerun_after_filtering), {

        showNotification("Scalling the data",
                         duration = NULL,
                         id = "m4")

        sc_data <- single_cell_data_norm()

        if ( input$filter_clusters == 1 ) {

            if (input$filter_clusters_opt == "select") {

                sc_data <- subset(sc_data,
                                  cells = to_filter())

            } else if (input$filter_clusters_opt == "exclude") {

                sc_data <- subset(sc_data,
                                  cells = to_filter())

            }

            data_sc <- FindVariableFeatures(sc_data,
                                            selection.method = input$most_var_method,
                                            nfeatures = input$n_of_var_genes)

            allgenes <- rownames(data_sc)
            ret_data <- ScaleData(data_sc, features = allgenes)

        } else if (input$filter_clusters == 0 ) {

            data_sc <- FindVariableFeatures(sc_data,
                                            selection.method = input$most_var_method,
                                            nfeatures = input$n_of_var_genes)

            allgenes <- rownames(data_sc)
            ret_data <- ScaleData(data_sc, features = allgenes)

        }


        removeNotification(id = "m4")
        ret_data
    })

    single_cell_data_pca <- eventReactive(c(input$run_pca, input$rerun_after_filtering), {

        if ( input$filter_clusters == 1 ) {

            data <- single_cell_data_scaled_filtered()

        } else if ( input$filter_clusters == 0 ) {

            data <- single_cell_data_scaled()

        } else {

            print("Something is wrong!")
        }

        showNotification("Running PCA",
                         duration = NULL,
                         id = "m5")
        #data_sc <- single_cell_data_scaled()

        ret_data <-  RunPCA(data,
                            features = VariableFeatures(object = data), verbose = F)

        removeNotification(id = "m5")
        ret_data
    })

    single_cell_data_pca2 <- eventReactive(c(input$run_pca, input$rerun_after_filtering), {

        if (input$normaliz_method == "SCTransform") {

            single_cell_data_pca2 <- single_cell_data_pca_SCT()

        } else if (input$normaliz_method == "LogNormalize") {

            single_cell_data_pca2 <- single_cell_data_pca()

        }

    })

    # Elbow plot (also triggers the PCA)
    observeEvent(c(input$run_pca, input$rerun_after_filtering), {

        # Generates the elbow plot showing the PCs #
        output$n_of_PCAs <- renderPlot({

            data_sc <- single_cell_data_pca2()

            ElbowPlot(data_sc, ndims = 50, reduction = "pca")

        })

    })

    output$p2_down <- downloadHandler(

        filename = function() {
            paste(input$dataset, ".", input$p2_format, sep = "")
        },
        content = function(file) {

            withProgress(message = "Please wait, preparing the data for download.",
                         value = 0.5, {

                             height <- as.numeric(input$p2_height)
                             width <- as.numeric(input$p2_width)
                             res <- as.numeric(input$p2_res)

                             p <- ElbowPlot(single_cell_data_pca2(), ndims = 50, reduction = "pca")

                             ggsave(file,
                                    p,
                                    height=height,
                                    width=width,
                                    units="cm",
                                    dpi=res)

                         })
        }

    )


    # to_filter <- reactive({
    #
    #
    #     if ( input$filter_clusters_opt == "select" ) {
    #
    #         to_filter <- subset(single_cell_data_reso_umap(),
    #                             #idents = c(0, 1)
    #                             idents = as.numeric(strsplit(input$cluster_list, ",")[[1]])
    #         )
    #
    #         colnames(to_filter@meta.data)
    #
    #     } else if ( input$filter_clusters_opt == "exclude" ) {
    #
    #         to_filter <- subset(single_cell_data_reso_umap(),
    #                             idents = as.numeric(strsplit(input$cluster_list, ",")[[1]]),
    #                             #idents = c(0, 1),
    #                             invert = TRUE)
    #
    #         colnames(to_filter@meta.data)
    #
    #     }
    #
    # })


    single_cell_data_neigh <- eventReactive(input$run_clustering, {
        showNotification("Running the clustering step",
                         duration = NULL,
                         id = "m6")

        data_sc <- single_cell_data_pca2()

        FindNeighbors(data_sc, dims = 1:input$n_of_PCs)

    })

    single_cell_data_reso <- eventReactive(input$run_clustering, {

        FindClusters(single_cell_data_neigh(), resolution = input$resolution_clust)

    })

    single_cell_data_reso_umap <- eventReactive(input$run_clustering, {

        sc_data <- RunUMAP(single_cell_data_reso(), dims = 1:input$n_of_PCs)
        sc_data <- RunTSNE(sc_data, dims = 1:input$n_of_PCs)

    })

    observeEvent(input$run_clustering, {

        output$tSNE <- renderPlot({

            DimPlot(single_cell_data_reso_umap(), reduction = "tsne", label = T, pt.size = .1)

        })

        output$umap <- renderPlot({

            DimPlot(single_cell_data_reso_umap(), reduction = "umap", label = T, pt.size = .1)

        })

        output$cluster_size <-  renderPrint({

            sing_cell_data <- single_cell_data_reso_umap()

            sc_meta <- as.data.frame(sing_cell_data[[]])
            sc_meta$cellcluster <- rownames(sc_meta)
            sc_meta <- sc_meta[, c("cellcluster", "seurat_clusters")]

            removeNotification(id = "m6")

            table(sc_meta$seurat_clusters)
        })

    })

    output$p3_down <- downloadHandler(

        filename = function() {
            paste(input$dataset, ".", input$p3_format, sep = "")
        },
        content = function(file) {

            withProgress(message = "Please wait, preparing the data for download.",
                         value = 0.5, {

                             height <- as.numeric(input$p3_height)
                             width <- as.numeric(input$p3_width)
                             res <- as.numeric(input$p3_res)

                             if ( input$p3_down_opt == "UMAP" ) {

                                 p <- DimPlot(single_cell_data_reso_umap(), reduction = "umap", label = T, pt.size = .1)

                             } else if ( input$p3_down_opt == "t-SNE") {

                                 p <- DimPlot(single_cell_data_reso_umap(), reduction = "tsne", label = T, pt.size = .1)

                             }

                             ggsave(file,
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
            paste(input$dataset, ".rds", sep = "")
        },
        content = function(file) {

            withProgress(message = "Please wait, preparing the data for download.",
                         value = 0.5, {

                             saveRDS(single_cell_data_reso_umap(), file)

                         })

        }

    )

    output$find_markers_clust_id_tab1_ui <- renderUI ({

        clusters <- clusters_single_cell_data_reso_umap()

        pickerInput(
            inputId = "find_markers_clust_id_tab1",
            label = "Select the cluster of interest",
            choices = sort(clusters),
            multiple = FALSE,
            options = list(`actions-box` = TRUE)
        )

    })

    output$find_markers_clust_ID1_tab1_ui <- renderUI ({

        clusters <- clusters_single_cell_data_reso_umap()

        pickerInput(
            inputId = "find_markers_clust_ID1_tab1",
            label = "Select the cluster of interest",
            choices = sort(clusters),
            multiple = FALSE,
            options = list(`actions-box` = TRUE)
        )

    })

    output$find_markers_clust_ID2_tab1_ui <- renderUI ({

        clusters <- clusters_single_cell_data_reso_umap()

        # Exclude the cluster already selected in the option 1, since the two cluster must be different

        clusters <- clusters[ !clusters == req(input$find_markers_clust_ID1_tab1)]

        pickerInput(
            inputId = "find_markers_clust_ID2_tab1",
            label = "Select the cluster(s) to compare",
            choices = sort(clusters),
            multiple = T,
            options = list(`actions-box` = TRUE)
        )

    })

    # Identification of markers (tab 1)
    markers_tab1 <- eventReactive( input$run_ident_markers_tab1, {

        showNotification("Identifing markers or D.E. genes",
                         duration = NULL,
                         id = "tab1_n1")

        sc_data <- single_cell_data_reso_umap()

        if ( input$find_markers_tab1_opt == 0 ) {

            if ( is.na(input$find_markers_tab1_return.thresh) ) {

                markers_tab1 <- FindAllMarkers(sc_data,
                                               logfc.threshold = input$find_markers_tab1_logfc.threshold,
                                               min.pct = input$find_markers_tab1_min.pct,
                                               test.use = input$find_markers_tab1_test.use
                )

            } else {

                markers_tab1 <- FindAllMarkers(sc_data,
                                               logfc.threshold = input$find_markers_tab1_logfc.threshold,
                                               min.pct = input$find_markers_tab1_min.pct,
                                               test.use = input$find_markers_tab1_test.use,
                                               return.thresh = input$find_markers_tab1_return.thresh
                )

            }

        } else if ( input$find_markers_tab1_opt == 1 ) {

            if ( is.na(input$find_markers_tab1_return.thresh) ) {

                markers_tab1 <- FindMarkers(sc_data,
                                            ident.1 = req(input$find_markers_clust_id_tab1),
                                            ident.2 = NULL,
                                            logfc.threshold = input$find_markers_tab1_logfc.threshold,
                                            min.pct = input$find_markers_tab1_min.pct,
                                            test.use = input$find_markers_tab1_test.use,
                                            only.pos = input$find_markers_tab1_filt_pos)

            } else {

                markers_tab1 <- FindMarkers(sc_data,
                                            ident.1 = req(input$find_markers_clust_id_tab1),
                                            ident.2 = NULL,
                                            logfc.threshold = input$find_markers_tab1_logfc.threshold,
                                            min.pct = input$find_markers_tab1_min.pct,
                                            test.use = input$find_markers_tab1_test.use,
                                            only.pos = input$find_markers_tab1_filt_pos,
                                            return.thresh = input$find_markers_tab1_return.thresh)

            }


        } else if ( input$find_markers_tab1_opt == 2 ) {

            if ( is.na(input$find_markers_tab1_return.thresh) ) {

                markers_tab1 <- FindMarkers(sc_data,
                                            ident.1 = req(input$find_markers_clust_ID1_tab1),
                                            ident.2 = req(input$find_markers_clust_ID2_tab1),
                                            logfc.threshold = input$find_markers_tab1_logfc.threshold,
                                            min.pct = input$find_markers_tab1_min.pct,
                                            test.use = input$find_markers_tab1_test.use,
                                            only.pos = input$find_markers_tab1_filt_pos,
                                            return.thresh = input$find_markers_tab1_return.thresh)

            } else {

                markers_tab1 <- FindMarkers(sc_data,
                                            ident.1 = req(input$find_markers_clust_ID1_tab1),
                                            ident.2 = req(input$find_markers_clust_ID2_tab1),
                                            logfc.threshold = input$find_markers_tab1_logfc.threshold,
                                            min.pct = input$find_markers_tab1_min.pct,
                                            test.use = input$find_markers_tab1_test.use,
                                            only.pos = input$find_markers_tab1_filt_pos,
                                            return.thresh = input$find_markers_tab1_return.thresh)
            }


        }

        markers_tab1$geneID <- rownames(markers_tab1)
        markers_tab1 <- markers_tab1[, c( ncol(markers_tab1), 1:( ncol(markers_tab1)-1 ) ) ]

        markers_tab1

    })

    output$markers_tab1_react <- renderReactable({

        markers_tab1 <- markers_tab1()

        removeNotification(id = "tab1_n1")

        # markers_tab1 <- markers_tab1 %>%
        #     dplyr::select( geneID, avg_log2FC, pct.1, pct.2, p_val_adj ) %>%
        #     dplyr::rename( geneID = "Gene ID",
        #             avg_log2FC = "Fold Change (log)",
        #             p_val_adj = "FDR (adj. p-value)")
        #
        reactable(markers_tab1,
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

    })

    output$download_markers_tab1 <- downloadHandler(

        filename = function() {
            paste(input$dataset, ".csv", sep = "")
        },
        content = function(file) {

            withProgress(message = "Please wait, preparing the data for download.",
                         value = 0.5, {

                             markers_tab1 <- markers_tab1()

                             # markers_tab1 <- markers_tab1 %>%
                             #     select( geneID, avg_logFC, pct.1, pct.2, p_val_adj ) %>%
                             #     rename( geneID = "Gene ID",
                             #             avg_logFC = "Fold Change (log)",
                             #             p_val_adj = "FDR (adj. p-value)")

                             write.csv(markers_tab1, file, row.names = FALSE)

                         })
        }

    )

    single_cell_data_heatmap <- eventReactive(input$run_heatmap, {

        # showNotification("Generating Heatmap",
        #                  duration = NULL,
        #                  id = "m7")
        #
        # if (input$normaliz_method == "SCTransform") {
        #
        #     data_sc <- single_cell_data_reso_umap()
        #
        #     data_sc <- NormalizeData(data_sc,
        #                              assay = "RNA",
        #                              normalization.method = "LogNormalize",
        #                              scale.factor = input$scale_factor)
        #
        #     data_sc <- FindVariableFeatures(data_sc,
        #                                     assay = "RNA",
        #                                     selection.method = input$most_var_method,
        #                                     nfeatures = input$n_of_var_genes)
        #
        #     allgenes <- rownames(data_sc)
        #
        #     ScaleData(data_sc,
        #               assay = "RNA",
        #               features = allgenes)
        #
        #     #    data_sc <- FindNeighbors(data_sc, dims = 1:input$n_of_PCs)
        #     #    data_sc <- FindClusters(data_sc, resolution = input$resolution_clust)
        #
        # } else if (input$normaliz_method == "LogNormalize") {
        #
        #     data_sc <- single_cell_data_reso_umap()
        #
        # }

        data_sc <- single_cell_data_reso_umap()
        data_sc

    })

    ## Gene expression
    features <- reactive({

        inFile <- input$markers_list

        if (is.null(inFile))
            return(NULL)

        features <- read.csv(inFile$datapath, header = F)

        if ( ncol(features) == 2 ) {

            features <- features %>%
                dplyr::rename(GeneID = V1,
                              Group = V2)

        } else if ( ncol(features) > 2 ) {

            features <- features %>%
                dplyr::rename(GeneID = V1,
                              Group = V2,
                              Name = V3)

        }
        #    sc_data_av <- AverageExpression(sing_cell_data,
        #                                    # features = features,
        #                                    slot = "scale.data", #'counts', 'data', and 'scale.data'

        features

    })

    # Load the file and offers the parameters for heatmap
    observeEvent(input$load_markers, {

        # Painel that will apper after loading the list of markers containing the filter options
        output$marker_group_selec = renderUI({

            groups_names <- features()
            groups_names <- unique(groups_names$Group)

            div(class = "option-group",
                pickerInput(
                    inputId = "features_group",
                    label = "Select the group of markers to test",
                    choices = sort( as.character(groups_names) ),
                    multiple = TRUE,
                    options = list(`actions-box` = TRUE)
                ))


        })

        output$marker_filter_genes_q = renderUI({

            div(class = "option-group",
                radioButtons("filter_genes_q",
                             "Visualization options",
                             choices = list("Show all genes" = 1,
                                            "Select genes to show" = 0),
                             selected = 1))
        })

        # Define if using names of v5 or common names
        output$marker_genes_ids = renderUI({

            div(class = "option-group",
                radioButtons("genes_ids",
                             "ID options",
                             choices = list("Use IDs" = "v5",
                                            "Use Name" = "name"),
                             selected = "v5"))
        })

        output$marker_genes_selec = renderUI({

            groups_names <- features()
            genes_names <- dplyr::filter(groups_names,
                                         Group %in% req(input$features_group))

            if (input$genes_ids == "v5") {

                genes_names <- unique(genes_names$GeneID)

            } else if (input$genes_ids == "name") {

                genes_names <- unique(genes_names$Name)

            }

            genes_names <- genes_names[!is.na(genes_names)]

            # div(class = "option-group",
            #     selectizeInput('selected_genes',
            #                    'Select the genes to show',
            #                    choices = sort(genes_names),
            #                    selected = NULL,
            #                    multiple = T))

            div(class = "option-group",
                pickerInput(
                    inputId = "selected_genes",
                    label = "Select the genes to show",
                    choices = sort(as.character(genes_names)),
                    multiple = TRUE,
                    options = list(`actions-box` = TRUE)
                ))

        })

    })

    # Filter accordly with the parameters selected by the user
    filt_features <- eventReactive(input$run_heatmap, {

        features_f <- features()

        features_f <- dplyr::filter(features_f,
                                    Group %in% req(input$features_group))

        if (input$filter_genes_q == 0) {

            if (input$genes_ids == "v5") {

                features_f <- features_f[features_f$GeneID %in% req(input$selected_genes), ]

            } else if (input$genes_ids == "name") {

                features_f <- features_f[features_f$Name %in% req(input$selected_genes), ]

            }

        }

        features_f

    })

    sc_data_av_react <- eventReactive(input$run_heatmap, {

        # features <- filt_features()

        # get the expression values
        sc_data_av <- AverageExpression(single_cell_data_heatmap(),
                                        assays = "RNA",
                                        # features = features,
                                        slot = input$slot_selection_heatmap)


        sc_data_av <- as.matrix(sc_data_av[[1]])

        removeNotification(id = "m7")
        sc_data_av

    })

    observeEvent(input$run_heatmap, {

        # This calculate the height of the plot based on the n of genes so the plot can be ajusted o fit all genes
        heatmap_n_genes <- reactive({

            features <- filt_features()
            #    features_selec <- paste0("gene:", unique(features$GeneID))
            features_selec <- as.data.frame(unique(features$GeneID))

            # get the expression values
            sc_data_av <- sc_data_av_react()

            if (nrow(features) == 1) {

                sc_data_av_feat <- sc_data_av[rownames(sc_data_av) %in% features_selec[, 1], ]
                sc_data_av_feat <- as.data.frame(t(sc_data_av_feat))

                size_hetmap <- as.numeric(nrow(sc_data_av_feat))

            } else {

                sc_data_av_feat <- sc_data_av[rownames(sc_data_av) %in% features_selec[, 1], ]

                size_hetmap <- as.numeric(nrow(sc_data_av_feat))

            }
        })
        heatmap_Height <- reactive( 150 + ( 20 * heatmap_n_genes() ) )


        heat_map_prep <- reactive({

            features <- filt_features()

            # get the expression values
            sc_data_av <- sc_data_av_react()

            features_selec <- as.data.frame(unique(features$GeneID))

            if (nrow(features) == 1) {

                sc_data_av_feat <- sc_data_av[rownames(sc_data_av) %in% features_selec[, 1], ]
                sc_data_av_feat <- as.matrix(t(sc_data_av_feat))
                rownames(sc_data_av_feat) <- features_selec[1, 1]

            } else {

                sc_data_av_feat <- sc_data_av[rownames(sc_data_av) %in% features_selec[, 1], ]

            }
            sc_data_av_feat
        })

        # Generates the heatmap plot
        output$heat_map <- renderPlot({

            ### Drawing heat map ###

            heat_map_prep <- heat_map_prep()

            Heatmap(heat_map_prep, border = TRUE,
                    rect_gp = gpar(col = "white", lwd = 2),
                    column_title = "Clusters",
                    column_title_side = "bottom",
                    name = "Expression",
                    show_row_dend = T)


        }, height = heatmap_Height())

        output$heat_map_ui <- renderUI({

            plotOutput("heat_map", height = heatmap_Height())

        })

        output$marker_to_feature_plot = renderUI({

            sc_data_av_react <- sc_data_av_react()
            genes_names <- filt_features()

            #    genes_names <- genes_names[paste0("gene:", genes_names$GeneID) %in% rownames(sc_data_av_react), ]
            genes_names <- genes_names[ genes_names$GeneID %in% rownames(sc_data_av_react), ]

            if (input$genes_ids == "v5") {

                genes_names <- unique(genes_names$GeneID)

            } else if (input$genes_ids == "name") {

                genes_names <- unique(genes_names$Name)

            }

            # div(class = "option-group",
            #     selectizeInput('selected_genes_for_feature_plot',
            #                    'Select the genes to feature plots (must be in the heatmap)',
            #                    choices = genes_names,
            #                    selected = NULL,
            #                    multiple = T))

            div(class = "option-group",
                pickerInput(
                    inputId = "selected_genes_for_feature_plot",
                    label = "Select the genes to feature plots",
                    choices = sort(as.character(genes_names)),
                    multiple = TRUE,
                    options = list(`actions-box` = TRUE)
                ))

        })

        output$p4_down <- downloadHandler(

            filename = function() {
                paste(input$dataset, ".", input$p4_format, sep = "")
            },
            content = function(file) {

                withProgress(message = "Please wait, preparing the data for download.",
                             value = 0.5, {

                                 height <- as.numeric(input$p4_height)
                                 width <- as.numeric(input$p4_width)
                                 res <- as.numeric(input$p4_res)

                                 ### Drawing heat map ###

                                 p <- Heatmap(heat_map_prep(), border = TRUE,
                                              rect_gp = gpar(col = "white", lwd = 2),
                                              column_title = "Clusters",
                                              column_title_side = "bottom",
                                              name = "Expression",
                                              show_row_dend = T)

                                 # Complex heatmap does not work well with ggsave. So, using grid.grap to make it compatible
                                 gb = grid.grabExpr(draw(p))

                                 ggsave(file,
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

        features <- filt_features()

        if (input$genes_ids == "v5") {

            features_f <- dplyr::filter(features,
                                        GeneID %in% req(input$selected_genes_for_feature_plot))

        }
        else if (input$genes_ids == "name") {

            features_f <- dplyr::filter(features,
                                        Name %in% req(input$selected_genes_for_feature_plot))

        }

        #    paste0("gene:", unique(features_f$GeneID))
        as.character(unique(features_f$GeneID))

    })

    # output$test <- renderPrint({
    #
    #     features_selec()
    #
    # })

    observeEvent(input$run_feature_plot, {

        ### FEATURE PLOTS ####

        output$feature_plot <- renderUI({

            feat_length <- length(features_selec())

            plot_output_list <- lapply(1:feat_length, function(i) {
                plotname <- paste("plot1", i, sep="")
                plotOutput(plotname, height = 300)
            })

            # Convert the list to a tagList - this is necessary for the list of items
            # to display properly.
            do.call(tagList, plot_output_list)
        })

        #### FEATURE PLOTS DARK THEME####
        output$feature_plot_dark <- renderUI({

            feat_length <- length(features_selec())

            plot_output_list <- lapply(1:feat_length, function(i) {
                plotname <- paste("plot2", i, sep="")
                plotOutput(plotname, height = 300)
            })

            # Convert the list to a tagList - this is necessary for the list of items
            # to display properly.
            do.call(tagList, plot_output_list)
        })

        output$umap2 <- renderUI({

            feat_length <- length(features_selec())

            plot_output_list <- lapply(1:feat_length, function(i) {
                plotname <- paste("plot3", i, sep="")
                plotOutput(plotname, height = 300)
            })

            # Convert the list to a tagList - this is necessary for the list of items
            # to display properly.
            do.call(tagList, plot_output_list)
        })

        output$run_vln_plot <- renderUI({

            feat_length <- length(features_selec())

            plot_output_list <- lapply(1:feat_length, function(i) {
                plotname <- paste("plot4", i, sep="")
                plotOutput(plotname, height = 300)
            })

            # Convert the list to a tagList - this is necessary for the list of items
            # to display properly.
            do.call(tagList, plot_output_list)
        })

        output$run_dot_plot <- renderUI({

            feat_length <- length(features_selec())

            plot_output_list <- lapply(1:feat_length, function(i) {
                plotname <- paste("plot5", i, sep="")
                plotOutput(plotname, height = 300)
            })

            # Convert the list to a tagList - this is necessary for the list of items
            # to display properly.
            do.call(tagList, plot_output_list)
        })

        sc_data <- single_cell_data_reso_umap()
        for (i in 1:length(features_selec())) {

            features <- features_selec()
            local({

                my_i <- i

                plotname <- paste("plot1", my_i, sep="")
                output[[plotname]] <- renderPlot({

                    minimal <- min(sc_data[['RNA']]@data[features[my_i], ])
                    maximal <- max(sc_data[['RNA']]@data[features[my_i], ])

                    suppressMessages( FeaturePlot(sc_data,
                                                  cols = c("lightgrey", "red"),
                                                  #assay = "RNA",
                                                  features = features[my_i],
                                                  slot = input$slot_selection_feature_plot,
                                                  reduction = "umap") +
                                          scale_colour_gradient2(limits=c(minimal, maximal),
                                                                 midpoint = maximal / 2,
                                                                 low = "gray80",
                                                                 mid = "gold",
                                                                 high = "red") ) # +
                    #plot_annotation(title = paste("Gene name:", features[my_i]),
                    #               theme = theme(plot.title = element_text(size = 16,
                    #                                                      face = "bold",
                    #                                                     hjust = 0.5)) )

                })

                plotname <- paste("plot2", my_i, sep="")
                output[[plotname]] <- renderPlot({

                    minimal <- min(sc_data[['RNA']]@data[features[my_i], ])
                    maximal <- max(sc_data[['RNA']]@data[features[my_i], ])

                    suppressMessages( FeaturePlot(sc_data,
                                                  cols = c("lightgrey", "red"),
                                                  #assay = "RNA",
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

                    DimPlot(single_cell_data_reso_umap(), reduction = "umap", label = T, pt.size = .1)

                })

                plotname <- paste("plot4", my_i, sep="")
                output[[plotname]] <- renderPlot({

                    VlnPlot(sc_data,
                            features = features[my_i],
                            slot = input$slot_selection_feature_plot)

                })

                plotname <- paste("plot5", my_i, sep="")
                output[[plotname]] <- renderPlot({

                    DotPlot(sc_data,
                            features = features[my_i],
                            cols = c("lightgrey", "red")) #+ RotatedAxis()

                })

            })

        }

        showNotification("Generating additional plots",
                         duration = 25,
                         id = "m8")

        output$select_genes_add_plot_to_down_ui <- renderUI({

            filt_features <- req(input$selected_genes_for_feature_plot)

            div(class = "option-group",
                pickerInput(
                    inputId = "select_genes_add_plot_to_down",
                    label = "Select the genes that you want to download",
                    choices = sort(filt_features),
                    multiple = TRUE,
                    options = list(`actions-box` = TRUE)
                ))

        })

    })

    observeEvent(input$start_down_add_plots_tab1, {

        removeNotification(id = "m8")

        # showNotification("Downloading additional plots",
        #                  duration = NULL,
        #                  id = "mm1")

        withProgress(message = "Please wait, preparing the data for download.",
                     value = 0.5, {


                         if (file.exists("./images") == F ) {
                             dir.create('./images')
                         }

                         # Creates new folder to keep the plots organized. Adding date and time helps to avoid overwriting the data
                         path_new <- paste0("./images/one_sample_plots", format(Sys.time(),'_%Y%-%m%-%d%__%H%M%S'))

                         dir.create(path_new)
                         dir.create( paste0(path_new,"/feature_plots") )
                         dir.create( paste0(path_new,"/violin_plots") )
                         dir.create( paste0(path_new, "/dot_plots") )

                         sc_data <- single_cell_data_reso_umap()

                         genes <- req(input$select_genes_add_plot_to_down)

                         for( i in 1:length(genes) ){

                             # Saves the feature plots

                             minimal <- min(sc_data[['RNA']]@data[genes[i], ])
                             maximal <- max(sc_data[['RNA']]@data[genes[i], ])

                             p <- suppressMessages(FeaturePlot(sc_data,
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

                             ggsave(file,
                                    p,
                                    height=input$add_p_tab1_feat_height,
                                    width=input$add_p_tab1_feat_width,
                                    units="cm",
                                    dpi=as.numeric(input$add_p_tab1_feat_res))

                             # Saves the violin plots

                             file <- paste0(path_new, "/violin_plots/", genes[i], ".", input$add_p_tab1_violin_format)

                             p <- VlnPlot(sc_data,
                                          features = genes[i],
                                          slot = input$slot_selection_feature_plot)

                             ggsave(file,
                                    p,
                                    height=input$add_p_tab1_violin_height,
                                    width=input$add_p_tab1_violin_width,
                                    units="cm",
                                    dpi=as.numeric(input$add_p_tab1_violin_res))

                             # Saves the dot plots

                             file <- paste0(path_new, "/dot_plots/", genes[i], ".", input$add_p_tab1_dot_format)

                             p <- DotPlot(sc_data,
                                          features = genes[i],
                                          cols = c("lightgrey", "red")) #+ RotatedAxis()

                             ggsave(file,
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

    ## Reads what folder we have to offer the option of samples to load

    output$load_integrated_ui <- renderUI({

        rds_list <- list.files('./RDS_files/', pattern = "*.rds")

        div(class = "option-group",
            pickerInput(
                inputId = "load_integrated",
                label = "Select the file containing the integrated data",
                choices = sort(rds_list),
                multiple = FALSE,
                options = list(`actions-box` = TRUE)
            ))
        # radioButtons("load_integrated",
        #              "Select the file containing the integrated data",
        #              choices = rds_list))
    })

    samples_list_integration <- reactive({

        inFile <- input$samples_list_integration

        if (is.null(inFile))
            return(NULL)

        samples_list <- read.csv(inFile$datapath, header = T)

        samples_list

    })

    output$select_sample_tab2 <- renderUI({

        dir_list <- samples_list_integration()
        dir_list <- unique(dir_list[, 1])

        div(class = "option-group",
            pickerInput(
                inputId = "sample_folder_tab2",
                label = "Select the samples to use",
                choices = sort(as.character(dir_list)),
                multiple = TRUE,
                options = list(`actions-box` = TRUE)
            ))

    })

    single_cell_data_reac_tab2 <- eventReactive(input$load_rds_file, {

        if ( input$integration_options == 1) {

            ## Transform this option in a load file box.

            showNotification("Loading the integrated data",
                             id = "m9",
                             duration = NULL)

            sc_data <- readRDS(paste0("./RDS_files/", req(input$load_integrated)) )

            DefaultAssay(sc_data) <- "integrated"

            removeNotification(id = "m9")

            sc_data

        } else if (input$integration_options == 0 ) {

            showNotification("Loading data",
                             id = "m10",
                             duration = NULL)

            config_csv <- samples_list_integration()
            #config_csv <- config_csv[ config_csv[, 1] %in% input$sample_folder_tab2, ]

            sing_cell_list <- list()
            for( i in 1:nrow(config_csv) ) {

                data_10x_raw <- Read10X( data.dir = paste0("./data/", config_csv[ i, 1 ] ) )

                data_10x <- CreateSeuratObject(counts = data_10x_raw,
                                               project = input$int_project_name,
                                               min.cells = as.numeric(config_csv[ i, 3 ]),
                                               min.features = as.numeric(config_csv[ i, 4 ]))

                data_10x <- AddMetaData(data_10x,
                                        as.character(config_csv[ i, 2 ]), # name of sample
                                        col.name = 'treat')

                data_10x <- PercentageFeatureSet(data_10x,
                                                 pattern = input$int_regex_mito,
                                                 col.name = "percent.mt")

                data_10x <- subset(data_10x,
                                   subset = percent.mt < as.numeric(config_csv[ i, 6 ]) &
                                       nFeature_RNA < as.numeric(config_csv[ i, 5 ]))

                data_10x <- NormalizeData(data_10x, verbose = T)

                data_10x <- FindVariableFeatures(data_10x,
                                                 selection.method = input$most_var_method_integration,
                                                 nfeatures = input$n_of_var_genes_integration,
                                                 verbose = T)

                sing_cell_list[[i]] <- data_10x

            }

            removeNotification(id = "m10")

            showNotification("Integrating the data. Please wait, it can take a few minutes.",
                             id = "m12",
                             duration = NULL)

            sc_data.anchors <- FindIntegrationAnchors( object.list = sing_cell_list,
                                                       dims = 1:input$n_of_PCs_integration )
            sc_data.anchors <- IntegrateData( anchorset = sc_data.anchors, dims = 1:input$n_of_PCs_integration )

            removeNotification(id = "m12")

            DefaultAssay(sc_data.anchors) <- "integrated"

            sc_data.anchors

        }

    })

    output$download_int_data <- downloadHandler(

        filename = function() {
            paste(input$dataset, ".rds", sep = "")
        },
        content = function(file) {

            withProgress(message = "Please wait, preparing the data for download.",
                         value = 0.5, {

                             saveRDS(single_cell_data_reac_tab2(), file)

                         })

            #removeNotification(id = "D1")
        }

    )

    # Visualize QC metrics as a violin plot
    output$VlnPlot_tab2 <- renderPlot({

        data_set <- single_cell_data_reac_tab2()

        VlnPlot(data_set,
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                ncol = 3,
                assay = "RNA",
                split.plot = F)

    })

    # Filtering and removing cells based on counts and % of mitochondrial
    single_cell_data_filt_tab2 <- reactive({

        data_sc <- single_cell_data_reac_tab2()

        DefaultAssay(data_sc) <- "RNA"

        if ( !is.na(input$min_count_tab2) ) {

            # Filtering features and cells based on the counts and % of mito contamination.
            data_sc <- subset(data_sc,
                              subset = nFeature_RNA > input$min_count_tab2)

        }

        if ( !is.na(input$max_count_tab2) ) {

            # Filtering features and cells based on the counts and % of mito contamination.
            data_sc <- subset(data_sc,
                              subset =  nFeature_RNA < input$max_count_tab2)

        }

        if ( !is.na(input$max_mito_perc_tab2) ) {

            # Filtering features and cells based on the counts and % of mito contamination.
            data_sc <- subset(data_sc,
                              subset =  percent.mt < input$max_mito_perc_tab2)

        }

        DefaultAssay(data_sc) <- "integrated"

        data_sc

    })

    observeEvent(input$run_vinplot_tab2, {

        output$VlnPlot_filt_tab2 <- renderPlot({

            VlnPlot(single_cell_data_filt_tab2(),
                    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    ncol = 3,
                    assay = "RNA",
                    split.plot = F)

        })

        output$p5_down <- downloadHandler(

            filename = function() {
                paste(input$dataset, ".", input$p5_format, sep = "")
            },
            content = function(file) {

                withProgress(message = "Please wait, preparing the data for download.",
                             value = 0.5, {

                                 height <- as.numeric(input$p5_height)
                                 width <- as.numeric(input$p5_width)
                                 res <- as.numeric(input$p5_res)


                                 p <- VlnPlot(single_cell_data_filt_tab2(),
                                              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                                              ncol = 3,
                                              assay = "RNA",
                                              split.plot = F)

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

    # Integrate assay
    single_cell_data_scaled_tab2 <- eventReactive(input$run_pca_tab2, {

        single_cell_data_filt_tab2 <- single_cell_data_filt_tab2()

        # DefaultAssay(single_cell_data_filt_tab2) <- "RNA"
        DefaultAssay(single_cell_data_filt_tab2) <- "integrated"

        showNotification("Scalling the data",
                         duration = NULL,
                         id = "tab2_m4")

        data_sc <- FindVariableFeatures(single_cell_data_filt_tab2,
                                        selection.method = input$most_var_method_tab2,
                                        nfeatures = input$n_of_var_genes_tab2)

        allgenes <- rownames(data_sc)
        data_sc <- ScaleData(data_sc, features = allgenes)

        removeNotification(id = "tab2_m4")
        data_sc
    })

    single_cell_data_pca_tab2 <- eventReactive(c(input$run_pca_tab2, input$rerun_after_filtering_tab2), {

        if (input$filter_clusters_tab2 == 1) {

            data <- single_cell_data_scaled_tab2_filtered()

        } else if ( input$filter_clusters_tab2 == 0 ) {

            data <- single_cell_data_scaled_tab2()

        } else {

            print("Something is wrong!")
        }

        showNotification("Running PCA",
                         duration = NULL,
                         id = "m5")
        #data_sc <- single_cell_data_scaled()

        ret_data <-  RunPCA(data,
                            features = VariableFeatures(object = data), verbose = F)

        removeNotification(id = "m5")
        ret_data
    })

    # single_cell_data_pca2_tab2 <- eventReactive(c(input$run_pca_tab2, input$rerun_after_filtering_tab2), {
    #
    #     if (input$normaliz_method_tab2 == "SCTransform") {
    #
    #         single_cell_data_pca2 <- single_cell_data_pca_SCT_tab2()
    #
    #     } else if (input$normaliz_method_tab2 == "LogNormalize") {
    #
    #         single_cell_data_pca2 <- single_cell_data_pca_tab2()
    #
    #     }
    #
    # })

    # Elbow plot (also triggers the PCA)
    observeEvent(c(input$run_pca_tab2, input$rerun_after_filtering_tab2), {

        # Generates the elbow plot showing the PCs #
        output$n_of_PCAs_tab2 <- renderPlot({

            data_sc <- single_cell_data_pca_tab2()

            ElbowPlot(data_sc, ndims = 50, reduction = "pca")

        })

        output$p6_down <- downloadHandler(

            filename = function() {
                paste(input$dataset, ".", input$p6_format, sep = "")
            },
            content = function(file) {

                withProgress(message = "Please wait, preparing the data for download.",
                             value = 0.5, {

                                 height <- as.numeric(input$p6_height)
                                 width <- as.numeric(input$p6_width)
                                 res <- as.numeric(input$p6_res)

                                 data_sc <- single_cell_data_pca_tab2()

                                 p <- ElbowPlot(data_sc, ndims = 50, reduction = "pca")

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

    ## clustering tab2

    single_cell_data_neigh_tab2 <- eventReactive(input$run_clustering_tab2, {
        showNotification("Running the clustering step",
                         duration = NULL,
                         id = "m6")

        data_sc <- single_cell_data_pca_tab2()

        FindNeighbors(data_sc, dims = 1:input$n_of_PCs_tab2)

    })

    single_cell_data_reso_tab2 <- eventReactive(input$run_clustering_tab2, {

        sc_data <- single_cell_data_neigh_tab2()

        sc_data <-  FindClusters(sc_data, resolution = input$resolution_clust_tab2)

        # After finishing all steps that depends on the integrated data, we can swich for the RNA assay
        DefaultAssay(sc_data) <- "RNA"

        allgenes <- rownames(sc_data)
        sc_data <- ScaleData(sc_data, features = allgenes)

    })

    single_cell_data_reso_umap_tab2 <- eventReactive(input$run_clustering_tab2, {

        RunUMAP(single_cell_data_reso_tab2(), dims = 1:input$n_of_PCs_tab2)

    })

    observeEvent(input$run_clustering_tab2, {

        # output$tSNE_tab2 <- renderPlot({
        #
        #     #data_sc <- single_cell_data_pca()
        #
        #     tSNE_plot <- RunTSNE(single_cell_data_reso_tab2(), dims = 1:input$n_of_PCs_tab2)
        #     DimPlot(tSNE_plot, reduction = "tsne", label = T, pt.size = .1)
        #
        # })

        output$umap_tab2 <- renderPlot({

            #data_sc <- single_cell_data_pca()

            #umap_plot <- RunUMAP(single_cell_data_reso(), dims = 1:input$n_of_PCs)
            DimPlot(single_cell_data_reso_umap_tab2(), reduction = "umap", label = T, pt.size = .1)

        })

        output$umap_three_samples_comb <- renderPlot({

            DimPlot(single_cell_data_reso_umap_tab2(),
                    reduction = "umap",
                    label = F,
                    group.by = "treat",
                    pt.size = .1
            )
        })

        output$umap_three_samples <- renderPlot({

            DimPlot(single_cell_data_reso_umap_tab2(),
                    reduction = "umap",
                    label = T,
                    split.by = "treat", pt.size = .1
            )
        })

        output$p7_down <- downloadHandler(

            filename = function() {
                paste(input$dataset, ".", input$p7_format, sep = "")
            },
            content = function(file) {

                withProgress(message = "Please wait, preparing the data for download.",
                             value = 0.5, {

                                 height <- as.numeric(input$p7_height)
                                 width <- as.numeric(input$p7_width)
                                 res <- as.numeric(input$p7_res)

                                 if ( input$p7_down_opt == "UMAP" ) {

                                     p <- DimPlot(single_cell_data_reso_umap_tab2(),
                                                  reduction = "umap",
                                                  label = T,
                                                  pt.size = .1)

                                 } else if ( input$p7_down_opt == "UMAP1" ) {

                                     p <- DimPlot(single_cell_data_reso_umap_tab2(),
                                                  reduction = "umap",
                                                  label = F,
                                                  group.by = "treat",
                                                  pt.size = .1
                                     )

                                 } else if ( input$p7_down_opt == "UMAP2" ) {

                                     p <- DimPlot(single_cell_data_reso_umap_tab2(),
                                                  reduction = "umap",
                                                  label = T,
                                                  split.by = "treat", pt.size = .1
                                     )

                                 }

                                 ggsave(file,
                                        p,
                                        height=height,
                                        width=width,
                                        units="cm",
                                        dpi=res)

                             })
            }

        )


        output$cluster_size_tab2 <- renderPrint({

            sing_cell_data <- single_cell_data_reso_umap_tab2()

            sc_meta <- as.data.frame(sing_cell_data[[]])
            sc_meta$cellcluster <- rownames(sc_meta)
            sc_meta <- sc_meta[, c("cellcluster", "seurat_clusters")]

            # returns the number of cells in each cluster
            cell_per_cluster <- sc_meta$seurat_clusters

            removeNotification(id = "m6")

            table(sc_meta$seurat_clusters)


        })

    } )

    clusters_single_cell_data_reso_umap_tab2 <- reactive({

        sc_data <- single_cell_data_reso_umap_tab2()
        as.numeric( unique( as.character(sc_data@meta.data$seurat_clusters ) ) )


    })

    treat_single_cell_data_reso_umap_tab2 <- reactive({

        sc_data <- single_cell_data_reso_umap_tab2()
        unique( as.character(sc_data@meta.data$treat ) )

    })

    output$cluster_list_tab2_ui <- renderUI ({

        clusters <- clusters_single_cell_data_reso_umap_tab2()

        pickerInput(
            inputId = "cluster_list_tab2",
            label = "Choose clusters to select or exclude",
            choices = sort(clusters),
            multiple = TRUE,
            options = list(`actions-box` = TRUE)
        )

    })

    to_filter_tab2 <- reactive({

        if ( input$filter_clusters_opt_tab2 == "select" ) {

            to_filter <- subset(single_cell_data_reso_umap_tab2(),
                                idents = as.numeric(req(input$cluster_list_tab2))
            )

            to_filter_ch <- to_filter@meta.data
            to_filter_ch <- rownames(to_filter_ch)

            to_filter_ch

        } else if ( input$filter_clusters_opt_tab2 == "exclude" ) {

            to_filter <- subset(single_cell_data_reso_umap_tab2(),
                                idents = as.numeric(req(input$cluster_list_tab2)),
                                invert = TRUE)

            to_filter_ch <- to_filter@meta.data
            to_filter_ch <- rownames(to_filter_ch)

            to_filter_ch
        }


    })

    single_cell_data_scaled_tab2_filtered <- eventReactive(input$rerun_after_filtering_tab2, {

        showNotification("Scalling the data",
                         duration = NULL,
                         id = "tab2_m4")

        sc_data <- single_cell_data_filt_tab2()

        DefaultAssay(sc_data) <- "integrated"

        if ( input$filter_clusters_tab2 == 1 ) {

            sc_data <- subset(sc_data,
                              cells = to_filter_tab2())

            # data_sc <- FindVariableFeatures(sc_data,
            #                                 selection.method = input$most_var_method_tab2,
            #                                 nfeatures = input$n_of_var_genes_tab2)
            #
            allgenes <- rownames(sc_data)
            sc_data <- ScaleData(sc_data, features = allgenes)

            sc_data
        }

        # } else if (input$filter_clusters_tab2 == 0 ) {
        #
        #
        #     # data_sc <- FindVariableFeatures(sc_data,
        #     #                                 selection.method = input$most_var_method_tab2,
        #     #                                 nfeatures = input$n_of_var_genes_tab2)
        #     #
        #     # allgenes <- rownames(data_sc)
        #     # sc_data <- ScaleData(data_sc, features = allgenes)
        #
        # }

        removeNotification(id = "tab2_m4")
        sc_data
    })


    # output$range <- renderText({to_filter_tab2()})

    output$downloadRDS_tab2 <- downloadHandler(

        filename = function() {
            paste(input$dataset, ".rds", sep = "")
        },
        content = function(file) {

            withProgress(message = "Please wait, preparing the data for download.",
                         value = 0.5, {

                             saveRDS(single_cell_data_reso_umap_tab2(), file)

                         })

            #removeNotification(id = "D1")
        }

    )

    output$find_markers_clust_id_tab2_ui <- renderUI ({

        clusters <- clusters_single_cell_data_reso_umap_tab2()

        pickerInput(
            inputId = "find_markers_clust_id_tab2",
            label = "Select the cluster of interest",
            choices = sort(clusters),
            multiple = FALSE,
            options = list(`actions-box` = TRUE)
        )

    })

    output$find_markers_clust_ID1_tab2_ui <- renderUI ({

        clusters <- clusters_single_cell_data_reso_umap_tab2()

        pickerInput(
            inputId = "find_markers_clust_ID1_tab2",
            label = "Select the cluster of interest",
            choices = sort(clusters),
            multiple = FALSE,
            options = list(`actions-box` = TRUE)
        )

    })

    output$find_markers_clust_ID2_tab2_ui <- renderUI ({

        clusters <- clusters_single_cell_data_reso_umap_tab2()

        # Exclude the cluster already selected in the option 1, since the two cluster must be different

        clusters <- clusters[ !clusters == req(input$find_markers_clust_ID1_tab2)]

        pickerInput(
            inputId = "find_markers_clust_ID2_tab2",
            label = "Select the cluster(s) to compare",
            choices = sort(clusters),
            multiple = T,
            options = list(`actions-box` = TRUE)
        )

    })

    output$find_markers_or_DE_tab2_cluster_ui <- renderUI ({

        clusters <- clusters_single_cell_data_reso_umap_tab2()

        pickerInput(
            inputId = "find_markers_or_DE_tab2_cluster",
            label = "Select the cluster of interest",
            choices = sort(clusters),
            multiple = FALSE,
            options = list(`actions-box` = TRUE)
        )

    })

    output$find_markers_or_DE_tab2_treat1_ui <- renderUI({

        treat <- treat_single_cell_data_reso_umap_tab2()

        pickerInput(
            inputId = "find_markers_or_DE_tab2_treat1",
            label = "Select the sample/treatment of interest",
            choices = treat,
            multiple = FALSE,
            options = list(`actions-box` = TRUE)
        )

    })

    output$find_markers_or_DE_tab2_treat2_ui  <- renderUI({

        treat <- treat_single_cell_data_reso_umap_tab2()

        treat <- treat[!treat %in% req(input$find_markers_or_DE_tab2_treat1)]

        pickerInput(
            inputId = "find_markers_or_DE_tab2_treat2",
            label = "Select the sample/treatment of interest",
            choices = treat,
            multiple = T,
            options = list(`actions-box` = TRUE)
        )

    })


    # Identification of markers (tab 2)
    markers_tab2 <- eventReactive( input$run_ident_markers_tab2, {

        showNotification("Identifing markers or D.E. genes",
                         duration = NULL,
                         id = "tab2_n1")

        sc_data <- single_cell_data_reso_umap_tab2()

        if (input$find_markers_or_DE_tab2 == 0 ) {

            if ( input$find_markers_tab2_opt == 0 ) {

                clusters <- as.numeric(as.character(unique(sc_data@meta.data$seurat_clusters)))
                markers_tab2 <- data.frame()

                for (i in clusters) {

                    markers <- FindConservedMarkers(sc_data,
                                                    ident.1 = i,
                                                    grouping.var = "treat",
                                                    verbose = T,
                                                    assay = "RNA")

                    markers_tab2 <- rbind(markers_tab2, markers)

                }

                markers_tab2

            } else if ( input$find_markers_tab2_opt == 1 ) {

                markers_tab2 <- FindConservedMarkers(sc_data,
                                                     ident.1 = req(input$find_markers_clust_id_tab2),
                                                     grouping.var = "treat",
                                                     assay = "RNA")

            } else if ( input$find_markers_tab2_opt == 2 ) {

                markers_tab2 <- FindConservedMarkers(sc_data,
                                                     ident.1 = req(input$find_markers_clust_ID1_tab2),
                                                     ident.2 = req(input$find_markers_clust_ID2_tab2),
                                                     grouping.var = "treat",
                                                     assay = "RNA")


            }

        } else if (input$find_markers_or_DE_tab2 == 1 ) {

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
                                        verbose = FALSE,
            )

            markers_tab2 <- markers_tab2[markers_tab2$p_val_adj < input$find_markers_or_DE_tab2_pvalue, ]

        }

        markers_tab2$geneID <- rownames(markers_tab2)
        markers_tab2

    })

    output$markers_tab2_react <- renderReactable({

        markers_tab2 <- markers_tab2()
        markers_tab2 <- markers_tab2[ , c( ncol(markers_tab2), 1: (ncol(markers_tab2) -1) ) ]

        removeNotification(id = "tab2_n1")

        reactable(markers_tab2,
                  defaultColDef = colDef(align = "center"),
                  filterable = TRUE,
                  bordered = TRUE,
                  highlight = TRUE,
                  searchable = TRUE,
                  showPageSizeOptions = TRUE,
                  defaultPageSize = 10,
                  pageSizeOptions = c(5, 10, 25, 50, 100, 200),
                  showPagination = TRUE,
                  rownames = FALSE)

    })

    output$download_markers_tab2 <- downloadHandler(

        filename = function() {
            paste(input$dataset, ".csv", sep = "")
        },
        content = function(file) {

            withProgress(message = "Please wait, preparing the data for download.",
                         value = 0.5, {

                             markers_tab2 <- markers_tab2()
                             markers_tab2 <- markers_tab2[ , c( ncol(markers_tab2), 1: (ncol(markers_tab2) -1) ) ]

                             write.csv(markers_tab2, file, row.names = FALSE)

                         })
        }

    )

    single_cell_data_heatmap_tab2 <- eventReactive(input$run_heatmap_tab2, {

        showNotification("Generating Heatmap",
                         duration = NULL,
                         id = "m7")

        if (input$normaliz_method_tab2 == "SCTransform") {

            # data_sc <- single_cell_data_reso_umap_tab2()
            #
            # data_sc <- NormalizeData(data_sc,
            #                          assay = "RNA",
            #                          normalization.method = "LogNormalize",
            #                          scale.factor = input$scale_factor_tab2)
            #
            # data_sc <- FindVariableFeatures(data_sc,
            #                                 assay = "RNA",
            #                                 selection.method = input$most_var_method,
            #                                 nfeatures = input$n_of_var_genes)
            #
            # allgenes <- rownames(data_sc)
            #
            # ScaleData(data_sc,
            #           assay = "RNA",
            #           features = allgenes)

        } else if (input$normaliz_method_tab2 == "LogNormalize") {

            sc_data <- single_cell_data_reso_umap_tab2()

            DefaultAssay(sc_data) <- "RNA"

            #    allgenes <- rownames(data_sc)

            #    ScaleData(data_sc,
            #              assay = "RNA",
            #              features = allgenes)

            sc_data

        }

    })

    ## Gene expression
    features_tab2 <- reactive({

        inFile <- input$markers_list_tab2

        if (is.null(inFile))
            return(NULL)

        features_tab2 <- read.csv(inFile$datapath, header = F)

        if ( ncol(features_tab2) == 2 ) {

            features_tab2 <- features_tab2 %>%
                dplyr::rename(GeneID = V1,
                              Group = V2)

        } else if ( ncol(features_tab2) > 2 ) {

            features_tab2 <- features_tab2 %>%
                dplyr::rename(GeneID = V1,
                              Group = V2,
                              Name = V3)

        }

        features_tab2

    })

    # Load the file and offers the parameters for heatmap
    observeEvent(input$load_markers_tab2, {

        # Painel that will apper after loading the list of markers containing the filter options
        output$marker_group_selec_tab2 = renderUI({

            groups_names <- features_tab2()
            groups_names <- unique(groups_names$Group)

            div(class = "option-group",
                pickerInput(
                    inputId = "features_group_tab2",
                    label = "Select the group of markers to test",
                    choices = sort( as.character(groups_names) ),
                    multiple = TRUE,
                    options = list(`actions-box` = TRUE)
                ))

        })

        output$marker_filter_genes_q_tab2 = renderUI({

            div(class = "option-group",
                radioButtons("filter_genes_q_tab2",
                             "Visualization options",
                             choices = list("Show all genes" = 1,
                                            "Select genes to show" = 0),
                             selected = 1))
        })

        # Define if using names of v5 or common names
        output$marker_genes_ids_tab2 = renderUI({

            div(class = "option-group",
                radioButtons("genes_ids_tab2",
                             "ID options",
                             choices = list("Use IDs" = "v5",
                                            "Use Name" = "name"),
                             selected = "v5"))
        })

        output$marker_genes_selec_tab2 = renderUI({

            groups_names <- features_tab2()

            genes_names <- dplyr::filter(groups_names,
                                         Group %in% req(input$features_group_tab2))

            if (input$genes_ids_tab2 == "v5") {

                genes_names <- unique(genes_names$GeneID)

            } else if (input$genes_ids_tab2 == "name") {

                genes_names <- unique(genes_names$Name)

            }

            genes_names <- genes_names[!is.na(genes_names)]

            div(class = "option-group",
                pickerInput(
                    inputId = "selected_genes_tab2",
                    label = "Select the genes to show",
                    choices = sort(as.character(genes_names)),
                    multiple = TRUE,
                    options = list(`actions-box` = TRUE)
                ))

        })

    })

    # Filter accordly with the parameters selected by the user
    filt_features_tab2 <- eventReactive(input$run_heatmap_tab2, {

        features_f <- features_tab2()

        features_f <- dplyr::filter(features_f,
                                    Group %in% req(input$features_group_tab2))

        if (input$filter_genes_q_tab2 == 0) {

            if (input$genes_ids_tab2 == "v5") {

                features_f <- features_f[features_f$GeneID %in% req(input$selected_genes_tab2), ]

            } else if (input$genes_ids == "name") {

                features_f <- features_f[features_f$Name %in% req(input$selected_genes_tab2), ]

            }

        }

        features_f

    })

    sc_data_av_react_tab2 <- eventReactive(input$run_heatmap_tab2, {

        features <- filt_features_tab2()
        sc_data <- single_cell_data_heatmap_tab2()

        # get the expression values
        sc_data_av <- AverageExpression(sc_data,
                                        assays = "RNA",
                                        # features = features,
                                        slot = input$slot_selection_heatmap_tab2)


        sc_data_av <- as.matrix(sc_data_av[[1]])

        removeNotification(id = "m7")
        sc_data_av

    })

    observeEvent(input$run_heatmap_tab2, {

        # This calculate the height of the plot based on the n of genes
        heatmap_n_genes_tab2 <- reactive({

            features <- filt_features_tab2()
            # get the expression values
            sc_data_av <- sc_data_av_react_tab2()
            features_selec <- unique(features$GeneID)

            if (nrow(features) == 1) {

                sc_data_av_feat <- sc_data_av[rownames(sc_data_av) %in% features_selec, ]
                sc_data_av_feat <- as.data.frame(t(sc_data_av_feat))

                size_hetmap <- as.numeric(nrow(sc_data_av_feat))

            } else {

                sc_data_av_feat <- sc_data_av[rownames(sc_data_av) %in% features_selec, ]

                size_hetmap <- as.numeric(nrow(sc_data_av_feat))

            }

            size_hetmap

        })

        heatmap_Height_tab2 <- reactive( 150 + ( 20 * heatmap_n_genes_tab2() ) )

        heat_map_prep_tab2 <- reactive({

            features <- filt_features_tab2()
            # get the expression values
            sc_data_av <- sc_data_av_react_tab2()
            features_selec <- unique(features$GeneID)

            if (nrow(features) == 1) {

                sc_data_av_feat <- sc_data_av[rownames(sc_data_av) %in% features_selec, ]
                sc_data_av_feat <- as.data.frame(t(sc_data_av_feat))
                rownames(sc_data_av_feat) <- features_selec[1, 1]

            } else {

                sc_data_av_feat <- sc_data_av[rownames(sc_data_av) %in% features_selec, ]

            }

            sc_data_av_feat

        })

        # Generates the heatmap plot
        output$heat_map_tab2 <- renderPlot({

            heat_map_prep_tab2 <- heat_map_prep_tab2()

            ### Drawing heat map ###
            Heatmap(heat_map_prep_tab2, border = TRUE,
                    rect_gp = gpar(col = "white", lwd = 2),
                    column_title = "Clusters",
                    column_title_side = "bottom",
                    name = "Expression",
                    show_row_dend = T)

        }, height = heatmap_Height_tab2() )

        output$heat_map_ui_tab2 <- renderUI({

            plotOutput("heat_map_tab2", height = heatmap_Height_tab2())

        })

        output$marker_to_feature_plot_tab2 = renderUI({

            sc_data_av_react <- sc_data_av_react_tab2()
            genes_names <- filt_features_tab2()

            #    genes_names <- genes_names[paste0("gene:", genes_names$GeneID) %in% rownames(sc_data_av_react), ]

            genes_names <- genes_names[genes_names$GeneID %in% rownames(sc_data_av_react), ]

            if (input$genes_ids_tab2 == "v5") {

                genes_names <- unique(genes_names$GeneID)

            } else if (input$genes_ids_tab2 == "name") {

                genes_names <- unique(genes_names$Name)

            }

            div(class = "option-group",
                pickerInput(
                    inputId = "selected_genes_for_feature_plot_tab2",
                    label = "Select the genes to feature plots",
                    choices = sort(as.character(genes_names)),
                    multiple = TRUE,
                    options = list(`actions-box` = TRUE)
                ))

        })

        output$p8_down <- downloadHandler(

            filename = function() {
                paste(input$dataset, ".", input$p8_format, sep = "")
            },
            content = function(file) {

                withProgress(message = "Please wait, preparing the data for download.",
                             value = 0.5, {

                                 height <- as.numeric(input$p8_height)
                                 width <- as.numeric(input$p8_width)
                                 res <- as.numeric(input$p8_res)

                                 ### Drawing heat map ###

                                 heat_map_prep_tab2 <- heat_map_prep_tab2()

                                 ### Drawing heat map ###
                                 p <- Heatmap(heat_map_prep_tab2, border = TRUE,
                                              rect_gp = gpar(col = "white", lwd = 2),
                                              column_title = "Clusters",
                                              column_title_side = "bottom",
                                              name = "Expression",
                                              show_row_dend = T)

                                 # Complex heatmap does not work well with ggsave. So, using grid.grap to make it compatible
                                 gb = grid.grabExpr(draw(p))

                                 ggsave(file,
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

        features <- filt_features_tab2()

        if (input$genes_ids_tab2 == "v5") {

            features_f <- dplyr::filter(features,
                                        GeneID %in% req(input$selected_genes_for_feature_plot_tab2))

        }
        else if (input$genes_ids_tab2 == "name") {

            features_f <- dplyr::filter(features,
                                        Name %in% req(input$selected_genes_for_feature_plot_tab2))

        }

        #    paste0("gene:", unique(features_f$GeneID))
        as.character(unique(features_f$GeneID))

    })

    observeEvent(input$run_feature_plot_tab2, {

        ### FEATURE PLOTS ####

        output$feature_plot_tab2 <- renderUI({

            feat_length <- length(features_selec_tab2())

            plot_output_list <- lapply(1:feat_length, function(i) {
                plotname <- paste("plot1_tab2", i, sep="")
                plotOutput(plotname, height = 300)
            })

            # Convert the list to a tagList - this is necessary for the list of items
            # to display properly.
            do.call(tagList, plot_output_list)
        })

        # #### FEATURE PLOTS DARK THEME####
        # output$feature_plot_dark_tab2 <- renderUI({
        #
        #     feat_length <- length(features_selec_tab2())
        #
        #     plot_output_list <- lapply(1:feat_length, function(i) {
        #         plotname <- paste("plot2_tab2", i, sep="")
        #         plotOutput(plotname, height = 300)
        #     })
        #
        #     # Convert the list to a tagList - this is necessary for the list of items
        #     # to display properly.
        #     do.call(tagList, plot_output_list)
        # })

        output$umap2_tab2 <- renderUI({

            feat_length <- length(features_selec_tab2())

            plot_output_list <- lapply(1:feat_length, function(i) {
                plotname <- paste("plot3_tab2", i, sep="")
                plotOutput(plotname, height = 300)
            })

            # Convert the list to a tagList - this is necessary for the list of items
            # to display properly.
            do.call(tagList, plot_output_list)
        })

        output$run_vln_plot_tab2 <- renderUI({

            feat_length <- length(features_selec_tab2())

            plot_output_list <- lapply(1:feat_length, function(i) {
                plotname <- paste("plot4_tab2", i, sep="")
                plotOutput(plotname, height = 300)
            })

            # Convert the list to a tagList - this is necessary for the list of items
            # to display properly.
            do.call(tagList, plot_output_list)
        })

        output$run_dot_plot_tab2 <- renderUI({

            feat_length <- length(features_selec_tab2())

            plot_output_list <- lapply(1:feat_length, function(i) {
                plotname <- paste("plot5_tab2", i, sep="")
                plotOutput(plotname, height = 300)
            })

            # Convert the list to a tagList - this is necessary for the list of items
            # to display properly.
            do.call(tagList, plot_output_list)
        })

        sc_data_tab2 <- single_cell_data_reso_umap_tab2()
        DefaultAssay(sc_data_tab2) <- "RNA"

        for (i in 1:length(features_selec_tab2())) {

            features <- features_selec_tab2()
            local({

                my_i <- i

                plotname <- paste("plot1_tab2", my_i, sep="")
                output[[plotname]] <- renderPlot({

                    # FeaturePlot(sc_data_tab2,
                    #             cols = c("lightgrey", "red"),
                    #             split.by = "treat",
                    #             #assay = "RNA",
                    #             features = features[my_i],
                    #             #slot = input$slot_selection_feature_plot_tab2,
                    #             reduction = "umap",
                    #             min.cutoff = "q10",
                    #             max.cutoff = "q90",
                    #             pt.size = .1)
                    #

                    p_list <- FeaturePlotSingle(sc_data_tab2,
                                                feature = features[my_i],
                                                metadata_column = "treat",
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

                plotname <- paste("plot2_tab2", my_i, sep="")
                output[[plotname]] <- renderPlot({

                    suppressMessages( FeaturePlot(sc_data_tab2,
                                                  cols = c("lightgrey", "red"),
                                                  features = features[my_i],
                                                  #assay = "RNA",
                                                  #slot = input$slot_selection_feature_plot_tab2,
                                                  reduction = "umap",
                                                  min.cutoff = "q10",
                                                  max.cutoff = "q90",
                                                  pt.size = .1) +
                                          DarkTheme() )

                })

                plotname <- paste("plot3_tab2", my_i, sep="")
                output[[plotname]] <- renderPlot({

                    DimPlot(single_cell_data_reso_umap_tab2(), reduction = "umap", label = T, pt.size = .1)

                })

                plotname <- paste("plot4_tab2", my_i, sep="")
                output[[plotname]] <- renderPlot({

                    VlnPlot(sc_data_tab2,
                            features = features[my_i],
                            split.by = "treat",
                            assay = "RNA"
                            #slot = input$slot_selection_feature_plot_tab2
                    )

                })

                plotname <- paste("plot5_tab2", my_i, sep="")
                output[[plotname]] <- renderPlot({

                    DotPlot(sc_data_tab2,
                            features = features[my_i],
                            cols = c("lightgrey", "red"),
                            split.by = "treat",
                            assay = "RNA") #+ RotatedAxis()

                })

            })

        }
        showNotification("Generating additional plots",
                         duration = 30,
                         id = "m8")

        output$select_genes_add_plot_to_down_tab2_ui <- renderUI({

            filt_features <- req(input$selected_genes_for_feature_plot_tab2)

            div(class = "option-group",
                pickerInput(
                    inputId = "select_genes_add_plot_to_down_tab2",
                    label = "Select the genes that you want to download",
                    choices = sort(filt_features),
                    multiple = TRUE,
                    options = list(`actions-box` = TRUE)
                ))

        })
    })


    observeEvent(input$start_down_add_plots_tab2, {

        removeNotification(id = "m8")

        # showNotification("Downloading additional plots",
        #                  duration = NULL,
        #                  id = "tab2_mm1")

        withProgress(message = "Please wait, preparing the data for download.",
                     value = 0.5, {


                         if (file.exists("./images") == F ) {
                             dir.create('./images')
                         }

                         # Creates new folder to keep the plots organized. Adding date and time helps to avoid overwriting the data
                         path_new <- paste0("./images/integrated_sample_plots", format(Sys.time(),'_%Y%-%m%-%d%__%H%M%S'))

                         dir.create(path_new)
                         dir.create( paste0(path_new,"/feature_plots") )
                         dir.create( paste0(path_new,"/violin_plots") )
                         dir.create( paste0(path_new, "/dot_plots") )

                         sc_data <- single_cell_data_reso_umap_tab2()

                         genes <- req(input$select_genes_add_plot_to_down_tab2)

                         for( i in 1:length(genes) ){

                             # Saves the feature plots

                             minimal <- min(sc_data[['RNA']]@data[genes[i], ])
                             maximal <- max(sc_data[['RNA']]@data[genes[i], ])

                             p <- suppressMessages(FeaturePlot(sc_data,
                                                               cols = c("lightgrey", "red"),
                                                               features = genes[i],
                                                               slot = input$slot_selection_feature_plot_tab2,
                                                               reduction = "umap") +
                                                       scale_colour_gradient2(limits=c(minimal, maximal),
                                                                              midpoint = maximal / 2,
                                                                              low = "gray80",
                                                                              mid = "gold",
                                                                              high = "red"))

                             file <- paste0(path_new, "/feature_plots/", genes[i], ".", input$add_p_tab2_feat_format)

                             ggsave(file,
                                    p,
                                    height=input$add_p_tab2_feat_height,
                                    width=input$add_p_tab2_feat_width,
                                    units="cm",
                                    dpi=as.numeric(input$add_p_tab2_feat_res))

                             # Saves the violin plots

                             file <- paste0(path_new, "/violin_plots/", genes[i], ".", input$add_p_tab2_violin_format)

                             p <- VlnPlot(sc_data,
                                          features = genes[i],
                                          slot = input$slot_selection_feature_plot_tab2)

                             ggsave(file,
                                    p,
                                    height=input$add_p_tab2_violin_height,
                                    width=input$add_p_tab2_violin_width,
                                    units="cm",
                                    dpi=as.numeric(input$add_p_tab2_violin_res))

                             # Saves the dot plots

                             file <- paste0(path_new, "/dot_plots/", genes[i], ".", input$add_p_tab2_dot_format)

                             p <- DotPlot(sc_data,
                                          features = genes[i],
                                          cols = c("lightgrey", "red")) #+ RotatedAxis()

                             ggsave(file,
                                    p,
                                    height=input$add_p_tab2_dot_height,
                                    width=input$add_p_tab2_dot_width,
                                    units="cm",
                                    dpi=as.numeric(input$add_p_tab2_dot_res))


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
            pickerInput(
                inputId = "load_integrated_tab3",
                label = "Select the file containing the clustered data",
                choices = sort(rds_list),
                multiple = FALSE,
                options = list(`actions-box` = TRUE)
            ))
    })

    sc_data_traj_inf <- eventReactive(input$run_ti_model, {

        showNotification("Loading the clustered data",
                         id = "tab3_m1",
                         duration = NULL)

        # Test if using a data with clustering or not
        # use_processed_data = "yes"
        # if ( use_processed_data == "yes" ) {
        #
        #     sc_data <- readRDS(paste0("./RDS_files/", input$load_integrated_tab3) )
        #
        # } else if ( use_processed_data == "no" ) {
        #
        #     ## Implement a way to do the analysis without using the other tabs
        #
        # }

        sc_data <- readRDS(paste0("./RDS_files/", req(input$load_integrated_tab3)) )

        removeNotification(id = "tab3_m1")

        sc_data

    } )

    expressed_genes <- reactive({

        sce <- as.SingleCellExperiment(sc_data_traj_inf())
        sce@assays@data@listData$counts@Dimnames[[1]]

    })

    # extract the expression matrix from the Seurat object
    object_expression <- reactive({

        sing_cell_data <- sc_data_traj_inf()

        Matrix::t(as(as.matrix(sing_cell_data@assays$RNA@data), 'sparseMatrix'))

    })

    object_counts <- reactive({

        sing_cell_data <- sc_data_traj_inf()

        Matrix::t(as(as.matrix(sing_cell_data@assays$RNA@counts), 'sparseMatrix'))

    })

    dataset_inf <- reactive({

        wrap_expression(
            counts = object_counts(),
            expression = object_expression()
        )

    })

    # Extracts meta data to fill the dyno object
    sc_meta <- reactive({

        sing_cell_data <- sc_data_traj_inf()

        sing_cell_data[[]] %>%
            mutate(cells = rownames(.))

    })

    ## Clustering
    sc_meta_cluster <- reactive({

        sc_meta <- sc_meta()

        data.frame(cell_id = sc_meta$cells,
                   group_id = sc_meta$seurat_clusters)

    })

    sc_cells_time_vec <- reactive ({

        sc_meta <- sc_meta()

        # if (input$ti_sample_number == 1 & input$ti_graphs_color_choice == 0 ) {
        #
        #     sc_cells_time_vec <- as.character(sc_meta$treat)
        #     names(sc_cells_time_vec) <- sc_meta$cells
        #
        # } else if (input$ti_sample_number == 1 & input$ti_graphs_color_choice == 1) {
        #
        #     sc_cells_time_vec <- as.numeric(sc_meta$treat)
        #     names(sc_cells_time_vec) <- sc_meta$cells
        #
        # }

        sc_cells_time_vec <- as.character(sc_meta$treat)
        names(sc_cells_time_vec) <- sc_meta$cells

        sc_cells_time_vec
    })

    output$ti_methods_list_ui <- renderUI ({

        ti_methods_list <- dynmethods::methods$method_id

        div(class = "option-group",
            pickerInput(
                inputId = "ti_methods_list",
                label = "Select the dynverse method to execute",
                choices = sort(as.character(ti_methods_list)),
                multiple = FALSE,
                selected = "slingshot",
                options = list(`actions-box` = TRUE)
            ))

    })

    model <- eventReactive(input$run_ti_model, { # change it to run model

        sing_cell_data <- sc_data_traj_inf()

        showNotification("Running trajectory inference model",
                         id = "tab3_m2",
                         duration = NULL)

        ## Extracts the dimension reduction info
        dimred <- sing_cell_data@reductions$pca@cell.embeddings

        sc_meta <- sc_meta()
        sc_meta_cluster <-  sc_meta_cluster()
        object_expression <- object_expression()
        object_counts <- object_counts()

        # test if the user set the initial cluster
        if ( !is.na(input$traj_init_clusters) ) {

            start_cells <- sc_meta[sc_meta$seurat_clusters == input$traj_init_clusters, ]
            #
            # filter(sc_meta, seurat_clusters == input$traj_init_clusters) %>%
            # select(cells)

            start_cells <- as.character(start_cells$cells)

        } else {

            start_cells <- NULL

        }

        # test if the user set the end cluster
        if (!is.na(input$traj_end_clusters)) {

            end_cells <- sc_meta[sc_meta$seurat_clusters == input$traj_end_clusters, ]

            # filter(sc_meta, seurat_clusters == input$traj_end_clusters) %>%
            # select(cells)

            end_cells <- as.character(end_cells$cells)

        } else {

            end_cells <- NULL

        }

        if ( input$ti_select_models == 0 ) { # slingshot locally

            if (input$ti_sample_number == 1) { # there are multiple samples

                # if (input$ti_timepoint == 1) { # the samples are timepoints
                #
                #     # Extract the treatment info to so it can be used as input for dynverse timepoints or discreat groups
                #     sc_cells_time_vec <- sc_cells_time_vec()
                #
                #     if ( input$ti_timepoint_cont_disct == 0 )  { # the sample names must be a numeric vector.
                #
                #         priors <- list(
                #             start_id = start_cells,
                #             end_id = end_cells,
                #             groups_id = sc_meta_cluster,
                #             dimred = dimred,
                #             timecourse_discrete = sc_cells_time_vec
                #
                #         )
                #
                #     } else if ( input$ti_timepoint_cont_disct == 1) {
                #
                #         priors <- list(
                #             start_id = start_cells,
                #             end_id = end_cells,
                #             groups_id = sc_meta_cluster,
                #             dimred = dimred,
                #             timecourse_continuous = sc_cells_time_vec
                #         )
                #     }
                #
                # } else if (input$ti_timepoint == 0) { # the samples are NOT timepoints

                # sc_cells_time_vec <- sc_cells_time_vec()

                priors <- list(
                    start_id = start_cells,
                    end_id = end_cells,
                    groups_id = sc_meta_cluster,
                    dimred = dimred
                )
                #}

            } else if (input$ti_sample_number == 0) { # there are only one sample

                priors <- list(
                    start_id = start_cells,
                    end_id = end_cells,
                    groups_id = sc_meta_cluster,
                    dimred = dimred
                )

            }

            sds <- run_fun_slingshot(expression = object_expression, priors = priors, verbose = T)

            removeNotification(id = "tab3_m2")

        } else if ( input$ti_select_models == 1 ) { # dynverse models - depends on docker

            dataset <- wrap_expression(
                counts = object_counts,
                expression = object_expression
            )

            # fetch newest version of the method
            method_id <- paste0("dynverse/ti_", req(input$ti_methods_list), ":latest")
            methods_selected <- create_ti_method_container(method_id)

            if ( input$ti_sample_number == 1 ) { # there are multiple samples

                # if ( input$ti_timepoint == 1 ) { # the samples are timepoints
                #
                #     sc_cells_time_vec <- sc_cells_time_vec()
                #
                #     if ( input$ti_timepoint_cont_disct == 0 )  { # test if the treatment is a timecourse (continuos o discret). For that, the sample name must be a numeric vector.
                #
                #         dataset <- add_prior_information(
                #             dataset,
                #             start_id = start_cells,
                #             end_id = end_cells,
                #             groups_id = sc_meta_cluster,
                #             dimred = dimred,
                #             timecourse_discrete = sc_cells_time_vec
                #         )
                #
                #         sds <- infer_trajectory(dataset,
                #                                 methods_selected(),
                #                                 verbose = T,
                #                                 give_priors = c("start_id",
                #                                                 "groups_id",
                #                                                 "end_id",
                #                                                 "dimred",
                #                                                 "timecourse_discrete"))
                #
                #     } else if ( input$ti_timepoint_cont_disct == 1 ) {
                #
                #         dataset <- add_prior_information(
                #             dataset,
                #             start_id = start_cells,
                #             end_id = end_cells,
                #             groups_id = sc_meta_cluster,
                #             dimred = dimred,
                #             timecourse_continuous = sc_cells_time_vec)
                #
                #         sds <- infer_trajectory(dataset,
                #                                 methods_selected(),
                #                                 verbose = T,
                #                                 give_priors = c("start_id",
                #                                                 "groups_id",
                #                                                 "end_id",
                #                                                 "dimred",
                #                                                 "timecourse_continuous"))
                #     }
                #
                # } else if (input$ti_timepoint == 0) { # the samples are NOT timepoints

                dataset <- add_prior_information(

                    dataset,
                    start_id = start_cells,
                    end_id = end_cells,
                    groups_id = sc_meta_cluster,
                    dimred = dimred
                )

                sds <- infer_trajectory(dataset,
                                        methods_selected(),
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
                                        methods_selected(),
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

        model <- model()
        sc_meta_cluster <- sc_meta_cluster()

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

                sc_cells_time_vec <- sc_cells_time_vec()

                plot_dimred(model,
                            grouping = sc_cells_time_vec,
                            color_density = "grouping") +
                    ggtitle("Cell grouping")#,

            }
        }

    })

    output$ti_traject <- renderPlot({

        model <- model()
        sc_meta_cluster <- sc_meta_cluster()

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
                sc_cells_time_vec <- sc_cells_time_vec()

                plot_dendro(model,
                            grouping = sc_cells_time_vec) +
                    ggtitle("Trajectory")#,

            }
        }

    })

    output$ti_graph <- renderPlot({

        model <- model()
        sc_meta_cluster <- sc_meta_cluster()

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
                sc_cells_time_vec <- sc_cells_time_vec()

                plot_graph(model,
                           grouping = sc_cells_time_vec,
                           expression_source = dataset) +
                    ggtitle("Trajectory represented as a graph")

            }
        }

    })

    output$p9_down <- downloadHandler(

        filename = function() {
            paste(input$dataset, ".", input$p9_format, sep = "")
        },
        content = function(file) {

            withProgress(message = "Please wait, preparing the data for download.",
                         value = 0.5, {

                             height <- as.numeric(input$p9_height)
                             width <- as.numeric(input$p9_width)
                             res <- as.numeric(input$p9_res)

                             if ( input$p9_down_opt == 0 ) {

                                 model <- model()
                                 sc_meta_cluster <- sc_meta_cluster()

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

                                         sc_cells_time_vec <- sc_cells_time_vec()

                                         p <-  plot_dimred(model,
                                                           grouping = sc_cells_time_vec,
                                                           color_density = "grouping") +
                                             ggtitle("Cell grouping")#,

                                     }
                                 }

                             } else if ( input$p9_down_opt == 1 ) {

                                 model <- model()
                                 sc_meta_cluster <- sc_meta_cluster()

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
                                         sc_cells_time_vec <- sc_cells_time_vec()

                                         p <-  plot_dendro(model,
                                                           grouping = sc_cells_time_vec) +
                                             ggtitle("Trajectory")#,

                                     }
                                 }

                             } else if ( input$p9_down_opt == 2 ) {

                                 model <- model()
                                 sc_meta_cluster <- sc_meta_cluster()

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

                                         sc_cells_time_vec <- sc_cells_time_vec()

                                         p <- plot_graph(model,
                                                         grouping = sc_cells_time_vec,
                                                         expression_source = dataset) +
                                             ggtitle("Trajectory represented as a graph")

                                     }
                                 }

                             }

                             ggsave(file,
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

        sds <- model()
        sds$lineages

    })


    features_tab3 <- reactive({

        inFile <- input$markers_list_tab3

        if (is.null(inFile))
            return(NULL)

        features <- read.csv(inFile$datapath, header = F)

        if ( ncol(features) == 2 ) {

            features <- features %>%
                dplyr::rename(GeneID = V1,
                              Group = V2)

        } else if ( ncol(features) > 2 ) {

            features <- features %>%
                dplyr::rename(GeneID = V1,
                              Group = V2,
                              Name = V3)

        }

        features

    })

    # Load the file and offers the parameters for heatmap
    observeEvent(input$load_markers_tab3, {

        # Painel that will apper after loading the list of markers containing the filter options
        output$marker_group_selec_tab3 = renderUI({

            groups_names <- features_tab3()
            groups_names <- unique(groups_names$Group)

            div(class = "option-group",
                pickerInput(
                    inputId = "features_group_tab3",
                    label = "Select the group of markers to test",
                    choices = sort( as.character(groups_names) ),
                    multiple = TRUE,
                    options = list(`actions-box` = TRUE)
                ))
        })

        output$marker_filter_genes_q_tab3 = renderUI({

            div(class = "option-group",
                radioButtons("filter_genes_q_tab3",
                             "Visualization options",
                             choices = list("Show all genes" = 1,
                                            "Select genes to show" = 0),
                             selected = 1))
        })

        # Define if using names of v5 or common names
        output$marker_genes_ids_tab3 = renderUI({

            div(class = "option-group",
                radioButtons("genes_ids_tab3",
                             "ID options",
                             choices = list("Use IDs" = "v5",
                                            "Use Name" = "name"),
                             selected = "v5"))
        })

        output$marker_genes_selec_tab3 = renderUI({

            groups_names <- features_tab3()
            genes_names <- dplyr::filter(groups_names,
                                         Group %in% req(input$features_group_tab3))

            if (input$genes_ids_tab3 == "v5") {

                genes_names <- unique(genes_names$GeneID)

            } else if (input$genes_ids_tab3 == "name") {

                genes_names <- unique(genes_names$Name)

            }

            genes_names <- genes_names[!is.na(genes_names)]

            # div(class = "option-group",
            #     selectizeInput('selected_genes',
            #                    'Select the genes to show',
            #                    choices = sort(genes_names),
            #                    selected = NULL,
            #                    multiple = T))

            div(class = "option-group",
                pickerInput(
                    inputId = "selected_genes_tab3",
                    label = "Select the genes to show",
                    choices = sort(as.character(genes_names)),
                    multiple = TRUE,
                    options = list(`actions-box` = TRUE)
                ))

        })

    })

    # Filter accordly with the parameters selected by the user
    filt_features_tab3 <- eventReactive(input$run_heatmap_tab3, {

        features_f <- features_tab3()

        features_f <- dplyr::filter(features_f,
                                    Group %in% req(input$features_group_tab3))

        #    features_f <- features_f[paste0("gene:", features_f$GeneID) %in% expressed_genes(), ]
        features_f <- features_f[ features_f$GeneID %in% expressed_genes(), ]

        if (input$filter_genes_q_tab3 == 0) {

            if (input$genes_ids_tab3 == "v5") {

                features_f <- features_f[features_f$GeneID %in% req(input$selected_genes_tab3), ]

            } else if (input$genes_ids_tab3 == "name") {

                features_f <- features_f[features_f$Name %in% req(input$selected_genes_tab3), ]

            }

        }

        features_f <- unique(features_f)
        features_f
    })

    dynverse_genes_list <- eventReactive(input$dynverse_def_imp_genes, {

        showNotification("Defining the most relevant genes",
                         id = "tab3_m5",
                         duration = NULL)

        if ( input$dynverse_opt == 0 ) { # global

            feat_importances <- dynfeature::calculate_overall_feature_importance(model(),
                                                                                 expression_source = dataset_inf())


        } else if ( input$dynverse_opt == 1 ) { # branch/lineage

            feat_importances <- dynfeature::calculate_branch_feature_importance(model(),
                                                                                expression_source = dataset_inf())

        } else if ( input$dynverse_opt == 2 ) { # bifurcation

            if ( is.na(input$branching_milestone) ) {

                feat_importances <- calculate_branching_point_feature_importance(model(),
                                                                                 expression_source = dataset_inf())

            } else {

                feat_importances <- calculate_branching_point_feature_importance(model(),
                                                                                 milestones_oi = input$branching_milestone,
                                                                                 expression_source = dataset_inf())

            }

        }

        removeNotification(id = "tab3_m5")

        showNotification("Select the number of genes to show in the heatmap and trajectory plots",
                         id = "tab3_m6",
                         duration = 30)

        feat_importances

    })

    dynverse_genes_list_filt <- eventReactive( c(input$dynverse_def_imp_genes, input$run_heatmap_tab3_dynverse), {

        features_imp <- dynverse_genes_list()

        if ( input$dynverse_opt == 0 ) { # global

            features <- features_imp %>%
                top_n(input$dynverse_n_genes, importance)
            # %>%
            #     separate("feature_id", sep = "\\:", into = c("a", "b") )

        } else if ( input$dynverse_opt == 1 ) { # branch/lineage

            if ( is.na(req(input$dynverse_branch_from)) ) {

                features <- features_imp %>%
                    filter(to == req(input$dynverse_branch_to)) %>%
                    top_n(input$dynverse_n_genes, importance)
                # %>%
                #     separate("feature_id", sep = "\\:", into = c("a", "b") )

            } else if ( is.na(req(input$dynverse_branch_to)) ) {

                features <- features_imp %>%
                    filter( from  == req(input$dynverse_branch_from)) %>%
                    top_n(input$dynverse_n_genes, importance)
                # %>%
                #     separate("feature_id", sep = "\\:", into = c("a", "b") )

            } else {

                features <- features_imp %>%
                    filter( from  == req(input$dynverse_branch_from) & to == req(input$dynverse_branch_to)) %>%
                    top_n(input$dynverse_n_genes, importance)
                #
                # %>%
                #     separate("feature_id", sep = "\\:", into = c("a", "b") )

            }

        } else if ( input$dynverse_opt == 2 ) { # bifurcation

            features <- features_imp %>%
                top_n(input$dynverse_n_genes, importance)
            # %>%
            #     separate("feature_id", sep = "\\:", into = c("a", "b") )


        }

        # features <- features$b

        features
    })

    #tradseq_genes_list <- eventReactive(input$run_heatmap_tab3_tradseq, {})

    ##################################################
    ### Expression plots - Using markers as input ####
    ##################################################

    observeEvent( c(input$dynverse_def_imp_genes, input$run_heatmap_tab3), {

        # This calculate the height of the plot based on the n of genes
        heatmap_n_genes_tab3 <- reactive({

            # Test what input to use for the heatmap
            if (input$ti_expre_opt == 0) {

                features <- filt_features_tab3()
                features <- unique(features$GeneID)

            } else if (input$ti_expre_opt == 1) {

                features <- dynverse_genes_list_filt()
                features <- features$feature_id

            } else if (input$ti_expre_opt == 2) {

                features <- tradseq_genes_list()

            }

            as.numeric(length(features))
        })
        heatmap_Height_tab3 <- reactive( 250 + ( 20 * heatmap_n_genes_tab3() ) )

        # Generates the heatmap plot
        output$heat_map_tab3 <- renderPlot({

            showNotification("Generating the heatmap",
                             id = "tab3_m4",
                             duration = NULL)

            # Test what input to use for the heatmap
            if (input$ti_expre_opt == 0) {

                features <- filt_features_tab3()
                features <- unique(features$GeneID)

            } else if (input$ti_expre_opt == 1) {

                features <- dynverse_genes_list_filt()
                features <- features$feature_id

            } else if (input$ti_expre_opt == 2) {

                features <- tradseq_genes_list()

            }


            #    features <- paste0("gene:", as.character(features))
            features <- as.character(features)

            ### Drawing heat map ###
            heatmap <- plot_heatmap(
                model(),
                expression_source = dataset_inf(),
                grouping = sc_meta_cluster(),
                features_oi = features
            )

            removeNotification(id = "tab3_m4")

            heatmap

        }, height = heatmap_Height_tab3() )

        output$heat_map_ui_tab3 <- renderUI({

            plotOutput( "heat_map_tab3", height = heatmap_Height_tab3() )

        })

        output$marker_to_feature_plot_tab3 = renderUI({

            # Test what input to use for the heatmap
            if (input$ti_expre_opt == 0) {

                features <- filt_features_tab3()
                features <- unique(features$GeneID)

            } else if (input$ti_expre_opt == 1) {

                features <- dynverse_genes_list_filt()
                features <- features$feature_id

            } else if (input$ti_expre_opt == 2) {

                features <- tradseq_genes_list()

            }

            div(class = "option-group",
                pickerInput(
                    inputId = "selected_genes_for_feature_plot_tab3",
                    label = "Select the genes to feature plots",
                    choices = sort(as.character(features)),
                    multiple = TRUE,
                    options = list(`actions-box` = TRUE)
                ))

        })

        output$p10_down <- downloadHandler(

            filename = function() {
                paste(input$dataset, ".", input$p10_format, sep = "")
            },
            content = function(file) {

                withProgress(message = "Please wait, preparing the data for download.",
                             value = 0.5, {

                                 height <- as.numeric(input$p10_height)
                                 width <- as.numeric(input$p10_width)
                                 res <- as.numeric(input$p10_res)

                                 # Test what input to use for the heatmap
                                 if (input$ti_expre_opt == 0) {

                                     features <- filt_features_tab3()
                                     features <- unique(features$GeneID)

                                 } else if (input$ti_expre_opt == 1) {

                                     features <- dynverse_genes_list_filt()
                                     features <- features$feature_id

                                 } else if (input$ti_expre_opt == 2) {

                                     features <- tradseq_genes_list()

                                 }

                                 #    features <- paste0("gene:", as.character(features))
                                 features <- as.character(features)

                                 ### Drawing heat map ###
                                 p <- plot_heatmap(
                                     model(),
                                     expression_source = dataset_inf(),
                                     grouping = sc_meta_cluster(),
                                     features_oi = features
                                 )

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

    #######3 >>>>>>>>>sc_meta_cluster
    output$dynverse_branch_from_ui <- renderUI({

        sc_meta_cluster <- sc_meta_cluster()
        clust <- unique(as.character(sc_meta_cluster$group_id))

        div(class = "option-group",
            pickerInput(
                inputId = "dynverse_branch_from",
                label = "Set the start branch (from)",
                choices = sort(as.character(clust)),
                multiple = F,
                options = list(`actions-box` = TRUE)
            ))

    })

    output$dynverse_branch_to_ui <- renderUI({

        sc_meta_cluster <- sc_meta_cluster()

        clust <- unique(as.character(sc_meta_cluster$group_id))
        clust <- clust[!clust %in% req(input$dynverse_branch_from)]

        div(class = "option-group",
            pickerInput(
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

            model <- model()

            local({

                my_i <- i

                plotname <- paste("ti_order_exp", my_i, sep="")
                output[[plotname]] <- renderPlot({

                    plot_dimred(model,
                                #    feature_oi = paste0("gene:", features_to_plot[my_i]),
                                feature_oi = features_to_plot[my_i],
                                expression_source = dataset_inf()) +
                        ggtitle("Cell grouping")

                })

                plotname <- paste("ti_traject_expr", my_i, sep="")
                output[[plotname]] <- renderPlot({

                    plot_dendro(model,
                                #    feature_oi = paste0("gene:", features_to_plot[my_i]),
                                feature_oi = features_to_plot[my_i],
                                expression_source = dataset_inf()) +
                        ggtitle("Trajectory")#,

                })

                plotname <- paste("ti_graph_expr", my_i, sep="")
                output[[plotname]] <- renderPlot({

                    plot_graph(model,
                               #   feature_oi = paste0("gene:", features_to_plot[my_i]),
                               feature_oi = features_to_plot[my_i],
                               expression_source = dataset_inf()) +
                        ggtitle("Trajectory represented as a graph")#,

                })

            })

        }

        # removeNotification("tab3_m7")

    })


    output$select_genes_add_plot_to_down_tab3_ui <- renderUI({

        filt_features <- req(input$selected_genes_for_feature_plot_tab3)

        div(class = "option-group",
            pickerInput(
                inputId = "select_genes_add_plot_to_down_tab3",
                label = "Select the genes that you want to download",
                choices = sort(as.character(filt_features)),
                multiple = TRUE,
                options = list(`actions-box` = TRUE)
            ))

    })

    observeEvent(input$start_down_add_plots_tab3, {


        withProgress(message = "Please wait, preparing the data for download.",
                     value = 0.5, {


                         if (file.exists("./images") == F ) {
                             dir.create('./images')
                         }

                         # Creates new folder to keep the plots organized. Adding date and time helps to avoid overwriting the data
                         path_new <- paste0("./images/trajectories_plots", format(Sys.time(),'_%Y%-%m%-%d%__%H%M%S'))

                         dir.create(path_new)
                         dir.create( paste0(path_new,"/Dimension_reduction") )
                         dir.create( paste0(path_new,"/Dendrogram") )
                         dir.create( paste0(path_new, "/Graph") )

                         model <- model()

                         genes <- req(input$select_genes_add_plot_to_down_tab3)

                         for( i in 1:length(genes) ){

                             # Saves the feature plots

                             p <- plot_dimred(model,
                                              feature_oi = genes[i],
                                              expression_source = dataset_inf()) +
                                 ggtitle("Cell grouping")

                             file <- paste0(path_new, "/Dimension_reduction/", genes[i], ".", input$add_p_tab3_feat_format)

                             ggsave(file,
                                    p,
                                    height=input$add_p_tab3_feat_height,
                                    width=input$add_p_tab3_feat_width,
                                    units="cm",
                                    dpi=as.numeric(input$add_p_tab3_feat_res))

                             # Saves the violin plots

                             file <- paste0(path_new, "/Dendrogram/", genes[i], ".", input$add_p_tab3_violin_format)

                             p <- plot_dendro(model,
                                              feature_oi = genes[i],
                                              expression_source = dataset_inf()) +
                                 ggtitle("Trajectory")

                             ggsave(file,
                                    p,
                                    height=input$add_p_tab3_violin_height,
                                    width=input$add_p_tab3_violin_width,
                                    units="cm",
                                    dpi=as.numeric(input$add_p_tab3_violin_res))

                             # Saves the dot plots

                             file <- paste0(path_new, "/Graph/", genes[i], ".", input$add_p_tab3_dot_format)

                             p <- plot_graph(model,
                                             feature_oi = genes[i],
                                             expression_source = dataset_inf()) +
                                 ggtitle("Trajectory represented as a graph")#,

                             ggsave(file,
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

                features <- filt_features_tab3()
                features <- unique(features$GeneID)

            } else if (input$ti_expre_opt == 1) {

                features <- dynverse_genes_list_filt()
                features <- features$feature_id

            }
            # else if (input$ti_expre_opt == 2) {
            #
            #     #features <- tradseq_genes_list()
            #
            # }

            as.numeric(length(features))
        })
        heatmap_Height_tab3_dynv <- reactive( 250 + ( 20 * heatmap_n_genes_tab3_dynv() ) )

        # Generates the heatmap plot
        output$heat_map_tab3_dynv <- renderPlot({

            showNotification("Generating the heatmap",
                             id = "tab3_m8",
                             duration = NULL)

            # Test what input to use for the heatmap
            if (input$ti_expre_opt == 0) {

                features <- filt_features_tab3()
                features <- unique(features$GeneID)

            } else if (input$ti_expre_opt == 1) {

                features <- dynverse_genes_list_filt()
                features <- features$feature_id

            }
            # else if (input$ti_expre_opt == 2) {
            #
            #     features <- tradseq_genes_list()
            #
            # }

            #    features <- paste0("gene:", as.character(features))
            features <- as.character(features)

            ### Drawing heat map ###
            heatmap <- plot_heatmap(
                model(),
                expression_source = dataset_inf(),
                grouping = sc_meta_cluster(),
                features_oi = features
            )

            removeNotification(id = "tab3_m8")

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
                pickerInput(
                    inputId = "selected_genes_for_feature_plot_tab3_dynv",
                    label = "Select the genes to feature plots",
                    choices = sort(as.character(features)),
                    multiple = TRUE,
                    options = list(`actions-box` = TRUE)
                ))

        })

        output$p11_down <- downloadHandler(

            filename = function() {
                paste(input$dataset, ".", input$p11_format, sep = "")
            },
            content = function(file) {

                withProgress(message = "Please wait, preparing the data for download.",
                             value = 0.5, {

                                 height <- as.numeric(input$p11_height)
                                 width <- as.numeric(input$p11_width)
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

    observeEvent(input$run_feature_plot_tab3_dynv, {

        showNotification("Generating the expression plots",
                         id = "tab3_m8",
                         duration = 30)

        output$ti_order_express_dynv <- renderUI({

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
                                #    feature_oi = paste0("gene:", features_to_plot[my_i]),
                                feature_oi = features_to_plot[my_i],
                                expression_source = dataset_inf()) +
                        ggtitle("Cell grouping")

                })

                plotname <- paste("ti_traject_expr_dynv", my_i, sep="")
                output[[plotname]] <- renderPlot({

                    plot_dendro(model,
                                #    feature_oi = paste0("gene:", features_to_plot[my_i]),
                                feature_oi = features_to_plot[my_i],
                                expression_source = dataset_inf()) +
                        ggtitle("Trajectory")#,

                })

                plotname <- paste("ti_graph_expr_dynv", my_i, sep="")
                output[[plotname]] <- renderPlot({

                    plot_graph(model,
                               #   feature_oi = paste0("gene:", features_to_plot[my_i]),
                               feature_oi = features_to_plot[my_i],
                               expression_source = dataset_inf()) +
                        ggtitle("Trajectory represented as a graph")#,

                })

            })

        }

        #removeNotification("tab3_m8")

    })

    output$select_genes_add_plot_to_down_tab3_ui_dynv <- renderUI({

        filt_features <- input$selected_genes_for_feature_plot_tab3_dynv

        div(class = "option-group",
            pickerInput(
                inputId = "select_genes_add_plot_to_down_tab3_dynv",
                label = "Select the genes that you want to download",
                choices = sort(as.character(filt_features)),
                multiple = TRUE,
                options = list(`actions-box` = TRUE)
            ))

    })

    observeEvent(input$start_down_add_plots_tab3_dynv, {

        withProgress(message = "Please wait, preparing the data for download.",
                     value = 0.5, {

                         if (file.exists("./images") == F ) {
                             dir.create('./images')
                         }

                         # Creates new folder to keep the plots organized. Adding date and time helps to avoid overwriting the data
                         path_new <- paste0("./images/trajectories_plots", format(Sys.time(),'_%Y%-%m%-%d%__%H%M%S'))

                         dir.create(path_new)
                         dir.create( paste0(path_new,"/Dimension_reduction") )
                         dir.create( paste0(path_new,"/Dendrogram") )
                         dir.create( paste0(path_new, "/Graph") )

                         model <- model()

                         genes <- req(input$select_genes_add_plot_to_down_tab3_dynv)

                         for( i in 1:length(genes) ){

                             # Saves the feature plots

                             p <- plot_dimred(model,
                                              feature_oi = genes[i],
                                              expression_source = dataset_inf()) +
                                 ggtitle("Cell grouping")

                             file <- paste0(path_new, "/Dimension_reduction/", genes[i], ".", input$add_p_tab3_feat_format_dynv)

                             ggsave(file,
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

                             ggsave(file,
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

                             ggsave(file,
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
            paste(input$dataset, ".csv", sep = "")
        },
        content = function(file) {

            withProgress(message = "Please wait, preparing the data for download.",
                         value = 0.5, {

                             #genes <- dynverse_genes_list()

                             write.csv(dynverse_genes_list(), file, row.names = FALSE)

                         })

            #removeNotification(id = "D1")
        }

    )

    output$download_dynverse_genes_filt <- downloadHandler(

        filename = function() {
            paste(input$dataset, ".csv", sep = "")
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
        pickerInput("dbselection",
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
        pickerInput("datasetselection",
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

        pickerInput("filterselection",
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
        pickerInput("attpageselection", "Select BioMart attributes Page:",
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

        # render pickerInput
        pickerInput("attselection",
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
        pickerInput("geneselection",
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
                             
                             ggsave(
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

}
