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
