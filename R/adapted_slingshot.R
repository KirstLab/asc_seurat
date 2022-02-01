run_fun_slingshot <- function(expression, priors, verbose, seed) {
    
    start_id <- priors$start_id
    end_id <- priors$end_id
    dimred <- priors$dimred
    groups_id <- priors$groups_id
    
    #####################################
    ###        INFER TRAJECTORY       ###
    #####################################
    
    start_cell <- if (!is.null(start_id)) { sample(start_id, 1) } else { NULL }
    
    # TIMING: done with preproc
    checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))
    
    rd <- dimred
    labels <- groups_id %>% deframe()

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
        reducedDim = 'PCA',
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
            cell_ids = base::rownames(expression)
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
