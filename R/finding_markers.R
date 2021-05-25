finding_markers <- function(name,
                            sc_data = sc_data,
                            find_markers_tab1_opt = find_markers_tab1_opt,
                            find_markers_tab1_return.thresh = find_markers_tab1_return.thresh,
                            find_markers_tab1_logfc.threshold = find_markers_tab1_logfc.threshold,
                            find_markers_tab1_min.pct = find_markers_tab1_min.pct,
                            find_markers_tab1_test.use = find_markers_tab1_test.use,
                            assay_choice = assay_choice,
                            ident.1 = ident.1,
                            ident.2 = ident.2,
                            find_markers_tab1_filt_pos = find_markers_tab1_filt_pos
) {
    
    if ( find_markers_tab1_opt == 0 ) { #all clusters
        
        if ( is.na(find_markers_tab1_return.thresh) ) {
            
            markers_tab1 <- FindAllMarkers(sc_data,
                                           assay = assay_choice,
                                           logfc.threshold = find_markers_tab1_logfc.threshold,
                                           min.pct = find_markers_tab1_min.pct,
                                           test.use = find_markers_tab1_test.use
            )
            
        } else {
            
            markers_tab1 <- FindAllMarkers(sc_data,
                                           assay = assay_choice,
                                           logfc.threshold = find_markers_tab1_logfc.threshold,
                                           min.pct = find_markers_tab1_min.pct,
                                           test.use = find_markers_tab1_test.use,
                                           return.thresh = find_markers_tab1_return.thresh
            )
            
        }
        
    } else if ( find_markers_tab1_opt == 1 ) { # one specific cluster
        
        if ( is.na(find_markers_tab1_return.thresh) ) {
            
            markers_tab1 <- FindMarkers(sc_data,
                                        ident.1 = ident.1,
                                        ident.2 = NULL,
                                        logfc.threshold = find_markers_tab1_logfc.threshold,
                                        min.pct = find_markers_tab1_min.pct,
                                        test.use = find_markers_tab1_test.use,
                                        only.pos = find_markers_tab1_filt_pos,
                                        assay = assay_choice)
            
        } else {
            
            markers_tab1 <- FindMarkers(sc_data,
                                        ident.1 = ident.1,
                                        ident.2 = NULL,
                                        logfc.threshold = find_markers_tab1_logfc.threshold,
                                        min.pct = find_markers_tab1_min.pct,
                                        test.use = find_markers_tab1_test.use,
                                        only.pos = find_markers_tab1_filt_pos,
                                        return.thresh = find_markers_tab1_return.thresh,
                                        assay = assay_choice)
            
        }
        
    } else if ( find_markers_tab1_opt == 2 ) { # distinguishing a cluster from other(s) cluster(s)
        
        if ( is.na(find_markers_tab1_return.thresh) ) {
            
            markers_tab1 <- FindMarkers(sc_data,
                                        ident.1 = ident.1,
                                        ident.2 = ident.2,
                                        logfc.threshold = input$find_markers_tab1_logfc.threshold,
                                        min.pct = find_markers_tab1_min.pct,
                                        test.use = find_markers_tab1_test.use,
                                        only.pos = find_markers_tab1_filt_pos,
                                        return.thresh = find_markers_tab1_return.thresh,
                                        assay = assay_choice)
            
        } else {
            
            markers_tab1 <- FindMarkers(sc_data,
                                        ident.1 = ident.1,
                                        ident.2 = ident.2,
                                        logfc.threshold = find_markers_tab1_logfc.threshold,
                                        min.pct = find_markers_tab1_min.pct,
                                        test.use = find_markers_tab1_test.use,
                                        only.pos = find_markers_tab1_filt_pos,
                                        return.thresh = find_markers_tab1_return.thresh,
                                        assay = assay_choice)
        }
        
    }
    
    
    if ( find_markers_tab1_opt == 0 ) { 
        
        markers_tab1 <- markers_tab1 %>% 
            dplyr::rename(geneID = gene) %>%
            relocate("geneID")
        
    } else {
        
        markers_tab1 <- rownames_to_column(markers_tab1,
                                           var = "geneID")
    }
    
    markers_tab1
    
}