FeaturePlotSingle <- function(obj = obj,
                              feature = feature,
                              metadata_column = metadata_column,
                              assay_id = assay_id,
                              pt.size = pt.size,
                              order = order,
                              reduction = reduction,
                              label = label) {
    
    all_cells <- colnames(obj)
    groups <- unique(obj@meta.data[, metadata_column])
    
    # the minimal and maximal of the value to make the legend scale the same.
    minimal <- min(obj[[assay_id]]@data[feature, ])
    maximal <- max(obj[[assay_id]]@data[feature, ])
    
    ps <- list()
    
    for (group in groups) {
        
        subset_indx <- obj@meta.data[, metadata_column] == group
        subset_cells <- all_cells[subset_indx]
        
        p <- suppressMessages(Seurat::FeaturePlot(obj,
                                                  features = feature,
                                                  cells = subset_cells) +
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