## Adapted from https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/

modify_vlnplot<- function(obj,
                          feature,
                          pt.size = pt.size,
                          flip_opt = T,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
    p<- VlnPlot(obj,
                features = feature,
                flip = flip_opt,
                pt.size = pt.size,
                group.by = "seurat_clusters2",
                ... )  +
        xlab("") + ylab(feature) + ggtitle("") +
        theme(legend.position = "none",
              plot.title= element_blank(),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.y = element_text(size = rel(1), angle = 0),
              axis.text.y = element_text(size = rel(1)),
              plot.margin = plot.margin)

    return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
    ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
    return(ceiling(ymax))
}

## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = pt.size,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {

    plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, pt.size = pt.size, ...))

    plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
        theme(axis.text.x=element_text(), axis.ticks.x = element_line())

    # change the y-axis tick to only max value
    ymaxs<- purrr::map_dbl(plot_list, extract_max)
    plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x +
                                scale_y_continuous(breaks = c(y)) +
                                expand_limits(y = y))

    p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
    return(p)
}

stacked_violin <- function(name,
                           rds_file = rds_file,
                           selected_genes = selected_genes,
                           pt.size = 0){

    StackedVlnPlot(obj = rds_file,
                   features = selected_genes,
                   pt.size = pt.size,
                   assay = "RNA") # Set 0 if you want to remove the dots
    

}