##############################################
#########  function:  html rgsea #############
##############################################
htmlrgsa <- function(obj, htmlpath = "", htmlname = "file.html", plotpath ="",
    plotstats = TRUE, plotgsea = TRUE, indheatmap = TRUE, ploteffsize = TRUE,
    links_plots = list(stats=NULL, gsea=NULL, heatmap=NULL, effsize=NULL),
    y, whplots = NULL, geneDEhtmlfiles  = NULL, tit = "", margins = c(5,16),
    sizesHeatmap = c(1200, 800), typeheatmap = c("heatmap.2","ggplot2"), intvar,
    adj.var = NULL, mycol, varrot, psel = NULL, sorttable, dragtable,  ... ){

    if(!inherits(obj,"roastgsa")) stop("not a roastgsa object")
    if(ploteffsize) if(missing(varrot)) stop("varrot is missing")

    x <- data.frame(geneset  = rownames(obj$res),obj$res)
    index <- obj$index[rownames(x)]
    psel2  <- psel

    if(plotstats | plotgsea |  indheatmap){
        dir.create(paste0(htmlpath,plotpath))
        if(is.null(whplots)) whplots <- names(index)
        if(!is.na(whplots[1])){
            stats <- sort(obj$stats)
            index <- lapply(obj$index, function(x) which(names(stats)%in%x))
            for(k in whplots){
                if(plotstats){
                    png(paste0(htmlpath, plotpath,
                        gsub("[[:punct:]]"," ",k),"_stats.png"))
                    plotStats(obj, whplot = k, ...)
                    dev.off()
                }
                if(plotgsea){
                    png(paste0(htmlpath,plotpath, gsub("[[:punct:]]"," ",k),
                        "_gsea.png"))
                    plotGSEA(obj, whplot = k, ...)
                    dev.off()
                }
                if(indheatmap){
                    png(paste0(htmlpath,plotpath, gsub("[[:punct:]]"," ",k),
                        "_heatmap.png"),width = sizesHeatmap[2],
                        height = sizesHeatmap[1])
                    if(typeheatmap[1] =="ggplot2")
                        heatmaprgsa_hm(obj, y, whplot = k, mycol =mycol,
                            intvar =intvar, adj.var= adj.var, psel = psel2, ...)
                    else
                        heatmaprgsa(obj, y, whplot = k, mycol =mycol,
                            intvar =intvar, adj.var= adj.var, psel = psel2, ...)
                    dev.off()
                }
                if(ploteffsize){
                    png(paste0(htmlpath,plotpath, gsub("[[:punct:]]"," ",k),
                        "_effsize.png"))
                    ploteffsignaturesize(obj, varrot,  whplot = k)
                    dev.off()
                }
            }
        }
    }
    if(!is.null(geneDEhtmlfiles)) x$geneDEinfo <- rep("view", dim(x)[1])


    if(plotstats) x$plot_stats <- NA
    if(plotgsea)  x$plot_gsea <- NA
    if(indheatmap) x$heatmap <- NA
    if(ploteffsize) x$plot_effsize <- NA

    links <- vector('list',length=ncol(x));
    names(links) <- colnames(x);
    plots <- links;

    if(!is.null(links_plots$stats))
        links$plot_stats <- plots$plot_stats <- links_plots$stats
    else{
        if(plotstats)
            links$plot_stats <- plots$plot_stats <-
                paste0(plotpath, gsub("[[:punct:]]"," ",rownames(x)),
                    "_stats.png")
    }

    if(!is.null(links_plots$gsea))
        links$plot_gsea <- plots$plot_gsea <- links_plots$gsea
    else{
        if(plotgsea)
            links$plot_gsea <- plots$plot_gsea <-
                paste0(plotpath, gsub("[[:punct:]]"," ",rownames(x)),
                    "_gsea.png")
    }

    if(!is.null(links_plots$heatmap))
        links$heatmap <- plots$heatmap <- links_plots$heatmap
    else{
        if(indheatmap)
            links$heatmap <- plots$heatmap <-
                paste0(plotpath, gsub("[[:punct:]]"," ",rownames(x)),
                    "_heatmap.png")
    }

    if(!is.null(links_plots$heatmap))
        links$effsize <- plots$effsize <- links_plots$effsize
    else{
        if(ploteffsize)
            links$plot_effsize <- plots$plot_effsize  <-
                paste0(plotpath, gsub("[[:punct:]]"," ",rownames(x)),
                    "_effsize.png")
    }

    if(!is.null(geneDEhtmlfiles)){
        links$geneDEinfo <- rep(NA, dim(x)[1])
        links$geneDEinfo[seq_len(length(geneDEhtmlfiles))] <- geneDEhtmlfiles
    }

    write.html.mod(x, file=paste0(htmlpath, htmlname), links=links,
        tiny.pic=plots, title=tit, sorttable = sorttable, dragtable = dragtable)
}

