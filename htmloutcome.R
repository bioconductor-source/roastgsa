##############################################
#########  function:  html agsea #############
##############################################
htmlagsea <- function(obj, htmlpath = "", htmlname = "file.html", plotpath ="",
                      plotmean = TRUE, plotgsea = TRUE, indheatmap = TRUE, ploteffsize = TRUE,
                      y, whplots = NULL, tit = "", margins = c(5,16), geneDEhtmlfiles  = NULL,
                      sizesHeatmap = c(1200, 800), varrot, intvar, adj.var = NULL, mycol, ... ){

      if(!inherits(obj,"roastgsa")) stop("not a roastgsa object")
      if(ploteffsize) if(missing(varrot)) stop("varrot is missing")

      x <- data.frame(geneset  = rownames(obj$res),obj$res)
      index <- obj$index[rownames(x)]

       if(plotmean | plotgsea |  indheatmap){
        dir.create(paste0(htmlpath,plotpath))
        if(is.null(whplots)) whplots <- names(index)
         if(!is.na(whplots[1])){
           stats <- sort(obj$stats)
           index <- sapply(obj$index, function(x) which(names(stats)%in%x))

           for(k in whplots){

                if(plotmean){
                    png(paste0(htmlpath, plotpath, gsub("[[:punct:]]"," ",k),"_mean.png"))
                    plotMean(obj, whplot = k, ...)
                    dev.off()
                }
                if(plotgsea){
                    png(paste0(htmlpath,plotpath, gsub("[[:punct:]]"," ",k),"_gsea.png"))
                    plotGSEA(obj, whplot = k, ...)
                    dev.off()
                }
                if(indheatmap){
                    png(paste0(htmlpath,plotpath, gsub("[[:punct:]]"," ",k),"_heatmap.png"),
                        width = sizesHeatmap[2], height = sizesHeatmap[1])
                    heatmapMean(obj, y, whplot = k, mycol =mycol, intvar, adj.var, ...)
                    dev.off()
                }
                if(ploteffsize){
                    png(paste0(htmlpath,plotpath, gsub("[[:punct:]]"," ",k),"_effsize.png"))
                    ploteffsignaturesize(obj, varrot,  whplot = k)
                    dev.off()
                }

           }
        }
       }
       if(!is.null(geneDEhtmlfiles)) x$geneDEinfo <- rep("view", dim(x)[1])



        if(plotmean) x$plot_mean <- NA
        if(plotgsea)  x$plot_gsea <- NA
        if(indheatmap) x$heatmap <- NA
        if(ploteffsize) x$plot_effsize <- NA

        links <- vector('list',length=ncol(x));
        names(links) <- colnames(x);
        plots <- links;

        if(plotmean)
            links$plot_mean <- plots$plot_mean <-  paste0(plotpath, gsub("[[:punct:]]"," ",rownames(x)),"_mean.png")
        if(plotgsea)
            links$plot_gsea <- plots$plot_gsea <- paste0(plotpath, gsub("[[:punct:]]"," ",rownames(x)),"_gsea.png")
        if(indheatmap)
            links$heatmap <- plots$heatmap <- paste0(plotpath, gsub("[[:punct:]]"," ",rownames(x)),"_heatmap.png")
        if(ploteffsize)
            links$plot_effsize <- plots$plot_effsize  <-
                paste0(plotpath, gsub("[[:punct:]]"," ",rownames(x)),"_effsize.png")

        if(!is.null(geneDEhtmlfiles))
        {
                links$geneDEinfo <- geneDEhtmlfiles
                if(dim(x)[1]>200) links$geneDEinfo[201:dim(x)[1]] <- NA
        }

        write.html.mod(x, file=paste0(htmlpath, htmlname), links=links, tiny.pic=plots,  title=tit)


}

