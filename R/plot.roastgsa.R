##############################################
#########  function:  plot.roastgsa  #########
##############################################
plot.roastgsa <- function(x, type = c("stats","GSEA"), whplot = 1,
        maintitle = "",gsainfo = TRUE, cex.sub = 0.8, lwd = 2, ...){

    if(!inherits(x, "roastgsa")) stop("x must be of class roastgsa")
    obj  <- x
    if("stats" %in% type)
        plotStats(obj, whplot = whplot, maintitle = maintitle,
            gsainfo = gsainfo, cex.sub = cex.sub, lwd = lwd, ...)
    if("GSEA" %in% type)
        plotGSEA(obj, whplot = whplot, maintitle = maintitle,
            gsainfo = gsainfo, cex.sub = cex.sub, lwd = lwd, ...)
}
