plot.ssGSA <- function(x, orderby, whplot = 1, col = "black", samplename =FALSE,
                maintitle = "", ssgsaInfo = TRUE, cex.sub = 0.8, ...){
    if(!inherits(x, "ssGSA")) stop("x must be of class ssGSA")

    if(is(orderby, 'factor'))
        orderby2  <- as.numeric(orderby)
    else
        orderby2  <- orderby

    wh2sum  <- diff(range(orderby2))
    ru  <- runif(length(orderby2), -0.3, 0.3)

    x$stats  <- x$stats[rownames(x$res),]
    if(is.character(whplot))
        whplot  <- which(rownames(x$res) %in% whplot)

    if(maintitle == ""){
        maintitle <- rownames(x$stats)[whplot[1]]
        if(ssgsaInfo)
            subtitle  <- paste(names(x$res[whplot,]),
                            round( as.numeric(x$res[whplot,]),3),
                            sep=" = ",collapse="; " )
    }
    else if(ssgsaInfo) subtitle <- ""

    if(samplename)
        {
            plot(orderby2 + ru, x$stats[whplot,], col = col, pch = "",
                xlim = c(min(orderby2)-wh2sum/3, max(orderby2)+wh2sum/3),
                xaxt = "n", main = maintitle, ...)
            text(orderby2 + ru, x$stats[whplot,], label = colnames(x$stats),
                col = col, ...)
        }
    else{
        plot(orderby2 + ru, as.numeric(x$stats[whplot,]), col = col,
            xlim = c(min(orderby2)-wh2sum/3, max(orderby2)+wh2sum/3),
            xaxt = "n", main = maintitle,...)
    }
    if(ssgsaInfo)
        mtext(subtitle, 3, line=0.5, col= "darkgray", cex = cex.sub)
    mtext(levels(orderby), side=1, at = unique(orderby2))
}


