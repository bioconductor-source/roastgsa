##############################################
#########  function:  plotstats    ###########
##############################################

plotStats <- function(obj, whplot = 1, maintitle = "", statistic = "mean",
                ylimAll = TRUE, ylim = NULL,  minpointsDens = 20,
                gsainfo = TRUE, cex.sub = 0.8, lwd = 2, ...) {
    res <- obj$res
    index <- obj$index[rownames(res[whplot,])]
    stats <- obj$stats
    modt2 <- stats[order(-stats)]

    if(is.list(index))
        index2 <- lapply(index, function(o) which(names(modt2)%in%o))
    else
        index2 <-  list(which(names(modt2)%in%index))
    def.par <- par(no.readonly = TRUE)

    if(obj$self.contained)   stt <- modt2
    else{
        if(statistic == "median") stt <- (modt2 - median(stats))/mad(stats)
        else stt <- (modt2 - mean(stats))/sd(stats)
    }
    es <- lapply(index2, function(o) stt[o])

    for(K in names(index2)){
        if(maintitle == ""){
            maintitle <- K
            if(gsainfo)
                subtitle  <- paste(names(res[K,]),
                    round( as.numeric(res[K,]),3),sep=" = ",collapse="; " )
        }
        else if(gsainfo) subtitle <- ""

        le <- length(es[[K]])
        if(is.null(ylim)){
            if(ylimAll) ylim <- c(min(stt),max(stt))
            else ylim <- c(min(es[[K]]), max(es[[K]]))
        }
        s <- rep(FALSE, length(stt)); s[index2[[K]]] <- TRUE
        if(sum(s) > minpointsDens)
            layout(t(as.matrix(c(1,2,3))),heights=c(4),widths=c(4,.5,1))
        else layout(t(as.matrix(c(1,2))),heights=c(4),widths=c(4,1))

        par(mar=c(5,5,5,1))
        plot(seq_len(le),es[[K]], type='l',axes=FALSE, ylab = "stat",
            xlab="ordered genes", oma=c(0,0,0,0),
            col='darkgreen',lwd=lwd, main = maintitle, ylim=ylim, ...)
        polygon( c(seq_len(le), le,0),c(es[[K]],0,0),col=4, border = 4)

        if(gsainfo)
            mtext(subtitle, 3, line=0.5, col= "darkgray", cex = cex.sub)

        axis(2)
        whchange <- ifelse(sum(es[[K]]<0)>0, which(es[[K]]<0)[1], le)
        axis(1, at = c(1,whchange,le))
        abline(v= whchange, col=1, lty=3)

        par(mar=c(5,1,5,1))
        s <- rep(FALSE, length(stt)); s[index2[[K]]] <- TRUE
        if (sum(s)<50) myCols <- "darkgreen"
        else myCols <- densCols(which(s[seq_len(length(s))]),
            colramp=colorRampPalette(brewer.pal(9, "Greens")[-c(seq_len(1))]))

        plot(NA,NA, ylim=ylim, xlim=c(0,1),axes=FALSE, xlab='',
            ylab='',sub="",col='green',lwd=lwd)
        abline(h=es[[K]],col=myCols)
        maxmin <- ylim[2] - ylim[1]
        polygon(c(0-0.05,0-0.05,1+0.05,1+0.05),
                c(ylim[1], ylim[1] + maxmin/15, ylim[1] + maxmin/15, ylim[1]),
                col= rgb(0, 0, 0.5,0.3), border = NA)
        polygon(c(0-0.05,0-0.05,1+0.05,1+0.05),
                c(ylim[2], ylim[2] - maxmin/15, ylim[2] - maxmin/15, ylim[2]),
                col= rgb(0.5, 0, 0,0.3), border = NA)

        text(0.5,ylim[1] + maxmin/30, "-", cex=1.5)
        text(0.5,ylim[2] - maxmin/30, "+", cex=1.5)

        if(sum(s)>20){
            par(mar=c(5,1,5,2))
            ds <- density(es[[K]])
            plot(ds$y,ds$x, ylim=ylim, type="l", axes=FALSE,
                oma=c(0,0,0,0), main="", ylab="", xlab ="", col="darkblue")
            abline(h=0, col=1,lty=3)
        }
    }
    par(def.par)
}
