
########################################################
#########  function:  plot adjheatmapMean  #############
########################################################

heatmapMean <- function(obj, y, intvar, adj.var = NULL, whplot = 1, toplot = TRUE, pathwaylevel = FALSE,
                           mycol = c("black","orange","green","white"), sample2zero = FALSE,
                           dendrogram =  "n", col=  bluered(100),trace='none',
                           margins=c(8,16),notecol='black',notecex=1, keysize=.9, cexCol=1.5,
                           Rowv = NULL, Colv = FALSE, las =2, fdrkey = FALSE, ...){
   # browser()
    if(is.numeric(whplot)) whplot <- rownames(obj$res)[whplot]

    if(missing(intvar)) intvar <- tail(colnames(attr(terms(as.formula(obj$form)), "factors")),1)
    whint <- which(colnames(obj$covar) == intvar)
    if(any(obj$contrast==-1)){
        wh1 <- rowSums(as.matrix(obj$design[,obj$contrast == 1])) > 0
        wh2 <- rowSums(as.matrix(obj$design[,obj$contrast == -1])) > 0
    }
    else{
        wh1 <- rowSums(as.matrix(obj$design[,obj$contrast == 1])) > 0
        wh2 <- rowSums(as.matrix(obj$design[,obj$contrast == 1])) <= 0
    }
    ##wh2 <- obj$covar[,whint] == levels(obj$covar[,whint])[2]
    ##wh1 <- obj$covar[,whint] == levels(obj$covar[,whint])[1]

    if(pathwaylevel){
    #    browser()
        if(length(whplot)<2) stop("whplot must be of lenght > 1")

        DAT4HM <- t(sapply(whplot, function(j){
            hm <- heatmapMean(obj, y, whplot = j, toplot = FALSE, sample2zero  = TRUE, adj.var = adj.var, ...)
            apply(hm,2,mean)
        }))

#        apply(DAT4HM[,which(wh1)+2],1,mean)-        apply(DAT4HM[,which(wh2)+2],1,mean)

        DAT4HM <- DAT4HM/(obj$res[whplot,"est"]/ obj$res[whplot,"nes"])
        colnames(DAT4HM)[1] <- "nes"
        DAT4HM[,1] <- obj$res[whplot,"nes"]
        if(fdrkey){
            DAT4HM2 <- DAT4HM
            rownames(DAT4HM) <- paste0(rownames(DAT4HM), formatApval(obj$res[whplot,]$adj.pval))
        }
    }
    else{
        lmf <- lmFit(y, obj$design)
        coefs <- lmf$coef
        coefs[is.na(coefs)] <- 0

        ##       if(ncol(obj$design)>2)
        if(!is.null(adj.var))
            yr <- y - coefs[, adj.var] %*% t(obj$design[,adj.var])
        else
            yr <- y

        ys <- yr[obj$index[[whplot]], ]

        modts <- obj$stats[obj$index[[whplot]]]
        sde <- coefs[names(modts), c(which(obj$contrast == 1), which(obj$contrast == -1))]

        if(class(sde)=="numeric") sde <- sde/modts
        if(is.matrix(sde))
            sde <- (sde[,1] - sde[,2])/modts

        if(obj$statistic == "maxmean")
            modts <- if(obj$res[whplot,"est"]>0) modts*(modts>0) else modts*(modts<0)

        if(!obj$self.contained ){
            if(obj$statistic == "maxmean"){
                if(obj$res[whplot,"est"]>0){
                    nuapos <- mean(obj$stats*(obj$stats > 0))
                    sdapos <- sd(obj$stats*(obj$stats > 0))
                    modts <-  (modts - nuapos)/sdapos
                }
                else{
                    nuaneg <- mean(-obj$stats*(obj$stats < 0))
                    sdaneg <- sd(obj$stats*(obj$stats < 0))
                    modts <-  (modts + nuaneg)/sdaneg
                }
            }
            else
                modts <- (modts - mean(obj$stats))/sd(obj$stats)
        }

        DAT4HM <- ys/sde

        if(!obj$self.contained){
            if(obj$statistic == "maxmean"){
                if(obj$res[whplot,"est"]>0){
                    modts[obj$stats[obj$index[[whplot]]]< 0] <- 0
                    if(sample2zero){
                        DAT4HM[obj$stats[obj$index[[whplot]]]< 0, ] <- 0
                        DAT4HM[,wh1] <-  (DAT4HM[,wh1] - nuapos/2)/sdapos
                        DAT4HM[,wh2] <-  (DAT4HM[,wh2] + nuapos/2)/sdapos
                        DAT4HM[,!(wh2 | wh1)] <-  (DAT4HM[,!(wh1 | wh2)])/sdapos
                    }
                    else
                        DAT4HM <- DAT4HM/sdapos
                }
                else{
                    modts[obj$stats[obj$index[[whplot]]] > 0] <- 0
                    if(sample2zero){
                        DAT4HM[obj$stats[obj$index[[whplot]]]> 0, ] <- 0
                        DAT4HM[,wh1] <-  (DAT4HM[,wh1] + nuaneg/2)/sdaneg
                        DAT4HM[,wh2] <-  (DAT4HM[,wh2] - nuaneg/2)/sdaneg
                        DAT4HM[,!(wh2 | wh1)] <-  (DAT4HM[,!(wh1 | wh2)])/sdaneg
                    }
                    else
                        DAT4HM <- DAT4HM/sdaneg
                }
            }
            else{
                if(sample2zero){
                    DAT4HM[,wh1] <-  (DAT4HM[,wh1] - mean(obj$stats)/2)/sd(obj$stats)
                    DAT4HM[,wh2] <-  (DAT4HM[,wh2] + mean(obj$stats)/2)/sd(obj$stats)
                }
                else
                    DAT4HM <- DAT4HM/sd(obj$stats)
            }
        }
##            browser()
        DAT4HM <- t(apply(DAT4HM, 1, function(x) (x - mean(x))))

        sdr <- apply(DAT4HM, 1, sd)
        sdq <- quantile(apply(DAT4HM, 1, sd), 0.95)
        if(!sample2zero)  DAT4HM[sdr > sdq,] <- t(apply(DAT4HM[sdr > sdq,, drop=F], 1, function(x) (x/sd(x))*sdq))

##        mean(apply(DAT4HM[,wh2],1,mean)-apply(DAT4HM[,wh1],1,mean))

        DAT4HM <- cbind(modts, rep(0,nrow(DAT4HM)), DAT4HM)
        colnames(DAT4HM) <- c("mod t","",colnames(ys))
        DAT4HM[DAT4HM[, "mod t"] == 0, "mod t"] <- NA
        rownames(DAT4HM) <- rownames(ys)
#        DAT4HM <- DAT4HM
    }
    if(toplot){

        or <- order(obj$covar[, intvar])

        aux = heatmap.2(DAT4HM[order(DAT4HM[,1], decreasing=T),
            c(1:2, match(rownames(obj$covar), colnames(DAT4HM)))][,c(1,2, or+2)],
            dendrogram=dendrogram,col=col,
            trace=trace, margins = margins, notecol = notecol, notecex = notecex, keysize = keysize,
            cexCol = cexCol,  ColSideColors = c("gray","white", mycol[obj$covar[or, intvar]]), Rowv= Rowv,
            Colv = Colv, las=las, colsep=1, sepcolor="white", sepwidth=0.1, lwid=c(1, 2), na.color="lightgray", ...)

        legend("left", fill=c("gray", mycol), legend=c("mean-t", levels(obj$covar[or, intvar])))
    }
    if(pathwaylevel &fdrkey) DAT4HM <- DAT4HM2

    return(DAT4HM)
}


