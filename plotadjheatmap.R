
########################################################
#########  function:  plot adjheatmapMean  #############
########################################################

heatmapMean <- function(obj, y, intvar, whplot = 1, toplot = TRUE, pathwaylevel = FALSE,
                           mycol = c("black","orange","green","white"),
                           dendrogram =  "n", col=  bluered(100),trace='none',
                           margins=c(8,16),notecol='black',notecex=1, keysize=.9, cexCol=1.5,
                           Rowv = TRUE, Colv = FALSE, las =2, fdrkey = FALSE, ...){

    if(is.numeric(whplot)) whplot <- rownames(obj$res)[whplot]

    if(missing(intvar)) intvar <- tail(colnames(attr(terms(as.formula(obj$form)), "factors")),1)
    whint <- which(colnames(obj$covar) == intvar)
    wh2 <- obj$covar[,whint] == levels(obj$covar[,whint])[2]
    wh1 <- obj$covar[,whint] == levels(obj$covar[,whint])[1]

    if(pathwaylevel){
        if(length(whplot)<2) stop("whplot must be of lenght > 1")
        DAT4HM <- t(sapply(whplot, function(j){
            hm <- heatmapMean(obj, y, whplot = j, toplot = FALSE)
            apply(hm,2,mean)
        }))
        DAT4HM <- DAT4HM/(obj$res[whplot,"est"]/ obj$res[whplot,"nes"])
        colnames(DAT4HM)[1] <- "nes"
        if(fdrkey){
            DAT4HM2 <- DAT4HM
            rownames(DAT4HM) <- paste0(rownames(DAT4HM), formatApval(obj$res[whplot,]$adj.pval))
        }
    }
    else{
        lmf <- lmFit(y, obj$design)
        coefs <- lmf$coef
        coefs[is.na(coefs)] <- 0

        if(ncol(obj$design)>2)
            yr <- y - coefs[,-obj$contrast]%*%t(obj$design[,-obj$contrast])
        else
            yr <- y

        ys <- yr[obj$index[[whplot]], ]

        modts <- obj$stats[obj$index[[whplot]]]
        sde <- (coefs[names(modts),obj$contrast]/modts)

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
                    DAT4HM[obj$stats[obj$index[[whplot]]]< 0, ] <- 0
                    DAT4HM[,wh1] <-  (DAT4HM[,wh1] + nuapos/2)/sdapos
                    DAT4HM[,wh2] <-  (DAT4HM[,wh2] - nuapos/2)/sdapos
                }
                else{
                    DAT4HM[obj$stats[obj$index[[whplot]]]> 0, ] <- 0
                    DAT4HM[,wh1] <-  (DAT4HM[,wh1] - nuaneg/2)/sdaneg
                    DAT4HM[,wh2] <-  (DAT4HM[,wh2] + nuaneg/2)/sdaneg
                }
            }
            else{
                DAT4HM[,wh1] <-  (DAT4HM[,wh1] + mean(obj$stats)/2)/sd(obj$stats)
                DAT4HM[,wh2] <-  (DAT4HM[,wh2] - mean(obj$stats)/2)/sd(obj$stats)
            }
        }
        DAT4HM <- t(apply(DAT4HM, 1, function(x) (x - mean(x))))
        mean(apply(DAT4HM[,wh2],1,mean)-apply(DAT4HM[,wh1],1,mean))

        DAT4HM <- cbind(modts, rep(0,nrow(DAT4HM)), DAT4HM)
        colnames(DAT4HM) <- c("mod t","",colnames(ys))
        rownames(DAT4HM) <- rownames(ys)

        DAT4HM <- DAT4HM[,c(1:2,which(wh1) +2,which(wh2)+2)]
    }
    if(toplot){
      aux = heatmap.2(DAT4HM,dendrogram=dendrogram,col=col, trace=trace, margins = margins,
      notecol = notecol, notecex = notecex, keysize = keysize, cexCol = cexCol,
        ColSideColors = mycol[c(3, 4, rep(1,sum(wh1)), rep(2,sum(wh2)))], Rowv= Rowv, Colv = Colv, las=las, ...)
    }
    if(pathwaylevel &fdrkey) DAT4HM <- DAT4HM2

    return(DAT4HM)
}


