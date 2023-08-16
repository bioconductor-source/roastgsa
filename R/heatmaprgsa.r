########################################################
#########  function:  plot adjheatmapMean  #############
########################################################

heatmaprgsa <- function(obj, y, intvar, adj.var = NULL, whplot = 1,
    toplot = TRUE, pathwaylevel = FALSE,
    mycol = c("black","orange","green","white"), sample2zero = FALSE,
    rgsa.like=FALSE, psel = NULL, dendrogram =  "n", col=  bluered(100),
    trace='none', margins=c(8,16),notecol='black', notecex=1, keysize=.9,
    cexCol=1.5, Rowv = NULL, Colv = FALSE, las =2, fdrkey = FALSE,
    quantile.sat = 0.95, order1= NULL, order2= NULL,...){

    if(!inherits(obj, "roastgsa")) stop("obj must be of class roastgsa")

    if(is.numeric(whplot)) whplot <- rownames(obj$res)[whplot]

    psel2  <- psel
    if(!is.null(psel2))
    {
        y <- y[psel2,]
        rownames(y) <- names(psel2)
        psel2 <- NULL
    }

    if(missing(intvar))
        intvar <- tail(colnames(attr(terms(as.formula(obj$form)), "factors")),1)

    if(any(obj$contrast < 0)){
        if(any(obj$contrast > 0)){
            wh1 <- rowSums(as.matrix(obj$design[,obj$contrast > 0])) > 0
            wh2 <- rowSums(as.matrix(obj$design[,obj$contrast < 0])) > 0
        }
        else{
            whaux  <- which(regexpr(intvar,names(obj$contrast))>0)
            wh1 <- rowSums(as.matrix(obj$design[,whaux])) == 0
            wh2 <- rowSums(as.matrix(obj$design[,obj$contrast < 0])) > 0
        }
    }
    else{
        wh1 <- rowSums(as.matrix(obj$design[,obj$contrast > 0])) > 0
        wh2 <- rowSums(as.matrix(obj$design[,obj$contrast > 0])) <= 0
    }
    sample2zeroaux  <- sample2zero

    if(pathwaylevel){
        if(length(whplot)<2) stop("whplot must be of lenght > 1")

        DAT4HM <- do.call(rbind,lapply(whplot, function(j){
            hm <- heatmaprgsa(obj, y, whplot = j, toplot = FALSE,
                    sample2zero  = TRUE, adj.var = adj.var,
                    psel = psel2, rgsa.like = obj$statistic != "mean.rank", ...)
            apply(hm,2,mean)
        }))

        if(obj$statistic != "mean.rank"){
            DAT4HM <- DAT4HM/(obj$res[whplot,"est"]/ obj$res[whplot,"nes"])
            colnames(DAT4HM)[1] <- "nes"
        }

        sdr <- apply(DAT4HM[,-seq_len(2)], 1, sd)
        sdq <- quantile(apply(DAT4HM[,-seq_len(2)], 1, sd), quantile.sat)
        DAT4HM[sdr > sdq,-seq_len(2)] <-
            t(apply(DAT4HM[sdr > sdq,-seq_len(2), drop=FALSE],1,
                function(x) (x/sd(x))*sdq))

        if(obj$statistic != "mean.rank"){
            if ("nes"%in%colnames(obj$res))
            {
                DAT4HM[,1] <- obj$res[whplot,"nes"]
            }else if ("est"%in%colnames(obj$res))
            {
                DAT4HM[,1] <- obj$res[whplot,"est"]
            }
        }
        if(fdrkey){
            DAT4HM2 <- DAT4HM
            rownames(DAT4HM) <- paste0(rownames(DAT4HM),
                        formatApval(obj$res[whplot,]$adj.pval))
            if(obj$statistic == "mean.rank")
                rownames(DAT4HM) <- paste0(rownames(DAT4HM),
                        formatApval(obj$res[whplot,]$adj.pval.diff))
        }

    }
    else{
        lmf <- lmFit(y, obj$design)
        coefs <- lmf$coef
        coefs[is.na(coefs)] <- 0

        if(!is.null(adj.var))
            yr <- y - coefs[, adj.var] %*% t(obj$design[,adj.var])
        else
            yr <- y

        ys <- yr[obj$index[[whplot]], ]


        modts <- obj$stats[obj$index[[whplot]]]
        if(obj$statistic == "mean.rank")
            sde <- rep(1,length(modts))
        else{
            sde <- (coefs%*%as.matrix(obj$contrast))[names(modts), 1]
            sde <- sde/modts
        }

        if(obj$statistic == "maxmean")
            modts <- if(obj$res[whplot,"est"]>0) modts*(modts>0) else
                modts*(modts<0)

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

        if(any(sde==0)|any(is.nan(sde))) sde[which(sde==0 | is.nan(sde))] <- Inf
            DAT4HM <- ys/sde

        if(rgsa.like){
            if(obj$statistic == "maxmean" & !obj$self.contained){
                if(obj$res[whplot,"est"]>0){
                    if(sample2zero){
                        DAT4HM[obj$stats[obj$index[[whplot]]]< 0, ] <- 0
                        DAT4HM[,wh1] <-  (DAT4HM[,wh1] - nuapos/2)/sdapos
                        DAT4HM[,wh2] <-  (DAT4HM[,wh2] + nuapos/2)/sdapos
                        DAT4HM[,!(wh2 | wh1)] <-  (DAT4HM[,!(wh1 | wh2)])/sdapos
                    }
                    else{
                        DAT4HM <- DAT4HM/sdapos
                        modts[obj$stats[obj$index[[whplot]]]< 0] <- 0
                    }
                }
                else{
                    if(sample2zero){
                        DAT4HM[obj$stats[obj$index[[whplot]]]> 0, ] <- 0
                        DAT4HM[,wh1] <-  (DAT4HM[,wh1] + nuaneg/2)/sdaneg
                        DAT4HM[,wh2] <-  (DAT4HM[,wh2] - nuaneg/2)/sdaneg
                        DAT4HM[,!(wh2 | wh1)] <-  (DAT4HM[,!(wh1 | wh2)])/sdaneg
                    }
                    else{
                        modts[obj$stats[obj$index[[whplot]]] > 0] <- 0
                        DAT4HM <- DAT4HM/sdaneg
                    }
                }
            }
            else{
                if(sample2zero){
                    DAT4HM[,wh1] <- (DAT4HM[,wh1] -
                        mean(obj$stats)/2)/sd(obj$stats)
                    DAT4HM[,wh2] <- (DAT4HM[,wh2] +
                        mean(obj$stats)/2)/sd(obj$stats)
                }
                else
                    DAT4HM <- DAT4HM/sd(obj$stats)
            }
        }
        DAT4HM <- t(apply(DAT4HM, 1, function(x) (x - mean(x))))

        sdr <- apply(DAT4HM, 1, sd)
        sdq <- quantile(apply(DAT4HM, 1, sd), quantile.sat)
        if(!sample2zero)
            DAT4HM[sdr > sdq,] <- t(apply(DAT4HM[sdr > sdq,, drop=FALSE],
                                    1, function(x) (x/sd(x))*sdq))

        DAT4HM <- cbind(modts, rep(0,nrow(DAT4HM)), DAT4HM)
        colnames(DAT4HM) <- c("mod t","",colnames(ys))
        DAT4HM[DAT4HM[, "mod t"] == 0, "mod t"] <- NA
        rownames(DAT4HM) <- rownames(ys)
    }
    if(toplot){
        if(is.null(order2)) or <- order(obj$covar[, intvar])
        else or <- order2
        if(is.null(order1)) or2  <- order(DAT4HM[,1], decreasing=TRUE)
        else or2  <- order1
        aux <- heatmap.2(DAT4HM[or2,
            c(seq_len(2), match(rownames(obj$covar),
            colnames(DAT4HM)))][,c(1,2, or+2)], dendrogram=dendrogram,col=col,
            trace=trace, margins = margins,  notecol = notecol,
            notecex = notecex, keysize = keysize, cexCol = cexCol,
            ColSideColors = c("gray","white", mycol[obj$covar[or, intvar]]),
            Rowv= Rowv, Colv = Colv, las=las, colsep=1, sepcolor="white",
            sepwidth=0.1, lwid=c(1,2), na.color="lightgray", ...)

        legend("left", fill=c("gray", mycol),
            legend=c("mean-t", levels(obj$covar[or, intvar])))
    }
    if(pathwaylevel &fdrkey) DAT4HM <- DAT4HM2

    return(DAT4HM)
}


