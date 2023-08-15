enrichmentScore <- function(fchr,sign) {
    nfchr <- length(fchr)
    p <- numeric(nfchr)
    p[] <- -1 / (nfchr - length(sign))
    absfchr <- abs(fchr[sign])
    p[sign] <- absfchr / sum(absfchr)
    cumsum(p)
}

froastgsa <- function(y, index, design = NULL, contrast = ncol(design),
    set.statistic = "maxmean", psel = NULL, weights = NULL,
    nrot = 9999, shrink.resid = TRUE, restand = TRUE,
    mccores = 1, executation.info = TRUE, normalizeScores = TRUE, ...)
{

    dots <- names(list(...))
    y <- as.matrix(y)
    ngenes <- nrow(y)
    n <- ncol(y)

    ## index
    if (!is.list(index))
        index <- list(set = index)
    nsets <- length(index)
    if (nsets == 0)
        stop("index is empty")
    SetSizes <- unlist(lapply(index, length))

    if(mode(index[[1]]) =="character")
        index <- lapply(index, function(o) which(rownames(y) %in% o))

    ## design
    if (is.null(design))
        stop("design matrix not specified")
    else{
        design <- as.matrix(design)
        if (mode(design) != "numeric")
            stop("design must be a numeric matrix")
    }
    if (nrow(design) != n)
    stop("row dimension of design matrix must match column dimension of data")
    ne <- nonEstimable(design)
    if (!is.null(ne))
        stop(paste0("Coefficients not estimable:", paste(ne, collapse = " "),
                    "\n"))
    p <- ncol(design)
    if (p < 2L)
        stop("design needs at least two columns")
    p0 <- p - 1L
    d <- n - p

    ## contrast
    if (length(contrast) == 1L){
        contrast <- round(contrast)
        if (contrast < p){
            i <- 1L:p
            design <- design[, c(i[-contrast], contrast)]
        }
    }else{
        design <- contrastAsCoef(design = design, contrast = contrast,
                        first = FALSE)$design
    }

    if(length(contrast) == 1L)
        contrast <- ifelse(seq_len(p) == contrast, 1, 0)

    ## moderated t
    qr <- qr(design)
    signc <- sign(qr$qr[p, p])
    effects <- qr.qty(qr, t(y))
    s2 <- colMeans(effects[-seq_len(p), , drop = FALSE]^2)
    sv <- squeezeVar(s2, df = d)
    d0 <- sv$df.prior
    s02 <- sv$var.prior
    sd.post <- sqrt(sv$var.post)
    Y <- effects[-seq_len(p0), , drop = FALSE]
    YY <- colSums(Y^2)
    B <- Y[1, ]
    modt <- signc * B/sd.post

    if(normalizeScores)
        modt <- limma::zscoreT(modt, df =  d0 + d, approx = FALSE)

    if (shrink.resid){
        p.value <- 2 * pt(abs(modt), df = d0 + d, lower.tail = FALSE)
        proportion <- 1 - propTrueNull(p.value)
        stdev.unscaled <- rep_len(1/abs(qr$qr[qr$rank, qr$rank]), ngenes)
        var.unscaled <- stdev.unscaled^2
        df.total <- rep_len(d, ngenes) + sv$df.prior
        stdev.coef.lim <- c(0.1, 4)
        var.prior.lim <- stdev.coef.lim^2/sv$var.prior
        var.prior <- tmixture.vector(modt, stdev.unscaled, df.total,
                        proportion, var.prior.lim)
        if (any(is.na(var.prior))){
            var.prior[is.na(var.prior)] <- 1/sv$var.prior
            warning("Estimation of var.prior failed - set to default value")
        }
        r <- (var.unscaled + var.prior)/var.unscaled
        if (sv$df.prior > 10^6){
            kernel <- modt^2 * (1 - 1/r)/2
        }else{
            kernel <- (1 + df.total)/2 *
                log((modt^2 + df.total)/(modt^2/r + df.total))
        }
        lods <- log(proportion/(1 - proportion)) - log(r)/2 + kernel
        ProbDE <- exp(lods)/(1 + exp(lods))
        ProbDE[is.na(ProbDE)] <- 1
        Y[1,] <- Y[1,] * sqrt(var.unscaled/(var.unscaled + var.prior * ProbDE))
    }

    if(!is.null(psel)){
        modt <- modt[psel]
        sd.post <- sd.post[psel]
        YY <- YY[psel]
        Y <- Y[,psel]
        s2 <- s2[psel]

        names(sd.post) <- colnames(Y) <- names(s2) <-
            names(modt) <- names(YY) <- names(psel)
        ngenes <- length(modt)
    }

    set.statistic <- match.arg(set.statistic, c("mean", "mean.rank", "absmean",
        "median", "maxmean", "ksmax", "ksmean"))
    if(executation.info)
        pb <- txtProgressBar(min = 0, max = nrot, style = 3)

    if (set.statistic == "maxmean"){
        if(is.null(weights)){
            estpos <- vapply(index, function(o, modt)
                sum((modt[o] + abs(modt[o]))/2)/length(o), numeric(1),modt)
            estneg <- vapply(index, function(o, modt)
                sum((abs(modt[o]) - modt[o])/2)/length(o), numeric(1),modt)
        }
        else{
            estpos <- vapply(seq_len(length(index)),
                function(o, modt) sum(weights[[o]]*(modt[index[[o]]] +
                    abs(modt[index[[o]]]))/2), numeric(1), modt)

            estneg <- vapply(seq_len(length(index)),
                function(o, modt) sum( weights[[o]]*(abs(modt[index[[o]]]) -
                    modt[index[[o]]])/2), numeric(1), modt)
            names(estpos)  <-names(estneg)  <- names(index)
        }

        nuapos <- mean(modt*(modt > 0))
        nuaneg <- mean(-modt*(modt < 0))
        sdapos <- sd(modt*(modt > 0))
        sdaneg <- sd(modt*(modt < 0))
        if (restand){
            estpos <- (estpos - nuapos)/sdapos
            estneg <- (estneg - nuaneg)/sdaneg
        }
        est <- pmax(estpos, estneg)
        est[estneg > estpos] = -1 * est[estneg > estpos]
        p.value <- matrix(0, nrow = nsets, ncol = 2)
        est0.all <- do.call(cbind,mclapply(seq_len(nrot),function(i){
            if(executation.info ) setTxtProgressBar(pb, i)
            R <- matrix(rnorm((d + 1)), 1, d + 1)
            R <- R/sqrt(rowSums(R^2))
            Br <- colSums(Y*as.numeric(R))#R %*% Y
            s2r <- (YY - Br^2)/d
            if (is.finite(d0)){
                sdr.post <- sqrt((d0 * s02 + d * s2r)/(d0 + d))
            }else {
                sdr.post <- sqrt(s02)
            }
            modtr <- (signc * Br/sdr.post)
            if(normalizeScores)
                modtr <- limma::zscoreT(modtr, df =  d0 + d, approx = FALSE)

            if(is.null(weights)){
                estpos0 <- vapply(index, function(o, modtr)
                    sum((modtr[o] + abs(modtr[o]))/2)/length(o),
                        numeric(1), modtr)
                estneg0 <- vapply(index, function(o, modtr)
                    sum((abs(modtr[o]) - modtr[o])/2)/length(o),
                        numeric(1), modtr)
            }
            else{
                estpos0 <- vapply(seq_len(length(index)),
                    function(o, modtr) sum(weights[[o]]*(modtr[index[[o]]] +
                        abs(modtr[index[[o]]]))/2), numeric(1), modtr)
                estneg0 <- vapply(seq_len(length(index)),
                    function(o, modtr) sum(weights[[o]]*(abs(modtr[index[[o]]])-
                        modtr[index[[o]]])/2), numeric(1), modtr)
            }

            nuapos0 <- mean(modtr*(modtr > 0))
            nuaneg0 <- mean(-modtr*(modtr < 0))
            sdapos0 <- sd(modtr*(modtr > 0))
            sdaneg0 <- sd(modtr*(modtr < 0))
            if (restand){
                estpos0 <- (estpos0 - nuapos0)/sdapos0
                estneg0 <- (estneg0 - nuaneg0)/sdaneg0
            }
            est0 <- pmax(estpos0, estneg0)
            est0[estneg0 > estpos0] <- -1 * est0[estneg0 > estpos0]
                                        # p.value <- est0 <= est
            est0
        }, mc.cores = mccores))
        p.value <- vapply(seq_len(dim(est0.all)[1]), function(o)
            sum(est[o]>=est0.all[o,]), numeric(1))
        pval <- as.matrix(2 * c(pmin(p.value, nrot-p.value)))
        nes <- est/ apply(est0.all,1,sd)
    }

    if (set.statistic == "mean"){
        if(is.null(weights))
            est <- vapply(index, function(o, modt) mean(modt[o]),
                    numeric(1), modt)
        else
            est <- vapply(seq_len(length(index)), function(o, modt)
                sum(modt[index[[o]]] * weights[[o]]), numeric(1), modt)

        nua <- mean(modt)
        sda <- sd(modt)
        if (restand){
            est <- (est - nua)/sda
        }
        p.value <- matrix(0, nrow = nsets, ncol = 2)
        est0.all <- do.call(cbind,mclapply(seq_len(nrot),function(i){
            if(executation.info ) setTxtProgressBar(pb, i)
            R <- matrix(rnorm((d + 1)), 1, d + 1)
            R <- R/sqrt(rowSums(R^2))
            Br <- colSums(Y*as.numeric(R))#R %*% Y
            which(is.na(Br))
            s2r <- (YY - Br^2)/d
            if (is.finite(d0)){
                sdr.post <- sqrt((d0 * s02 + d * s2r)/(d0 + d))
            }else{
                sdr.post <- sqrt(s02)
            }
            modtr <- signc * Br/sdr.post
            if(normalizeScores)
                modtr <- limma::zscoreT(modtr, df =  d0 + d, approx = FALSE)
            if(is.null(weights))
                est0 <- vapply(index, function(o, modtr) mean(modtr[o]),
                        numeric(1), modtr)
            else
                est0 <- vapply(seq_len(length(index)), function(o, modtr)
                    sum(modtr[index[[o]]] * weights[[o]]),
                        numeric(1), modtr)

            nua0 <- mean(modtr)
            sda0 <- sd(modtr)
            if (restand){
                est0 <- (est0 - nua0)/sda0
            }
            est0 #p.value <- est0 <= est #
        }, mc.cores = mccores))
        p.value <- vapply(seq_len(dim(est0.all)[1]), function(o)
            sum(est[o]>=est0.all[o,]), numeric(1))
        pval <- as.matrix(2 * c(pmin(p.value, nrot-p.value)))
        nes <- est/ apply(est0.all,1,sd)
    }

    if (set.statistic == "mean.rank"){
        obs.ranks <- matrix(0, ngenes, 3)
        obs.ranks[, 1] <- rank(modt)
        est <- vapply(index, function(o, modt) mean(modt[o]), numeric(1), modt)
        nua <- mean(modt)
        sda <- sd(modt)
        if (restand){
            est <- (est - nua)/sda
        }
        obs.ranks[, 2] <- ngenes - obs.ranks[, 1] + 1
        obs.ranks[, 3] <- rank(abs(modt))
        AllIndices <- unlist(index)
        Set <- rep(seq_len(nsets), SetSizes)
        obs.set.ranks <- rowsum(obs.ranks[AllIndices, ], group = Set,
                                reorder = FALSE)
        obs.set.ranks <- obs.set.ranks/SetSizes
        rot.ranks <- obs.ranks
        est <- mean(rank(modt)) - obs.set.ranks[,2]
        p.value <- matrix(0, nrow = nsets, ncol = 3)
        p.value <- mclapply(seq_len(nrot),function(i){
            if(executation.info ) setTxtProgressBar(pb, i)
            R <- matrix(rnorm((d + 1)), 1, d + 1)
            R <- R/sqrt(rowSums(R^2))
            Br <- colSums(Y*as.numeric(R))#R %*% Y
            s2r <- (YY - Br^2)/d
            if (is.finite(d0))
                sdr.post <- sqrt((d0 * s02 + d * s2r)/(d0 + d))
            else sdr.post <- sqrt(s02)
            modtr <- signc * Br/sdr.post
            if(normalizeScores)
                modtr <- limma::zscoreT(modtr, df =  d0 + d, approx = FALSE)
            rot.ranks[, 1] <- rank(modtr)
            rot.ranks[, 2] <- ngenes - rot.ranks[, 1] + 1
            rot.ranks[, 3] <- rank(abs(modtr))
            rot.set.ranks <- rowsum(rot.ranks[AllIndices, ],
                                    group = Set, reorder = FALSE)
            rot.set.ranks <- rot.set.ranks/SetSizes
            (rot.set.ranks <= obs.set.ranks)
        }, mc.cores = mccores)
        pval <- vapply(c(1,3), function(o)
            apply(vapply(p.value,function(x) x[,o],
                logical(length(index))),1,sum),numeric(length(index)))
        pval[,2] <- nrot - pval[,2]
        pval[,1] <- 2 * pmin(pval[,1],nrot - pval[,1])
        colnames(pval) <- c("pval.diff","pval.abs")
        rownames(pval) <- names(index)
        modt <- rank(modt) - mean(rank(modt))
    }

    if (set.statistic == "median"){
        est <- vapply(index, function(o, modt) median(modt[o]),
            numeric(1), modt)
        nua <- median(modt)
        sda <- mad(modt)
        if (restand)
            {
                est <- (est - nua)/sda
            }
        p.value <- matrix(0, nrow = nsets, ncol = 2)
        p.value <- apply(do.call(cbind,mclapply(seq_len(nrot),function(i){
            if(executation.info ) setTxtProgressBar(pb, i)
            R <- matrix(rnorm((d + 1)), 1, d + 1)
            R <- R/sqrt(rowSums(R^2))
            Br <- colSums(Y*as.numeric(R))#R %*% Y
            s2r <- (YY - Br^2)/d
            if (is.finite(d0)){
                sdr.post <- sqrt((d0 * s02 + d * s2r)/(d0 + d))
            }else{
                sdr.post <- sqrt(s02)
            }
            modtr <- signc * Br/sdr.post
            if(normalizeScores)
                modtr <- limma::zscoreT(modtr, df =  d0 + d, approx = FALSE)
            est0 <- vapply(index, function(o, modtr) median(modtr[o]),
                    numeric(1), modtr)
            nua0 <- median(modtr)
            sda0 <- mad(modtr)
            if (restand){
                est0 <- (est0 - nua0)/sda0
            }
            p.value <- est0 <= est
        }, mc.cores = mccores)),1,sum)
        pval <- as.matrix(2 * c(pmin(p.value, nrot-p.value)))
    }

    if (set.statistic == "absmean"){
        nua <- mean(modt)
        sda <- sd(modt)

        if (restand) modt1 <- (modt - nua)/sda
        else modt1 <- modt
        if(is.null(weights))
            est <- vapply(index, function(o, modt1) mean(abs(modt1[o])),
                    numeric(1), modt1)
        else
            est <- vapply(seq_len(length(index)), function(o, modt1)
                sum(modt1[index[[o]]] * weights[[o]]), numeric(1), modt1)

        p.value <- matrix(0, nrow = nsets, ncol = 2)
        p.value <- apply(do.call(cbind,mclapply(seq_len(nrot),function(i){
            if(executation.info ) setTxtProgressBar(pb, i)
            R <- matrix(rnorm((d + 1)), 1, d + 1)
            R <- R/sqrt(rowSums(R^2))
            Br <- colSums(Y*as.numeric(R))#R %*% Y
            s2r <- (YY - Br^2)/d
            if (is.finite(d0))
                {
                    sdr.post <- sqrt((d0 * s02 + d * s2r)/(d0 + d))
                }else
                    {
                        sdr.post <- sqrt(s02)
                    }
            modtr <- signc * Br/sdr.post
            if(normalizeScores)
                modtr <- limma::zscoreT(modtr, df =  d0 + d, approx = FALSE)

            nua0 <- mean(modtr)
            sda0 <- sd(modtr)

            if (restand) modtr1 <- (modtr - nua0)/sda0
            else modtr1 <- modtr
            if(is.null(weights))
                est0 <- vapply(index, function(o, modtr1)
                    mean(abs(modtr1[o])), numeric(1), modtr1)
            else
                est0 <- vapply(seq_len(length(index)), function(o, modtr1)
                    sum(modtr1[index[[o]]] * weights[[o]]), numeric(1), modtr1)

            p.value <- est0 >= est
        }, mc.cores = mccores)),1,sum)
        pval <- as.matrix(p.value)
    }

    if (set.statistic == "ksmax"){
        if (restand){
            nua <- mean(modt)
            sda <- sd(modt)
            modt <- (modt - nua)/sda
        }

        ord <- order(-modt)
        modt2 <- modt[ord]
        rankmod <- rank(-modt, ties.method = "first")
        index2 <- lapply(index, function(o) rankmod[o])

        es <- vapply(index2, function(o, modt2)
            enrichmentScore(modt2, o),numeric(length(modt2)),modt2)
        est <- apply(es, 2, function(o){
            es.range <- range(o)
            es.range[abs(es.range)==max(abs(es.range))][1]
        })
        p.value <- matrix(0, nrow = nsets, ncol = 2)
        est0.all <- do.call(cbind,mclapply(seq_len(nrot),function(i){
            if(executation.info ) setTxtProgressBar(pb, i)
            R <- matrix(rnorm((d + 1)), 1, d + 1)
            R <- R/sqrt(rowSums(R^2))
            Br <- colSums(Y*as.numeric(R))#R %*% Y
            (s2r <- (YY - Br^2)/d)

            if (is.finite(d0)) {
                sdr.post <- sqrt((d0 * s02 + d * s2r)/(d0 + d))
            }else{
                sdr.post <- sqrt(s02)
            }
            modtr <- signc * Br/sdr.post
            if(normalizeScores)
                modtr <- limma::zscoreT(modtr, df =  d0 + d, approx = FALSE)
            if (restand){
                nua <- mean(modt)
                sda <- sd(modt)
                modt <- (modt - nua)/sda
            }
            ordr <- order(-modtr)
            modtr2 <- modtr[ordr]
            rankmod <- rank(-modtr, ties.method = "first")
            index2 <- lapply(index, function(o) rankmod[o])
            es0 <- vapply(index2, function(o, modtr2)
                enrichmentScore(modtr2, o),numeric(length(modtr2)),modtr2)
            est0 <- apply(es0, 2, function(o){
                es.range <- range(o)
                es.range[abs(es.range)==max(abs(es.range))][1]
            })
            est0
        }, mc.cores = mccores))
        p.value <- vapply(seq_len(dim(est0.all)[1]), function(o)
            sum(est[o] >= as.numeric(est0.all[o,])), numeric(1))
        pval <- as.matrix(2 * c(pmin(p.value, nrot-p.value)))
        nes <- vapply(seq_len(dim(est0.all)[1]), function(o){
            est0.a <- est0.all[o,]
            est[o] / abs(mean(est0.a[sign(est0.a) == sign(est[o])]))
        },numeric(1))
    }
    if (set.statistic == "ksmean"){
        ord <- order(-modt)
        modt2 <- modt[ord]
        rankmod = rank(-modt, ties.method = "first")
        index2 <- lapply(index, function(o) rankmod[o])
        es <- vapply(index2, function(o, modt2)
            enrichmentScore(seq_len(length(modt2)) - length(modt2)/2, o),
                numeric(length(modt2)),modt2)
        estpos <- apply(es, 2, function(o)
            if (any(o > 0)) o[o==max(o)][1] else 0)
        estneg <- apply(es, 2, function(o)
            if (any(o < 0)) o[o==min(o)][1] else 0)
        est <- abs(estpos) - abs(estneg)
        est0.all <- do.call(cbind,mclapply(seq_len(nrot),function(i){
            if(executation.info ) setTxtProgressBar(pb, i)
            R <- matrix(rnorm((d + 1)), 1, d + 1)
            R <- R/sqrt(rowSums(R^2))
            Br <- colSums(Y*as.numeric(R))#R %*% Y
            s2r <- (YY - Br^2)/d
            if (is.finite(d0)){
                sdr.post <- sqrt((d0 * s02 + d * s2r)/(d0 + d))
            }else{
                sdr.post <- sqrt(s02)
            }
            modtr <- signc * Br/sdr.post
            if(normalizeScores)
                modtr <- limma::zscoreT(modtr, df =  d0 + d, approx = FALSE)
            ordr <- order(-modtr)
            modtr2 <- modtr[ordr]
            rankmod <- rank(-modtr, ties.method = "first")
            index2 <- lapply(index, function(o) rankmod[o])
            es0 <- vapply(index2, function(o, modtr2)
                enrichmentScore(seq_len(length(modtr2)) - length(modtr2)/2, o),
                    numeric(length(modtr2)),modtr2)
            estpos0 <- apply(es0, 2, function(o)
                if (any(o > 0)) o[o==max(o)][1] else 0)
            estneg0 <- apply(es0, 2, function(o)
                if (any(o < 0)) o[o==min(o)][1] else 0)
            est0 <- abs(estpos0) - abs(estneg0)
            est0
        }, mc.cores = mccores))
        p.value <- vapply(seq_len(dim(est0.all)[1]), function(o)
            sum(est[o] >= as.numeric(est0.all[o,])),numeric(1))
        pval <- as.matrix(2 * c(pmin(p.value, nrot-p.value)))
    }
    r <- (pval + 1)/(nrot + 1)
    rownames(r) <- names(index)
    if (set.statistic %in% c('ksmax',"mean", "maxmean")){
        r <- cbind(NGenes = SetSizes, est=est, nes = nes, r)
        r <- data.frame(r)
        colnames(r) <- c("ngenes", "est", "nes", "pval")
        r$adj.pval <- p.adjust(r$pval, method = "BH")
        r$sign <- sign(r$est)
        r <- r[order(-abs(r$nes)),c("ngenes", "est", "nes", "pval", "adj.pval")]
    }
    else{
        r <- cbind(NGenes = SetSizes, est=est, r)
        r <- data.frame(r)
        if (set.statistic!='mean.rank'){
            colnames(r) <- c("ngenes", "est", "pval")
            r$adj.pval <- p.adjust(r$pval, method = "BH")
            r$sign <- sign(r$est)
            r <- r[order(r$pval), c("ngenes", "est", "pval", "adj.pval")]
        }else{
            colnames(r)[1] <- "ngenes"
            colnames(r)[2] <- "est"
            colnames(r)[3] <- "pval.diff"
            colnames(r)[4] <- "pval.mixed"
            r$adj.pval.diff <- p.adjust(r$pval.diff, method = "BH")
            r$adj.pval.mixed <- p.adjust(r$pval.mixed, method = "BH")
            r$sign <- sign(r$est)
            r <- r[order(r$pval.diff),
                c("ngenes", "est", "pval.diff", "adj.pval.diff","pval.mixed",
                    "adj.pval.mixed")]
        }
    }
    return(list(res = r, stats = modt, contrast = contrast))
}



