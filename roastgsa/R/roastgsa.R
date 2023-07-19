
##############################################
#########  function:  roastgsa    ############
##############################################

roastgsa <- function(y, covar, form, contrast = NA, design = NULL, gsetsel,
    gspath,  index = NULL, self.contained = FALSE, set.statistic = "maxmean",
    psel = NULL, nrot = 9999, minsize = 10,
    maxsize = 500,  mccores = 1, executation.info = TRUE,
    weights = NULL, shrink.resid = TRUE, normalizeScores = TRUE, ...){

    fcall <-   match.call()

    # Matrix design
    if(is.null(design))
        design <- model.matrix(as.formula(form), data = covar)
    terms <- colnames(design)
    if(is.na(contrast[1])) contrast = ncol(design)
    if(is.character(contrast)) contrast <- which(colnames(design) == contrast)

    # probset level?
    if(is.null(psel)) y2 <- y
    else{
        y2 <- y[psel,]
        rownames(y2) <- names(psel)
    }

    # Index
    if(is.null(index)){
        index.o <- indexFunction(y2, gsetsel, gspath = gspath,
            minsize = minsize,maxsize = maxsize, geneRandom = FALSE,
            forPRgsea = TRUE)

        index <- indexFunction(y2, gsetsel, gspath = gspath, minsize = minsize,
            maxsize = maxsize, geneRandom = FALSE, forPRgsea = FALSE)

        index.o <- index.o[names(index)]
        index2 <- lapply(index.o, function(x) x[x%in%rownames(y2)])
    }
    else{
        if(mode(index[[1]]) == "numeric")
            index.o <- lapply(index, function(x) rownames(y2)[x])
        else
            index.o <- index
        index2 <- lapply(index.o, function(x) x[x%in%rownames(y2)])
        index <- lapply(index.o, function(x) which(rownames(y2)%in%x))
    }

    # froastgsa function
    res <- froastgsa(y = y, index = index, design = design, contrast = contrast,
            set.statistic = set.statistic, nrot = nrot, mccores = mccores,
            psel = psel, executation.info = executation.info,
            restand = !self.contained, weights = weights,
            normalizeScores = normalizeScores, ...)

    # results
    res$res <- data.frame(vapply(index.o,length,
        numeric(1))[rownames(res$res)],res$res)
    names(res$res)[c(1,2)] <- c("total_genes","measured_genes")
    res$index <- index2
    res$covar <- covar
    res$statistic <- set.statistic
    res$design <- design
    res$contrast <- res$contrast
    res$self.contained <- self.contained
    res$form <- form
    res$fcall <- fcall
    class(res) <- c("roastgsa")

    attr(res, "nrot") <- nrot
    attr(res, "mccores") <- mccores

    return(res)

}







