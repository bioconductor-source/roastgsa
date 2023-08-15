
########################################################
#########  function:  varrotrand     ###################
########################################################

varrotrand <- function(obj, y, testedsizes = c(3:30,seq(32,50, by=2),
    seq(55,200,by=5)), nrep = 200, nrot = NULL, mccores = NULL, psel = NULL){

    if(is.null(nrot)) nrot <- attr(obj,"nrot")

    if(is.null(mccores)) mccores <- attr(obj,"mccores")

    psel2  <- psel
    weigths <- NULL
    k <- 1
    indexauxall <- lapply(testedsizes, function(k){
        if(is.null(psel2))
            indexaux <- lapply(seq_len(nrep), function(o)
                sample(rownames(y), k))
        else
            indexaux <- lapply(seq_len(nrep), function(o)
                sample(names(psel2), k))
        indexaux
    })
    indexauxall <- do.call(c,indexauxall)

    obj2 <- roastgsa(y, form = obj$form, covar = obj$covar,
                contrast = obj$contrast, index = indexauxall,
                nrot = nrot, mccores = mccores, weights = weigths,
                set.statistic = obj$statistic,
                self.contained = obj$self.contained, psel = psel2)

    varrot <- vapply(testedsizes, function(k){
        wh <- which(obj2$res$measured_genes ==k)
        obj2$res[["est"]][wh]/obj2$res[["nes"]][wh]
    }, numeric(nrep))
    varrot <- list(varrot = varrot, testedsizes = testedsizes, nrep=nrep)
    class(varrot) <- "varrotrand"
    return(varrot)
}
