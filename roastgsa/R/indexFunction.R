## Index gene sets in the expression matrix used
indexFunction <- function(y, gspath, gsetsel ="Hallmarks",
    minsize=10, maxsize=300, geneRandom = FALSE, forPRgsea = FALSE){
    gsets <- paste0(gspath, list.files(gspath))
    gsets <- gsets[regexpr("\\.gmt", gsets) > 0]
    gsets <- gsets[regexpr(gsetsel,gsets)>0]

    auxg <- strsplit(readLines(gsets),"\t")
    GSET <- lapply(auxg,function(x) x[3:length(x)])
    names(GSET) <- vapply(auxg,function(x) x[1], character(1))
    if(forPRgsea)  index <- GSET
    else{
        index <- lapply(GSET, function(x) which(rownames(y)%in%x))
        if(geneRandom)
            index <- lapply(index, function(x)
                sample(seq_len(dim(y)[1]), length(x)))
        len <- vapply(index,length, numeric(1))
        index <- index[len > minsize & len < maxsize]
    }
    return(index)
}

