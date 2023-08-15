ssGSA <- function(y, obj=NULL, design = NULL, contrast = NULL, index = NULL,
                method = c("GScor","GSadj","zscore")){

    if(!is.null(obj)){
        if(!inherits(obj, "roastgsa")) stop("obj must be of class roastgsa")

        design  <- obj$design
        contrast  <- obj$contrast
        index  <- obj$index
    }
    if(length(contrast)==1 & ncol(design)>1){
        contrast2  <- rep(0,ncol(design))
        if(is.numeric(contrast))
            contrast2[contrast]  <- 1
        if(is.character(contrast))
            contrast2[which(names(design)==contrast)]  <- 1
        contrast  <- contrast2
    }

    y <- t(scale(t(y)))
    GS <- as.numeric( (rep(1,dim(y)[1]) %*% y)/dim(y)[1])
    gsvaD <- t(vapply(index, function(o)
        (as.numeric(rep(1,length(o))%*%y[o,])/length(o)),numeric(ncol(y))))
    colnames(gsvaD) <- colnames(y)

    if(method[1] == "GSadj"){
        design <- cbind(GS,design)
        contrast <- c(0, contrast)
    }
    if(method[1] == "GScor"){
        gsvaD  <- gsvaD  - GS
    }

    dsgs <- t(apply(gsvaD, 1, function(x){
        lm1 <- summary(lm(x ~ design))$coeff
        t(t(lm1)%*%contrast)
    }))

    if(method == "GSadj"){
        lmf2 <- lmFit(gsvaD, design)
        gsvaD <- gsvaD - lmf2$coef[,"GS"]%*%t(design[,"GS"])
    }
    colnames(dsgs) <- c("est","sd","t","pval")
    dsgs <- as.data.frame(dsgs)
    dsgs <- dsgs[order(-abs(dsgs$t)),]
    dsgs$adj.pval <- p.adjust(dsgs$pval,method="BH")

    dsgs$ngenes <- vapply(index[rownames(dsgs)],length, numeric(1))
    dsgs <- dsgs[,c("ngenes","est","pval","adj.pval","t")]

    x  <- list(res = dsgs, stats = gsvaD)
    class(x)  <- "ssGSA"
    return(x)
}


