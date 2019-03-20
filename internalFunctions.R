
## average correlation
avgcorfun <- function(x){
    x <- t(scale(t(x)))
    p <- dim(x)[1]
    n <- dim(x)[2]
    c1 <- ((rep(1,p) %*% x) %*% (t(x) %*% as.matrix(rep(1,p)))) / ((n-1)*p^2)
    return( (c1 - 1/p) * p/(p-1))
}


formatApval <- function(x, ns = "", nc = ""){
    st = rep(ns,length(x))
    st <- ifelse(x < 0.25, "+", st)
    st <- ifelse(x < 0.10, "*", st)
    st <- ifelse(x < 0.05, "**", st)
    st <- ifelse(x < 0.01, "***", st)
    if(any(is.na(st))) st[is.na(st)] <- nc
    return(st)
}


## design from formula and expressionSet
designFunction <- function(nd, form = NULL){
    if(is.null(form)) form <-  paste0("~ -1 + ev", sep='')
    design <- model.matrix(as.formula(form), data=nd)
    return(design)
}

## contrast matrix from formula and groups
contMatFunction <- function(pd, form = NULL, lgroups = NULL){
    lsigns <- c(sapply(1:length(lgroups), function(j) c(1, -1), simplify = F))
    conts <- sapply(1:length(lsigns), function(j)
                make.contrasts(d=pd, form=as.formula(form),
                               lsel=lgroups[[j]],
                               signs=lsigns[[j]]), simplify=F)
    cont.mat <- sapply(conts, function(o) o$cont)
    colnames(cont.mat) <- names(lgroups)
    return(cont.mat)
}

## differential expression using Limma
dsfunction <- function(nd, form = NULL, lgroups = NULL){
 design <- designFunction(nd, form = form)
 lmf <- lmFit(nd, design)
 betas <- coef(lmf)
 cont.mat <- contMatFunction(pData(nd), form = form, lgroups = lgroups)

 lmc <- contrasts.fit(lmf, cont.mat)
 eb <- eBayes(lmc)
 ds <- topTable(eb, adjust="BH", number=nrow(exprs(nd)), coef=colnames(cont.mat)[1], confint=T)

 return(ds)
}


