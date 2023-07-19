## batchcorlm
batchcorlm <- function(y, design, adj.var){
    lmf1 <- lmFit(y, design[,adj.var])
    y <- y - lmf1$coef[,adj.var]%*%t(design[,adj.var])
    return(y)
}

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
    lsigns <- c(lapply(seq_len(length(lgroups)), function(j) c(1, -1)))

    conts <- lapply(seq_len(length(lsigns)), function(j)
        make.contrasts(d=pd, form=as.formula(form),
            lsel=lgroups[[j]],
            signs=lsigns[[j]]))
    cont.mat <- vapply(conts, function(o)
        o$cont, numeric(length(conts[[1]]$cont)))
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
    ds <- topTable(eb, adjust.method = "BH", number=nrow(exprs(nd)),
            coef=colnames(cont.mat)[1], confint=TRUE)
    return(ds)
}

make.contrasts <- function(d, form, lsel, signs, type='balanced', ngroups=NULL)
{

    nv <- strsplit(paste(form)[2], split=' \\+ ')[[1]];
    nv <- nv[regexpr("\\:", nv) < 0];
    nv <- nv[nv!='-1']
    for (n in nv) if (is.character(d[, n])) d[, n] <- factor(d[, n])

    nvf <- nv[vapply(nv, function(v, d) is.factor(d[, v]),logical(1),d)];
    nvn <- nv[!nv%in%nvf];

    daf <- expand.grid(lapply(nvf, function(o) levels(d[, o])));
    dan <- matrix(rep(apply(d[, nvn, drop=FALSE], 2, mean),
                each=nrow(daf)), nrow=nrow(daf))
    colnames(dan) <- nvn;
    colnames(daf) <- nvf;
    dal <- cbind(daf, dan)[nv];

    dsg <- model.matrix(form, dal);

    lss <- lapply(lsel, function(sel, dal){
        ss <- lapply(seq_len(length(sel)), function(j, sel, dal){
            nm <- names(sel)[[j]];
            if (is.factor(dal[, nm])) paste("'", sel[[j]], "'", sep="");
        }, sel, dal);
        names(ss) <- names(sel);
        ss;
    }, dal);

    ds0 <- vapply(lss, function(ss, d){
        eval(parse(text=paste(vapply(seq_len(length(ss)), function(j, ss, d)
            {
                o <- ss[[j]];
                names(o) <- names(ss)[j];
                paste("d$", names(ss)[j], "%in%c(",
                    paste(o, collapse=', ', sep=''), ")", sep='');
            }, character(1), ss, dal), collapse=' & ')))
    }, logical(nrow(d)), d)


    ds <- vapply(lss, function(ss, dal){
        eval(parse(text=paste(vapply(seq_len(length(ss)), function(j, ss, dal)
            {
                o <- ss[[j]];
                names(o) <- names(ss)[j];
                paste("dal$", names(ss)[j], "%in%c(",
                    paste(o, collapse=', ', sep=''), ")", sep='');
            }, character(1), ss, dal), collapse=' & ')))
    }, logical(nrow(dal)), dal)

    grs <- t(apply(ds, 2, function(s, dsg)
        apply(dsg[s, , drop=FALSE], 2, mean), dsg));
    cont <- (t(grs)%*%signs)[, 1]


    if (!is.null(ngroups)) rownames(grs) <- colnames(ds0) <- ngroups;
    list(grs=grs, cont=cont, sels=ds0);
}


