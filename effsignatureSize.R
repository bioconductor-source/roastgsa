
########################################################
#########  function:  plot adjheatmapMean  #############
########################################################

effsignaturesize <- function(obj, y, testedsizes = c(3:30,seq(32,50, by=2), seq(55,200,by=5)),
                             nrep = 200){
    k <- 1
    varrot <- sapply(testedsizes, function(k){
        print(k)
        indexaux <- lapply(1:nrep, function(o) sample(rownames(y), k))
        obj2 <- roastgsa(y, form = obj$form, covar = obj$covar, contrast = obj$contrast,
                     index = indexaux, nrot = attr(obj,"nrot"), mccores = attr(obj,"mccores"),
                     set.statistic = obj$statistic, self.contained = obj$self.contained)
        obj2$res[["est"]]/obj2$res[["nes"]]
    })
    varrot <- list(varrot = varrot, testedsizes = testedsizes,nrep=nrep)
    class(varrot) <- "effsignaturesize"
    return(varrot)
}

ploteffsignaturesize <- function(obj, varrot,  whplot = 1){
    if(mode(whplot) =="character") whplot <- which(rownames(obj$res)%in%whplot)[1]

    if(class(varrot)!="effsignaturesize")
        stop("varrot must be of class effsignaturesize")
    varobs <- (obj$res[["est"]]/obj$res[["nes"]])[whplot]
    nrep <- varrot$nrep

    props <- apply(varrot$varrot,2,function(x) 2*min(sum(x<varobs),sum(x>varobs)))/nrep
    plot(varrot$testedsizes,props,type="p",ylab="approx-pval", ylim =c(0,1.05), xlab="eff.size", pch=16,
         xlim=c(0,max(max(varrot$testedsizes),obj$res$overlap_genes[whplot])+10))

    abline(v= obj$res$overlap_genes[whplot], col=2, lty=2)
    abline(v = varrot$testedsizes[which.max(props)], col=3, lty=3)
    text(varrot$testedsizes[which.max(props)]+1, 1, paste0(varrot$testedsizes[which.max(props)]), col=3)
    text(obj$res$overlap_genes[whplot], 1, paste0("nsig = ", obj$res$overlap_genes[whplot]), col=2)
}

