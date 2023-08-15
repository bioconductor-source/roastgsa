########################################################
#########  function:  ploteffsignaturesize     #########
########################################################

ploteffsignaturesize <- function(obj, varrot,  whplot = 1, ...){
    if(mode(whplot) =="character")
        whplot <- which(rownames(obj$res)%in%whplot)[1]

    if(!is(varrot, 'varrotrand'))
        stop("varrot must be of class varrotrand")
    varobs <- (obj$res[["est"]]/obj$res[["nes"]])[whplot]
    nrep <- varrot$nrep

    props <- apply(varrot$varrot,2,function(x)
        2*min(sum(x<varobs),sum(x>varobs)))/nrep
    plot(varrot$testedsizes,props,type="p",ylab="approx-pval",
        ylim =c(0,1.05), xlab="eff.size", pch=16,
        xlim=c(0,max(max(varrot$testedsizes),
        obj$res$measured_genes[whplot])+10), ...)

    abline(v= obj$res$measured_genes[whplot], col=2, lty=2)
    abline(v = varrot$testedsizes[which.max(props)], col=3, lty=3)
    text(varrot$testedsizes[which.max(props)]+1, 1,
        paste0(varrot$testedsizes[which.max(props)]), col=3)
    text(obj$res$measured_genes[whplot], 1,
        paste0("nsig = ", obj$res$measured_genes[whplot]), col=2)
}

