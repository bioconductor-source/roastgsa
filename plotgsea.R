
##############################################
#########  function:  plot GSEA  #############
##############################################

library(RColorBrewer)
plotGSEA <- function(obj, whplot = 1, maintitle = "", gsainfo = TRUE, cex.sub = 0.8, lwd = 2, ...){
  res <- obj$res
  index <- obj$index[rownames(res[whplot,])]
  stats <- obj$stats

  modt2 <- stats[order(-stats)]
  if(is.list(index)) index2 <- lapply(index, function(o) which(names(modt2)%in%o))
  else index2 <-  list(which(names(modt2)%in%index))
  def.par <- par(no.readonly = TRUE)

  for(K in names(index)){
      if(maintitle == ""){
          maintitle <- K
          if(gsainfo) subtitle  <- paste(names(res[K,]),round( as.numeric(res[K,]),3),sep=" = ",collapse="; " )
      }
      else if(gsainfo) subtitle <- ""
      par(mar=c(1,4,5,1))
      layout(c(1,2),heights=c(4,2))

      es <- sapply(index2, function(o, modt2) enrichmentScore(modt2, o), modt2)
      plot(es[,K],type='l',axes=FALSE, ylab="ES", xlab='', oma=c(0,0,0,0), col='darkgreen',lwd=2, main = maintitle)
      mtext(subtitle, 3, line=0.5, col= "darkgray", cex = cex.sub)
      abline(h=0)
      axis(2)
      abline(h=0)
      abline(v= which.max(abs(es[,K])), col=2,lty=3)

      par(mar=c(5,4,1,1))
      plot(NA,NA,xlim=c(1,length(es[,K])),ylim=c(0,1),axes=FALSE,ylab='',xlab='Gene list rank',
           sub="",col='green',lwd=2)
      myTicks <- axTicks(1); myTicks[myTicks==0] <- 1; axis(1,lty=0,at=myTicks)
      s <- rep(FALSE, length(es[,K])); s[index2[[K]]] <- TRUE
      if (sum(s)<50) myCols <- "darkgreen"
      else myCols <- densCols(which(s[1:length(s)]),colramp=colorRampPalette(brewer.pal(9, "Greens")[-c(1:3)]))
      abline(v=which(s[1:length(s)]),col=myCols)
      polygon(c(1,length(s)/15,length(s)/15,1), c(0-0.05,0-0.05,1+0.05,1+0.05),
              col= rgb(0.5, 0, 0,0.3), border = NA)
      polygon(c(length(s), length(s) - length(s)/15, length(s) - length(s)/15, length(s)),
              c(0-0.05,0-0.05,1+0.05,1+0.05), col= rgb(0,0,.5,0.3), border = NA)
      text(length(s)/30, 0.5, "+", cex=2)
      text(length(s)-length(s)/30, 0.5, "-", cex=1.5)
  }
  par(def.par)
}
