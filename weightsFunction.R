
############## function
weightgeneset <- function(ysel, index, method = c("corthr","downExc","pcaCont"), nite = 1000, thr = 0.4, Nsel = 20, mc.cores = 2){
  Weight.list <- sapply(1:length(index), function(o){
     yco <- ysel[index[[o]],]

     if(method[1] == "corthr")
     {
         cor1 <- cor(t(yco))
         pcvar <- 1 / apply(abs(cor1)>thr,1,sum)
#         pcvar <- 1 / apply(abs(cor(t(yco))),1,mean)
     }
     if(method[1] == "downExc")
     {
         yco <- t(scale(t(yco)))
         ll <- mclapply(1:nite,function(i){
             yn <- yco[,sample(1:dim(yco)[2],Nsel)]
             names(which(abs(apply(yn,1,mean)) * sqrt(Nsel) > thr))
         }, mc.cores=mc.cores)
         len <- sapply(ll,length)

         pcvar <- 1 / (sapply(rownames(yco), function(o) mean(len[sapply(1:nite, function(i) any(ll[[i]]==o))])))
     }
     if(method[1] == "pcaCont")
     {
         yco <- t(scale(t(yco)))
         pca1 <- prcomp(t(yco))
         pcvar <-as.numeric( 1 /( t(abs(cor(pca1$x,t(yco)))) %*% t(t(pca1$sdev))^2))
     }

     d2 <- as.matrix(pcvar/sum(pcvar))
  })
  return(Weight.list)
}

weightgeneset2 <- function(ysel, index, method = c("corthr","downExc","pcaCont"), nite = 1000, thr = 0.4, Nsel = 20, mc.cores = 2){
      if(method[1] == "corthr")
     {
         cor1 <- cor(t(ysel))
         pcvar <- 1 / apply(abs(cor1)>thr,1,sum)
#         pcvar <- 1 / apply(abs(cor(t(yco))),1,mean)
     }
     if(method[1] == "downExc")
     {
         yco <- t(scale(t(ysel)))
         ll <- mclapply(1:nite,function(i){
             yn <- yco[,sample(1:dim(ysel)[2],Nsel)]
             names(which(abs(apply(yn,1,mean)) * sqrt(Nsel) > thr))
         }, mc.cores=mc.cores)
         len <- sapply(ll,length)

         pcvar <- 1 / (sapply(rownames(ysel), function(o) mean(len[sapply(1:nite, function(i) any(ll[[i]]==o))])))
     }
     if(method[1] == "pcaCont")
     {
         yco <- t(scale(t(ysel)))
         pca1 <- prcomp(t(ysel))
         pcvar <-as.numeric( 1 /( t(abs(cor(pca1$x,t(ysel)))) %*% t(t(pca1$sdev))^2))
     }

    Weight.list <- sapply(1:length(index), function(o){
     d2 <- as.matrix(pcvar[index[[o]]]/sum(pcvar[index[[o]]]))
  })
  return(Weight.list)
}

dweightgeneset <- function(ysel,  method = c("clust","corthr","downExc","pcaCont"),
                           nite = 1000, thr = 0.4, Nsel = 20, mc.cores = 2, ncomp = 5){
      if(method[1] == "corthr")
     {
         cor1 <- cor(t(ysel))
         pcvar <- 1 / apply(abs(cor1)>thr,1,sum)
#         pcvar <- 1 / apply(abs(cor(t(yco))),1,mean)
     }
     if(method[1] == "downExc")
     {
         yco <- t(scale(t(ysel)))
         ll <- mclapply(1:nite,function(i){
             yn <- yco[,sample(1:dim(ysel)[2],Nsel)]
             names(which(abs(apply(yn,1,mean)) * sqrt(Nsel) > thr))
         }, mc.cores=mc.cores)
         len <- sapply(ll,length)

         pcvar <- 1 / (sapply(rownames(ysel), function(o) mean(len[sapply(1:nite, function(i) any(ll[[i]]==o))])))
     }
     if(method[1] == "pcaCont")
     {
         yco <- t((t(ysel)))
         pca1 <- prcomp(t(ysel))
         pcvar <-as.numeric( 1 /( t(abs(cor(pca1$x,t(ysel)))) %*% t(t(pca1$sdev))^2))
     }
     if(method[1] == "clust"){

       if(nrow(ysel) < 20)  pcvar <- rep(1,nrow(ysel))
       else{
         al <- sapply(5:round(nrow(ysel)/3), function(s){
            kms <- kmeans(ysel,s)$clust
            sapply(1:length(kms), function(k) sum(kms==kms[k]))
         })
         pcvar <- (1/apply(al,1,mean))
     }
   }

   if(method[1] == "lm"){
         ncomp <- min(c(ncomp, ncol(ysel)-2, nrow(ysel)-2))
         al  <- sapply(1:dim(ysel)[2], function(k) (summary(lm(ysel[,k] ~ -1 + prcomp(ysel[,-k])$x[,1:ncomp])))$r.squared)
         pcvar <- (1-al)
   }

    Weight.list <- as.matrix(pcvar/sum(pcvar))
    return(Weight.list)
}
