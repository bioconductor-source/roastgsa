

###########################################
################# example #################
###########################################

### load code (in its last version) and data
rdir <- "/Volumes/biostats/acaballe/GSEA_significance/code/agsea_to_share/"
lsd <- dir(rdir)
lsd <- paste0(rdir,lsd[regexpr("~",lsd)<0])
for(A in lsd) source(A)

load(paste0("/Volumes/biostats/acaballe/GSEA_significance/testing_code/", "data_files.Rdata"))

### genesets of interest
gspath <- '/Volumes/biostats/databases/BroadGSEA/genesets/hsapiens/'
gsetsel <- "Hallmarks"

# test with two samples
#covar2 <- covar[order(rownames(covar))[1:4],]
#covar2$Codigo.Muestra <- droplevels(covar2$Codigo.Muestra)
#y2 <- y[,order(rownames(covar))[1:4]]
covar2 <- covar
y2 <- y

### formula and covariates, design matrix can also be given
form <- "~ -1 + Codigo.Muestra + Time"
a2 <- roastgsa(y2, form = form, covar = covar2, contrast = "Time6meses", gspath = gspath,
               gsetsel = gsetsel, nrot = 1000, mccores = 2, set.statistic = "maxmean",
               self.contained = FALSE)

### plotMean
plotMean(a2, whplot = rownames(a2$res)[1], maintitle = "", statistic = "maxmean", ylimAll = TRUE, ylim = NULL,
         maxDensPlot = 20, cex.lab= 1.4)

### plotGSEA
plotGSEA(a2, whplot = rownames(a2$res)[1], maintitle = "", gsainfo = TRUE, cex.sub = 0.5, lwd = 2)

### heatmap at sample level (for each geneset) according to set.statistic
aux <- heatmapMean(a2, y2, whplot = rownames(a2$res)[1], toplot = TRUE,
               mycol = c("black","orange","green","white"),
               dendrogram =  "n", col=  bluered(100),trace='none',
               margins=c(8,16),notecol='black',notecex=1, keysize=.9, cexCol=1.5,
               Rowv = TRUE, Colv = FALSE,las =2)

### heatmap at sample level (for all geneset) according to set.statistic
aux <- heatmapMean(a2, y2, whplot = rownames(a2$res)[1:30], toplot = TRUE, pathwaylevel = TRUE, fdrkey = TRUE,
               mycol = c("black","orange","green","white"),
               dendrogram =  "n", col=  bluered(100),trace='none',
               margins=c(8,20),notecol='black',notecex=1, keysize=.9, cexCol=1.5,
                   Rowv = TRUE, Colv = FALSE,las =2)

### Finding effective signature size
varrot <- effsignaturesize(a2, y2, testedsizes = c(3:30,seq(32,50, by=2), seq(55,200,by=5)), nrep = 200) # can take some time
ploteffsignaturesize(a2, varrot,  whplot = 2)

### html output (DE at gene level for each pathway can be provided at geneDEhtmlfiles)
htmlagsea(a2, htmlpath = "/Volumes/biostats/acaballe/GSEA_significance/testing_code/",
           htmlname = "fileprova.html", plotpath ="plots/",
           plotmean = TRUE, plotgsea = TRUE, indheatmap = TRUE, ploteffsize = TRUE,
           y = y2, whplots = rownames(a2$res)[1:30], tit = "", margins = c(5,16),
           geneDEhtmlfiles  = NULL, varrot = varrot)

