\name{plotStats}
\alias{plotStats}

\title{ General GSA plot  }
\description{
    General gene set analysis plot showing the ordered moderated-t
    statistics for the selected pathway
}
\usage{

plotStats(obj, whplot = 1, maintitle = "", statistic = "mean",
    ylimAll = TRUE, ylim = NULL,  minpointsDens = 20,
    gsainfo = TRUE, cex.sub = 0.8, lwd = 2, ...)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{obj}{an object of class 'roastgsa'}
    \item{whplot}{selected pathway. If integer vector, the pathways are
    selected in the same order as the table in \code{obj$res} }

    \item{maintitle}{plot main title. If \code{maintitle = ""}, the name
    of the pathway in  \code{obj} is printed}
    \item{statistic}{to be selected from 'mean' or 'median'}
    \item{ylimAll}{y limits are found using data from all genesets (if
    TRUE) or using data from only the plotted geneset (if FALSE). Only
    if \code{ylim = NULL}}
    \item{ylim}{vector of size two with y limits}
    \item{minpointsDens}{minum number of genes needed to draw the density plot}

    \item{gsainfo}{if TRUE, the subtitle shows the enrichment results}
    \item{cex.sub}{cex for subtitle}
    \item{lwd}{line width}
    \item{\dots}{
    Arguments passed to or from other methods to the low level.}
}
\details{
The \code{statistic} argument is used for competitive testing computations of
restandardized moderated-t statistics. If \code{"median"}, the median
of all stats is used for centering and the median absolute deviation
is used for scaling. If \code{"mean"}, standard normalization
applies.

It shows the ordered moderated t-statistics in various formats, area
for up- and down- expressed genes, barcode plot for these ordered
values and density.
}
\value{
plot object with a general representation of roastgsa results
}

\author{ Adria Caballe Mestres }

\references{
[1] A comparison of rotation-based scores for gene set analysis
Adria Caballe Mestres, Antonio Berenguer Llergo, Camille Stephan-Otto Attolini
bioRxiv 2021.03.23.436604; 
}

\seealso{\code{\link{roastgsa}} and \code{\link{plotGSEA}} }
\examples{
y <- array(rnorm(10000),dim = c(1000,10))
covar <- data.frame(voi = factor(c(rep(0,5),rep(1,5))))
colnames(y) <- rownames(covar) <- paste0("sample",1:10)
rownames(y) <- paste0("gene",1:1000)
form <- as.formula("~ voi")
index <- lapply(1:10, function(o) sample(1:1000,50))
names(index) <- paste0("gset",1:10)


roastgsa1 <- roastgsa(y, covar, form = form, self.contained = TRUE,
    set.statistic = "maxmean", index = index, nrot = 200,
    mccores = 1, executation.info = FALSE)
plotStats(roastgsa1, whplot = 1, maintitle = "general plot", statistic =
"mean")
}


