\name{plot.roastgsa}
\alias{plot.roastgsa}

\title{ roastgsa plot  }
\description{
    Plot for \code{roastgsa} objects
}
\usage{

\method{plot}{roastgsa}(x, type = c("stats","GSEA"), whplot = 1,
    maintitle = "", gsainfo = TRUE, cex.sub = 0.8, lwd = 2, ...)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{x}{an object of class \code{'roastgsa'}}
    \item{type}{plot type, either 'stats' or 'GSEA'}
    \item{whplot}{selected pathway. If integer vector, the pathways are
    selected in the same order as observed in the \code{obj$res} table}
    \item{maintitle}{plot main title. If \code{maintitle == ""}, the
    name of the pathway in
    \code{obj} is printed}
    \item{gsainfo}{if \code{TRUE}, the subtitle shows the GSA main results}
    \item{cex.sub}{cex for subtitle}
    \item{lwd}{line width}
    \item{\dots}{
    Arguments passed to or from other methods to the low level.
    }
}
\details{
Details for using 'type = stats' in the plot are given in
\code{\link{plotStats}}.
Details for using 'type = GSEA' in the plot are given in \code{\link{plotGSEA}}.
}

\value{
    plot object with the graphical representation of \code{roastgsa} results.
}

\author{ Adria Caballe Mestres }

\references{

[1] A comparison of rotation-based scores for gene set analysis
Adria Caballe Mestres, Antonio Berenguer Llergo, Camille Stephan-Otto Attolini
bioRxiv 2021.03.23.436604;

}

\seealso{\code{\link{roastgsa}} and \code{\link{plotStats}} }
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
plot(roastgsa1, type = "stats", whplot = 1, gsainfo =TRUE, maintitle =
    "", statistic = "mean")

}


