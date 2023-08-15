\name{plotGSEA}
\alias{plotGSEA}

\title{  GSEA plot  }
\description{
    GSEA plot for \code{roastgsa} objects
}
\usage{

plotGSEA(obj, whplot = 1, maintitle = "", gsainfo = TRUE, cex.sub = 0.8,
    lwd = 2, ...)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{obj}{an object of class \code{'roastgsa'}}
    \item{whplot}{selected pathway. If integer vector, the pathways are
    selected in the same order as observed in the \code{obj$res} table}
    \item{maintitle}{plot main title. If \code{maintitle == ""}, the
    name of the pathway in
    \code{obj} is printed}
    \item{gsainfo}{if \code{TRUE}, the subtitle shows the GSA main results}
    \item{cex.sub}{cex for subtitle}
    \item{lwd}{line width}
    \item{\dots}{
    Arguments passed to or from other methods to the low level.}
}
\details{
Standard representation of Kolmogorov-Smirnov GSEA enrichment score.

}
\value{
plot object with the GSEA representation of roastgsa results
}

\author{ Adria Caballe Mestres }

\references{
Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, 
Paulovich A, Pomeroy SL, Golub TR, Lander ES, Mesirov JP. Gene set 
enrichment analysis: A knowledge-based approach for interpreting 
genome-wide expression profiles. PNAS. 2005;102(43);15545-15550.  
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
plotGSEA(roastgsa1, whplot = 1, gsainfo =TRUE, maintitle =
    "", statistic = "mean")
}


