\name{ploteffsignaturesize}
\alias{ploteffsignaturesize}

\title{roastgsa effective signature size }
\description{
    Approximation of effective signature size under gene randomization
}
\usage{
ploteffsignaturesize(obj, varrot,  whplot = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{obj}{an object of class 'roastgsa'}
    \item{varrot}{an object of class 'varrotrand'
    (see \code{\link{varrotrand}}) with estimated rotation score
    variances for randomly selected genesets of several sizes. }
    \item{whplot}{selected pathway. If integer vector, the pathways are
    selected in the same order as the table in \code{obj$res}}
    \item{\dots}{
    Arguments passed to or from other methods to the low level.}
}
\details{
    The plot shows the approximated probability of obtaining a test statistic 
    variance (under rotations of the residual space of the data) as extreme as 
    the observed when generating randomly gene sets of several sizes.
}
\value{
    plot object with the effective signature size representation of roastgsa 
    results
}

\author{ Adria Caballe Mestres }

\references{
[1] A comparison of rotation-based scores for gene set analysis
Adria Caballe Mestres, Antonio Berenguer Llergo, Camille Stephan-Otto Attolini
bioRxiv 2021.03.23.436604; 

}

\seealso{\code{\link{varrotrand}} and \code{\link{roastgsa}}}
\examples{
y <- array(rnorm(10000),dim = c(1000,10))
covar <- data.frame(voi = factor(c(rep(0,5),rep(1,5))))
colnames(y) <- rownames(covar) <- paste0("sample",1:10)
rownames(y) <- paste0("gene",1:1000)
form <- as.formula("~ voi")
index <- lapply(1:10, function(o) sample(1:1000,50))
names(index) <- paste0("gset",1:10)


roastgsa1 <- roastgsa(y, covar, form = form, self.contained = TRUE,
    set.statistic = "maxmean", index = index, nrot = 100,
    mccores = 1, executation.info = FALSE)

varrot <- varrotrand(roastgsa1, y,
    testedsizes = c(seq(5,50, by=5), seq(55,200,by=10)), 
	nrep = 50)

ploteffsignaturesize(roastgsa1, varrot,  whplot = 2)
}


