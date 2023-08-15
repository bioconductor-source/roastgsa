\name{plot.ssGSA}
\alias{plot.ssGSA}

\title{ Plot single sample Gene Set Analysis }
\description{
    Scatter plot of single sample z-score summarized data
}
\usage{

\method{plot}{ssGSA}(x, orderby, whplot = 1, col = "black", samplename =FALSE,
    maintitle = "", ssgsaInfo = TRUE, cex.sub = 0.8, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{x}{object of class 'ssGSA'}
    \item{orderby}{numeric or factor vector of the same size and order of
    data columns used  for ssGSA. It sets the x-axis of the plot}
    \item{whplot}{selected pathway. If integer vector, the pathways are
    selected in the same order as the table in \code{x$res} }
    \item{col}{color of scatterplot points}
    \item{samplename}{whether to show or not the names of the samples
    instead of points}
    \item{maintitle}{plot main title. If \code{maintitle = ""}, the name
    of the pathway in  \code{obj} is printed}
    \item{ssgsaInfo}{if TRUE, the subtitle shows the ssGSA results}
    \item{cex.sub}{cex for subtitle}

    \item{\dots}{
    Arguments passed to or from other methods to the low level.}
}
\details{
    This graphic is a great alternative to explore gene set variation at
    sample level. This is sometimes ignored when doing GSEA, where classic
    representations (e.g., {\code{\link{plotGSEA}}}) show gene variation
    after averaging out the sample differences within each experimental
    condition.
}
\value{
    plot object with the graphical representation of ssGSA results
}

\author{ Adria Caballe Mestres }

\references{
[1] Caballe Mestres A, Berenguer Llergo A and Stephan-Otto Attolini C.
Adjusting for systematic technical biases in risk assessment of gene signatures
in transcriptomic cancer cohorts. bioRxiv (2018).

}

\seealso{\code{\link{ssGSA}} }

\examples{
y <- array(rnorm(10000),dim = c(1000,10))
covar <- data.frame(voi = factor(c(rep(0,5),rep(1,5))))
colnames(y) <- rownames(covar) <- paste0("sample",1:10)
rownames(y) <- paste0("gene",1:1000)
form <- as.formula("~ voi")
index <- lapply(1:10, function(o) sample(1:1000,50))
names(index) <- paste0("gset",1:10)
design <- model.matrix(form, covar)


ssgsa1 <- ssGSA(y, obj=NULL, design = design, contrast = 2, index = index,
    method = c("GScor"))
plot(ssgsa1, orderby = covar$voi, whplot = 1 )
}


