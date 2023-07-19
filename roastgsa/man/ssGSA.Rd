\name{ssGSA}
\alias{ssGSA}
\alias{print.ssGSA}

\title{  Single sample Gene Set Analysis }
\description{
Single sample gene set analysis using z-score summarized data for
linear model hypothesis testing
}
\usage{
ssGSA(y, obj = NULL, design = NULL, contrast = NULL, index = NULL,
    method = c("GScor","GSadj","zscore"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{y}{expression matrix with columns indicating samples and rows
    indicating genes}
    \item{obj}{object of class 'roastgsa' used to extract the design, the
    contrast and the index arguments}
    \item{design}{the design matrix of the experiment. Considered only
    if \code{obj} is NULL}
    \item{contrast}{comparison to consider in the model. Considered only
    if \code{obj} is NULL}
    \item{index}{list with index vectors specifying which rows of \code{y}
    are in the testing sets. Either integer indexes with row positions or
    gene identifiers can be stated. Considered only if \code{obj} is NULL}
    \item{method}{If \code{"GSadj"}, a correction variable with the average 
    trend in the data enters in the model as confounding variable. 
    If \code{"GScor"}, gene signatures are adjusted a priori by subtracting 
    the correction variable values. Check details for more information. }
}
\details{
A correction by the overall tendency can be done a priori (GScor) or it
can be incorporated as a covariate in the linear model (GSadj).
The correction variable used here is what we have called the global signature
(GS) of the experiment, that for each sample can be calculated as the 
average z-score of all genes measured in \code{y}. This GS corrects or centers 
global technical / sampling directions in the data.

}
\value{
return an object of class  \code{ssGSA}  with  attributes

\item{"res"}{data.frame with main results obtained in hypothesis
testing. Total genes in the gene set, the average score, the test
statistic, p-value and adjusted pvalue.}
\item{"stats"}{adjusted z-scores matrix}
}

\author{ Adria Caballe Mestres }

\references{
[1] Caballe Mestres A, Berenguer Llergo A and Stephan-Otto Attolini C.
Adjusting for systematic technical biases in risk assessment of gene signatures
in transcriptomic cancer cohorts. bioRxiv (2018).

}

\seealso{\code{\link{plot.ssGSA}} }
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

}


