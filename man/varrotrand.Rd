\name{varrotrand}
\alias{varrotrand}

\title{roastgsa variance rotations under gene randomization }
\description{
Computation of the sample variance of rotation scores under gene randomization
}
\usage{

varrotrand(obj, y, testedsizes = c(3:30,seq(32,50, by=2),
    seq(55,200,by=5)), nrep = 200, nrot = NULL,
    mccores =  NULL, psel = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{obj}{an object of class 'roastgsa'}
    \item{y}{data used in \code{\link{roastgsa}} call}
    \item{testedsizes}{effective sizes to be tested}
    \item{nrep}{number of randomly selected gene sets created for each
    tested effective size}
    \item{nrot}{number of rotations used for hypothesis testing}
    \item{mccores}{the number of cores to use for parallel executions}
    \item{psel}{character vector with probesets (one per gene) to be used
    for roastgsa statistic in a microarray experiment }
}
\details{
When a specific gene that is highly correlated to the rest of the gene set
finds an extreme value, even under H\eqn{_0}, it is likely that many
other genes in the gene set follow it with large values as well. We define 
the concept  of effective signature size of a gene set by the number of 
randomly selected (not necessarily independent) genes that are needed to 
achieve comparable variability levels on rotation summary test statistics. 
This can be viewed as a realistic measure of the total number of independent 
variables that contribute to the power of the test. The function presented 
here computes the sample variance of the rotation scores in randomly generated 
signatures of several sizes. The comparison to the observed variances (using 
the testing gene sets in the \code{roastgsa} call) is done through the function 
\code{\link{ploteffsignaturesize}}. 
}
\value{
return an object of class  \code{varrotrand} with attributes

\item{"varrot"}{matrix \code{nrep} x \code{testedsizes} with the
estimated variance of the rotation scores using \code{nrot} rotations}
\item{"testedsizes"}{effective sizes being tested}
\item{"nrep"}{number of gene sets created for each tested effective size}
}

\author{ Adria Caballe Mestres }

\references{
[1] A comparison of rotation-based scores for gene set analysis
Adria Caballe Mestres, Antonio Berenguer Llergo, Camille Stephan-Otto Attolini
bioRxiv 2021.03.23.436604; 
}

\seealso{\code{\link{ploteffsignaturesize}} to visualize results and
\code{\link{roastgsa}} for gsa approach}
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
}


