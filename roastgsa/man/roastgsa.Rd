\name{roastgsa}
\alias{roastgsa}
\alias{print.roastgsa}

\title{  Rotation-based Gene Set Analysis }
\description{
Gene set analysis using rotations for hypothesis testing. Test statistic
options include KS-based statistics used in GSEA or GSVA as well as
summary statistics such as mean, maxmean, median, absmean and mean.rank
}

\usage{

roastgsa(y, covar, form, contrast = NA, design = NULL, gsetsel,
    gspath, index = NULL, self.contained = FALSE,
    set.statistic = "maxmean",  psel = NULL, nrot = 9999,
    minsize = 10, maxsize = 500, mccores = 1,
    executation.info = TRUE, weights = NULL, shrink.resid = TRUE,
    normalizeScores = TRUE, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{y}{expression matrix with columns indicating samples and rows
    indicating genes}
    \item{covar}{data frame with the covariates}
    \item{form}{description of the model to be fitted}
    \item{contrast}{comparison to consider in the model. If NA, the
    last column of the design matrix is used}
    \item{design}{the design matrix of the experiment. If null, this is
    calculated using the \code{form} and the \code{covar} arguments}
    \item{gsetsel}{character string with gene set database to be used in
    format \code{.gmt}. If missing, index argument has to be provided}
    \item{gspath}{path for the gene set database}
    \item{index}{list with index vectors specifying which rows of \code{y}
    are in the testing sets. Either integer indexes with row positions or
    gene identifiers can be stated. If NULL, the index is computed using
    information in the \code{gsetsel} and \code{gspath}
    arguments}
    \item{self.contained}{competitive test (\code{FALSE}) or self contained test
    (\code{TRUE})}
    \item{set.statistic}{to be chosen from \code{"maxmean"} (default),
    \code{"mean"},  \code{"mean.rank"}, \code{"median"},
    \code{"absmean"}, \code{"GSEA"} and \code{"GSVA"}}
    \item{psel}{character vector with probesets (one per gene) to be used
    for roastgsa statistic in a microarray experiment }

    \item{nrot}{number of rotations used for hypothesis testing }
    \item{minsize}{minimum size of the testing sets allowed for hypothesis
    testing }
    \item{maxsize}{maximum size of the testing sets allowed for hypothesis
    testing }
    \item{mccores}{the number of cores to use for parallel executions}
    \item{executation.info}{Show (if set to TRUE) the progress-bar of the
    iterative process }
    \item{weights}{list with the gene weights in each testing set. Only
    for set.statistic =  \code{"maxmean"} and \code{"mean"}. If
    \code{NULL}, weights are assumed to be constant}

    \item{shrink.resid}{if \code{TRUE}, the coefficients of the linear
    model are shrunk towards zero for rotations to increase the power}

    \item{normalizeScores}{transform the moderated t-statistics to
    z-scores}
    \item{\dots}{
    Arguments passed to or from other methods to the low level.}
}
\details{
We consider 7 different enrichment score functions which we refer by the
names of mean, maxmean, median, absmean, mean.rank, GSEA and GSVA.
The first four functions (mean, maxmean, median, absmean) are formulated
for the two type of testing problems (self-contained and
competitive). The mean.rank, GSEA and GSVA are exclusive scores for
the competitive approach. The absmean is a non-directional score that can
be used to give priority to gene sets with both activator and inhibitor
genes. The mean is a democratic score that gives priority to detecting
gene sets in which a large fraction of their genes present similar effect sizes
going at the same direction. The maxmean (default) falls in between the mean
and the absmean scores, being capable to recover both type of gene sets
consistently.

Some of the defined sets are composed by genes that interact together in
any particular biological condition, leading to intra-gene set
correlation structures with high levels of correlation. We encourage the
usage of effective signatures size, that can be a proxy for the number
of uncorrelated genes in the gene set used for GSA
(\code{\link{varrotrand}} and
\code{\link{ploteffsignaturesize}}). Through the argument
\code{weights}, we provide the possibility to redefining the gene set by
weighting the importance of each gene in the list.

GSEA and GSVA scores are computationally much more intensive than the
other scores.
}
\value{
return an object of class  \code{roastgsa} with attributes

\item{"res"}{data.frame with main results obtained in hypothesis
testing. Total genes in the geneset, the number of genes also in the
\code{y}, the test statistic, the normalized score and the
significance of the tests }
\item{"stats"}{Moderated t-statistics for all genes}
\item{"contrast"}{contrast used in a vector form}
\item{"index"}{list with gene set symbols }
}

\author{ Adria Caballe Mestres }

\references{
[1] A comparison of rotation-based scores for gene set analysis
Adria Caballe Mestres, Antonio Berenguer Llergo, Camille Stephan-Otto Attolini
bioRxiv 2021.03.23.436604; doi:

[2] E. Lim, D. Wu, G. K. Smyth, M.-L. Asselin-Labat, F. Vaillant, and
J. E. Visvader. ROAST: rotation gene set tests for complex microarray
experiments. Bioinformatics, 26(17):2176-2182, 2010.

}

\seealso{\code{\link{roast}} }
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
print(roastgsa1)
}


