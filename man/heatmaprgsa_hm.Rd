\name{heatmaprgsa_hm}
\alias{heatmaprgsa_hm}

\title{ Heatmap of roastgsa results  }
\description{
Heatmap showing sample variation for either genes
(in a particular gene set) or summarized gene signatures.
}
\usage{

heatmaprgsa_hm(obj, y, intvar, adj.var = NULL, whplot = 1, toplot = TRUE,
    pathwaylevel = FALSE, mycol =  c("black","orange","green","white"),
    sample2zero = FALSE, rgsa.like=FALSE, psel = NULL,
    dendrogram =  "n", col=  bluered(100), trace='none',
    notecol='black', notecex=1, keysize=.9,
    cexCol=1.5, Rowv = NULL, Colv = FALSE, las =2, fdrkey = FALSE,
    quantile.sat = 0.95, order1= NULL, order2 = NULL, sizex =8, sizey =5,  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{obj}{an object of class 'roastgsa'}
    \item{y}{data used for roastgsa}

    \item{intvar}{name of variable of
    interest in \code{obj$formula}. If
    missing, last term of \code{obj$formula} is used}

    \item{adj.var}{name of covariates in \code{obj$design} to
    adjust using  a linear model in the heatmap
    representation. If \code{NULL} no prior adjustment applies}

    \item{whplot}{selected pathway. If integer vector, the pathways are
    selected in the same order as the table in \code{obj$res} }

    \item{toplot}{whether to plot the heatmap or just return the adjusted
    expression matrix}

    \item{pathwaylevel}{If TRUE, the heatmap shows the variation at the
    pathway level. Otherwise,  the heatmap shows the variation
    of all genes in the selected pathways. }

    \item{mycol}{color for heatmap columns defining the groups of the
    variable of interest}

    \item{sample2zero}{Only applicable for \code{obj$statistic =
    "maxmean"}. If TRUE, expression of genes, whose moderated-t sign is
    contrary to the roastgsa score, is set to zero for all samples (as
    part of the maxmean strategy)}

    \item{rgsa.like}{apply roastgsa transformations of data
    (restandarization and set.statistic operations) samplewise (see
    details below).}

    \item{psel}{character vector with probesets (one per gene) to be used
    for roastgsa statistic in a microarray experiment }

    \item{dendrogram}{\code{\link{heatmap.2}} parameter. Character string
    indicating whether to draw 'none', 'row', 'column' or 'both'
    dendrograms. Defaults to 'n'}

    \item{col}{\code{\link{heatmap.2}} parameter. Colors used for the image}

    \item{trace}{\code{\link{heatmap.2}} parameter. Character string indicating
    whether a solid "trace" line should be drawn across 'row's or down
    'column's, 'both' or 'none'. The distance of the line from the center
    of each color-cell is proportional to the size of the measurement.}

    \item{notecol}{\code{\link{heatmap.2}} parameter. Color of note}

    \item{notecex}{\code{\link{heatmap.2}} parameter. Size of note}

    \item{keysize}{\code{\link{heatmap.2}} parameter. Numeric value
    indicating the size of the key}

    \item{cexCol}{\code{\link{heatmap.2}} parameter. Cex.axis in for the
    column axis labeling}

    \item{Rowv}{\code{\link{heatmap.2}} parameter. Determines if and how
    the row dendrogram should be reordered}

    \item{Colv}{\code{\link{heatmap.2}} parameter. Determines if and how
    the col dendrogram should be reordered}

    \item{las}{orientation of x axis}

    \item{fdrkey}{if \code{TRUE}, the BH adjusted p-value for every
    pathway tested is printed in the plot. Only considered when
    \code{pathwaylevel = TRUE} }

    \item{quantile.sat}{numeric between 0.5 and 1 used to saturate high
    values at such specified quantile (used to avoid extreme values in
    the visualization)}

    \item{order1}{
    genes order. If NULL its ordered based on the moderated-t statistics}

    \item{order2}{
    samples order. If NULL its ordered using the information of the
    variable of interest.}

    \item{sizex}{
    size of x axis}

    \item{sizey}{
    size of y axis}

    \item{\dots}{
    Arguments passed to or from other methods to the low level.}
}
\details{

This heatmap considers \eqn{n+1} columns (\eqn{n} being the sample size). The
first column represents the moderated-t statistic (or a
restandarization of the same in case of competitive testing). The
other columns confine the expression data scaled by the standard error of
the estimated coefficient in the model and centered (if
\code{rgsa.like = TRUE}). In such case, the cross product of
all data columns and the design matrix equals the first column of the
heatmap, and the  average of the first column of the heatmap equals
the observed roastgsa test statistic (at least when the set.statistic
used is either \code{mean} or \code{maxmean}).

}
\value{
a data.frame object with source data for heatmap representation
}

\author{ Adria Caballe Mestres }

\references{
[1] A comparison of rotation-based scores for gene set analysis
Adria Caballe Mestres, Antonio Berenguer Llergo, Camille Stephan-Otto Attolini
bioRxiv 2021.03.23.436604;  
}

\seealso{\code{\link{roastgsa}} and \code{\link{plotStats}} and
\code{\link{plotGSEA}} }
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


heatmaprgsa_hm(roastgsa1, y, intvar = "voi", whplot = 1, toplot = TRUE,
    pathwaylevel = FALSE, mycol = c("black","orange","green","white"),
    sample2zero = FALSE)

heatmaprgsa_hm(roastgsa1, y, intvar = "voi", whplot = 1:10, toplot = TRUE,
    pathwaylevel = TRUE, mycol = c("black","orange","green","white"),
    sample2zero = FALSE)
}


