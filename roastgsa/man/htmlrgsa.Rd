\name{htmlrgsa}
\alias{htmlrgsa}

\title{roastgsa results in html form}
\description{
    Writing html document with roastgsa output
}
\usage{

htmlrgsa(obj, htmlpath = "", htmlname = "file.html", plotpath ="",
    plotstats = TRUE, plotgsea = TRUE, indheatmap = TRUE, ploteffsize = TRUE,
    links_plots = list(stats= NULL, gsea = NULL, heatmap = NULL, effsize
    = NULL), y, whplots = NULL, geneDEhtmlfiles  = NULL, tit = "",
    margins = c(5,16), sizesHeatmap = c(1200, 800), typeheatmap =
    c("heatmap.2","ggplot2"), intvar, adj.var = NULL, mycol, varrot,
    psel = NULL, sorttable, dragtable,  ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{obj}{an object of class 'roastgsa'}

    \item{htmlpath}{path for html file to be placed}

    \item{htmlname}{name of html file}

    \item{plotpath}{added path from argument \code{htmlpath} where plots
    should be saved}

    \item{plotstats}{plots using \code{\link{plotStats}} are created}

    \item{plotgsea}{plots using \code{\link{plotGSEA}} are created}

    \item{indheatmap}{plots using \code{\link{heatmaprgsa_hm}} at the gene
    level are created}

    \item{ploteffsize}{plots using \code{\link{ploteffsignaturesize}} are
    created}

    \item{links_plots}{list with 4 elements  (stats, gsea, heatmap and
    effsize) specifying the path of all plots (paths set from htmlpath)
    in case these were already created. If NULL, links are obtained from
    \code{plotpath} if any of \code{plotstats}, \code{plotGSEA},
    \code{indheatmap} or \code{ploteffsize} is \code{TRUE}}

    \item{y}{data used for \code{\link{roastgsa}}}

    \item{whplots}{selected pathways. If integer vector, the pathways are
    selected in the same order as the table in obj$res. If null all
    tested pathways are selected }

    \item{geneDEhtmlfiles}{vector with links to html-tables showing the 
    differential expression results for the subsets of genes determined by 
    \code{whplots}}

    \item{tit}{title of the html file}

    \item{margins}{margins for the heatmap plots}

    \item{sizesHeatmap}{vector with two elements providing png sizes
    (width, height)}

    \item{typeheatmap}{either ggplot2 type or heatmap.2 type}

    \item{intvar}{for \code{\link{heatmaprgsa}}. Name of variable of
    interest in \code{obj$formula}. If
    missing, last term of \code{obj$formula} is used}

    \item{adj.var}{for \code{\link{heatmaprgsa}}. Name of covariates in
    \code{obj$design} to   adjust using  a linear model in the heatmap
    representation. If \code{NULL} no prior adjustment applies}

    \item{mycol}{color for heatmap columns defining the groups of the
    variable of interest}

    \item{varrot}{an object of class 'varrotrand'
    (see \code{\link{varrotrand}}) with estimated rotation score
    variances for randomly selected genesets of several sizes. Cannot be
    missing if \code{ploteffsize = TRUE}}

    \item{psel}{character vector with probesets (one per gene) to be used
    for roastgsa statistic in a microarray experiment }

    \item{sorttable}{internal data loaded with \code{roasgsa} package.
    Permits sorting columns in html tables.}
    \item{dragtable}{internal data loaded with \code{roasgsa} package.
    Permits dragging elements in html tables.}

    \item{\dots}{
    Arguments passed to or from other methods to the low level.}
}
\details{
    This function permits to explore a html-table with the statistical
    results and graphical representation of the top gene sets obtained
    from an object of class \code{roastgsa}.
    
    By default four plots are considered for each gene set of interest:
    \code{\link{plotStats}}, \code{\link{plotGSEA}},
    \code{\link{heatmaprgsa_hm}} and \code{\link{ploteffsignaturesize}}.
    The first three can be computed from the 'roastgsa' object, whereas
    for \code{ploteffsignaturesize}, an object of class 'varrotrand'
    (see \code{\link{varrotrand}}) with the estimated rotation score
    variances for randomly selected gene sets of several sizes has to
    be defined at first.

}
\value{
    It saves an html table with the main results of the roastgsa
    hypothesis tesing.
}

\author{ Adria Caballe Mestres }

\references{

[1] A comparison of rotation-based scores for gene set analysis
Adria Caballe Mestres, Antonio Berenguer Llergo, Camille Stephan-Otto Attolini
bioRxiv 2021.03.23.436604; 
}

\seealso{\code{\link{roastgsa}} }
\examples{
data(sorttable)
data(dragtable)

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


#htmlrgsa(roastgsa1, htmlpath = "", htmlname = "file.html", plotpath ="",
#          plotstats = FALSE, plotgsea = FALSE, indheatmap = FALSE,
#          ploteffsize = FALSE, links_plots = list(stats= NULL, gsea = NULL,
#          heatmap = NULL, effsize = NULL), y = y, sorttable = sorttable,
#          dragtable = dragtable)

}


