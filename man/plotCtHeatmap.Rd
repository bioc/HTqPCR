\name{plotCtHeatmap}
\Rdversion{1.1}
\alias{plotCtHeatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Heatmap of qPCR Ct values.}
\description{Function for drawing a heatmap of Ct values from high-throughput qPCR experiments such as using TaqMan Low Density Arrays.}

\usage{
plotCtHeatmap(q, main = NULL, col, col.range, dist = "pearson", zero.center, mar, gene.names, sample.names, ...)
}

\arguments{
  \item{q}{object of class qPCRset.}
  \item{main}{character string, plot title.}
  \item{col}{the colours to use. See details.}
  \item{col.range}{vector, the range of colours to use.}
  \item{dist}{character string, specifying whether to use "pearson" correlation (default) or "euclidean" distance for the clustering.}
  \item{zero.center}{logical, should the colours be shifted to be zero-centered. See details.}
  \item{mar}{vector, the size of the borrom and right hand side margins.}
  \item{gene.names}{character vector, names to replace the genes (rows) with. See details.}
  \item{sample.names}{character vector, names to replace the samples (columns) with. See details.}
  \item{\dots}{any other arguments will be passed to the \code{heatmap.2} function.}
}

\details{This function may be used to cluster the raw or normalized Ct values, and present the result in a heatmap. 
	
The color range is used to represent the range of values for the statistic. If \code{col==NULL} the colour will be set to a spectrum from red to blue/purple, unless there are negative values in which case it goes red-yellow-green to reflect up and down regulation of genes. If \code{zero.center=NULL} then zero.center will automatically be set to TRUE to make the colour scale symmetric around 0.

%If \code{mar=NULL} the margins will automatically be set to fit with the width of the column and row names of the data.

Especially gene names will often not be readable in a standard size plotting device, and might therefore be removed. If \code{gene.names} or \code{sample.names} is set to a single character (such as "" for no naming), then this character will be repeated for all rows or columns. 
	
A standard heatmap is drawn, but this can be modified extensively using the arguments available in the \code{heatmap.2} function.}

\value{A plot is created on the current graphics device.}

\author{Heidi Dvinge}

\seealso{\code{\link[gplots]{heatmap.2}}}

\examples{
# Load example data
data(qPCRraw)
# Some standard heatmaps
plotCtHeatmap(qPCRraw, gene.names="")
plotCtHeatmap(qPCRraw, gene.names="", dist="euclidean", col.range=c(10,35))
plotCtHeatmap(qPCRraw, gene.names="", dist="euclidean", col=colorRampPalette(rev(brewer.pal(9, "YlGnBu")))(20))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
