\name{clusterCt}
\Rdversion{1.1}
\alias{clusterCt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Clustering of qPCR Ct values}
\description{Hierarchical clustering of samples or genes from high-throuhgput qPCR experiments, such as the TaqMan Low Density Array platform. Individual clusters can be selected, and the features within them listed in the given order.}

\usage{
clusterCt(q, main = NULL, type = "genes", dist = "pearson", xlab = "Cluster dendrogram", n.cluster, h.cluster, select.cluster = FALSE, ...)
}

\arguments{
  \item{q}{object of class qPCRset.}
  \item{main}{character string, plot title.}
  \item{type}{character string, either "genes" (default) or "samples", indicating what is to be clustered.}
  \item{dist}{character string, specifying whether to use "pearson" correlation (default) or "euclidean" distance for the clustering.}
  \item{xlab}{character string, label for the x-axis.}
  \item{n.cluster}{integer, the number of cluster to divide the dendrogram into. See details.}
  \item{h.cluster}{numeric, the height at which to cut the dendrogram into clusters. See details.}
  \item{select.cluster}{logical, whether to select clusters interactively. See details.}
  \item{\dots}{any other arguments will be passed to the \code{plot} function.}
}

\details{This function may be used to cluster the Ct values and present the result as a dendrogram.

The \code{n.cluster} and \code{h.cluster} parameters are from the \code{rect.hclust} function and can be used to divide the dendrogram into subclusters based on either number of clusters or height of branch, drawing boxes around subclusters. The members of each cluster can be returned (see value). If \code{n.cluster} is specified \code{h.cluster} will be ignored.

If \code{select.cluster} is chosen individual subclusters can be selected and marked by a box by clicking on their highest comment branch with the (first) mouse buttom. Multiple clusters can be selected until any mouse button other than the first is pressed, and the function can be used in conjunction with either \code{n.cluster} or \code{h.cluster}. The members of each cluster will likewise be returned, in the order they were selected.}

\value{A plot is created on the current graphics device. If any subclusters are marked, these will be returned invisibly in a list, with one component for each subcluster. The individual slots in the list contain the names og the genes, and their position in the original input data (row number).}

\author{Heidi Dvinge}

\seealso{\code{\link{hclust}}, \code{\link{dist}}, \code{\link{rect.hclust}}, \code{\link{identify.hclust}}}

\examples{
# Load example data
data(qPCRraw)
# Clustering samples
clusterCt(qPCRraw, type="samples")
clusterCt(qPCRraw, type="samples", dist="euclidean")
# Clustering genes
clusterCt(qPCRraw, type="genes", cex=0.5)
clusterCt(qPCRraw, type="genes", h.cluster=1.5, cex=0.5)
cluster.list <- clusterCt(qPCRraw, type="genes", n.cluster=6, cex=0.5)
cluster.list[[1]]
}

\keyword{hplot}
