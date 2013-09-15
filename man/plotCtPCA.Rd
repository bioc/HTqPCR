\name{plotCtPCA}
\Rdversion{1.1}
\alias{plotCtPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{PCA for qPCR Ct values.}
\description{Perform and plot a principal component analysis for high-throughput qPCR data from any platform, for doing clustering.}

\usage{
plotCtPCA(q, s.names, f.names, scale = TRUE, features = TRUE, col, cex = c(1, 1))
}

\arguments{
  \item{q}{a matrix or an object of class qPCRset containing Ct values.}
  \item{s.names}{character vector, names of samples. See details.}
  \item{f.names}{character vector, names of features. See details.}
  \item{scale}{logical, should the variables be scaled to have unit variance. Passed on to \code{prcomp}.}
  \item{features}{logical, should the features be plotted. See details.}
  \item{col}{vector, the colours to use for the samples if \code{features=FALSE}.}
  \item{cex}{vector of length 2, the expansion to use for features and samples respectively if \code{features=FALSE}.}
}

\details{Per default the sample names from the qPCRset are used, however the feature names are replaced by "*" to avoid cluttering the plot.
	
If \code{features=TRUE} then a biplot including all features is produced, with samples  represented by vectors. I.e. both observations and variables are plotted, which can potentially be used to identify outliers among the features. For \code{features=FALSE} only the samples will be included in the plot. This might be more useful for clustering.
	
In case of high-throughput arrays, some samples may be all NAs. These are ignored during the PCA calculation.}

\value{A plot is created on the current graphics device.}

\author{Heidi Dvinge}

\section{Note }{This is still a work in progress, and the function is not particularly sophisticated. Suggestions/wishes are welcome though.}

\seealso{\code{\link{prcomp}}, \code{\link{biplot}}}

\examples{
# Load example data
data(qPCRraw)
# Plot
plotCtPCA(qPCRraw)
# Include feature names; make them smaller
plotCtPCA(qPCRraw, f.names=featureNames(qPCRraw), cex=c(0.5,1))
# Plot only the samples
plotCtPCA(qPCRraw, features=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
