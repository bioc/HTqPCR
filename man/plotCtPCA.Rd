\name{plotCtPCA}
\Rdversion{1.1}
\alias{plotCtPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{PCA for qPCR Ct values.}
\description{Perform and plot a principal component analysis for high-throughput qPCR data from any platform, for doing clustering.}

\usage{
plotCtPCA(q, s.names, g.names)
}

\arguments{
  \item{q}{a matrix or an object of class qPCRset containing Ct values.}
  \item{s.names}{character vector, names of samples. See details.}
  \item{g.names}{character vector, names of genes. See details.}
}

\details{Per default the sample names from the qPCRset are used, howeer the gene names are replaced by "*" to avoid cluttering the plot.}

\value{A plot is created on the current graphics device.}

\author{Heidi Dvinge}


\section{Note }{This is still a work in progress, and the function is not particularly sophisticated.}

\seealso{\code{\link{prcomp}}}

\examples{
# Load example data
data(qPCRraw)
# Plot
plotCtPCA(qPCRraw)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
