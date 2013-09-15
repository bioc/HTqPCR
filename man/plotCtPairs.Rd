\name{plotCtPairs}
\Rdversion{1.1}
\alias{plotCtPairs}
\title{Pairwise scatterplot of multiple sets of Ct values from qPCR data.}

\description{
Produces a plot of high-throughput qPCR Ct values from N number of samples plotted pairwise against each other in an N by N plot. The Ct values will be in the upper triangle, and the correlation between samples in the lower. Features can be marked based on for example feature class or type.
}

\usage{
plotCtPairs(q, cards = TRUE, lower.panel = panel.Ct.cor, upper.panel = panel.Ct.scatter, Ct.max = 35, col = "type", pch = 20, cex.cor = 2, cex.pch = 1, diag = TRUE, ...)
}

\arguments{
  \item{q}{object of class qPCRset. }
  \item{cards}{vector, the cards to plot against each other.}
  \item{lower.panel}{function, to use for plotting the lower triangle.}
  \item{upper.panel}{function, to use for plotting the upper triangle.}
  \item{Ct.max}{numeric, Ct values above this limit will be excluded when calculating the correlation.}
  \item{col}{vector with the colour(s) to use for the points, or a character string ("type" or "class") indicating whether points should be coloured according to \code{featureType} or \code{featureClass} of \code{q}.}
  \item{pch}{integer or single character, which plotting symbol to use for the points.}
  \item{cex.cor}{numeric, the expansion factor for the text in \code{panel.Ct.cor}.}
  \item{cex.pch}{numeric, the expansion factor for the points in \code{panel.Ct.scatter}.}
  \item{diag}{logical, should the diagonal line y=x be plotted.}
  \item{\dots}{any other arguments are passed to the panel function or \code{pairs}.}
}

\details{
Per default, the lower panels contain the correlations between data sets. For each correlation all complete pairs are used, i.e. NAs are ignored. If there are no complete observations between two samples the correlation will be set to NA.
}

\value{A figure is generated in the current graphics device.}

\author{Heidi Dvinge}

\seealso{\code{\link{pairs}} or \code{\link{plotCtScatter}} for plotting just two samples.
}

\examples{
# Load example data
data(qPCRraw)
# Various types of plot
plotCtPairs(qPCRraw, cards=1:4)
plotCtPairs(qPCRraw, col="black")
plotCtPairs(qPCRraw, Ct.max=40)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot }
