\name{plotCtScatter}
\Rdversion{1.1}
\alias{plotCtScatter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Scatterplot of two sets of Ct values from qPCR data.}

\description{Produces a plot of Ct values from two samples plotted against each other. Features can be marked based on for example feature class or type.}

\usage{
plotCtScatter(q, cards = c(1, 2), col = "class", pch = 20, diag = FALSE, cor = TRUE, Ct.max = 35, legend = TRUE, ...)
}

\arguments{
  \item{q}{object of class qPCRset. }
  \item{cards}{vector, the two cards to plot against each other.}
  \item{col}{vector with the colour(s) to use for the points, or a character string ("type" or "class") indicating whether points should be coloured according to \code{featureType} or \code{featureClass} of \code{q}.}
  \item{pch}{integer, the point type to use for the plot.}
  \item{diag}{logical, should the diagonal line y=x be plotted.}
  \item{cor}{logical, should information about the correlation between the two samples be included in the plot. The correlation is calculated both with and without removing Ct values above Ct.max.}
  \item{Ct.max}{numeric, all Ct values above this will be removed for calculating one of the correlations.}
  \item{legend}{logical, if \code{col} is either "type" or "class", should a colour legend for these be included.}
  \item{\dots}{any other arguments are passed to \code{plot}.}
}

%\details{}

\value{A figure is generated in the current graphics device.}

\author{Heidi Dvinge}

%\seealso{ }

\examples{
# Load example data
data(qPCRraw)
# Various types of plot
plotCtScatter(qPCRraw, cards=c(1,2))
plotCtScatter(qPCRraw, cards=c(1,4), col="type")
plotCtScatter(qPCRraw, cards=c(1,4), col="black", cor=FALSE, diag=TRUE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot }
