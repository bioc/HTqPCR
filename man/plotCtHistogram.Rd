\name{plotCtHistogram}
\Rdversion{1.1}
\alias{plotCtHistogram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Histrogram of Ct values from qPCR experiments.}
\description{The distribution of Ct values for a selected qPCR sample is shown in a histogram.}

\usage{
plotCtHistogram(q, card = 1, xlab = "Ct", col, main, n = 30, ...)
}

\arguments{
 \item{q}{an object of class qPCRset.}
 \item{card}{integer, the number of the card (sample) to plot.}
 \item{xlab}{character string, the label for the x-axis.}
 \item{col}{integer or character, the colour for the histogram.}
 \item{main}{character string, the plot title. Default is the name of the sample.}
 \item{n}{integer, number of bins to divide the Ct values into.}
 \item{\dots}{any other arguments are passed to \code{hist}.}
}

\value{A figure is generated in the current graphics device.}

\author{Heidi Dvinge}

\seealso{
\code{\link{plotCtDensity}} or \code{\link{plotCtBoxes}} for including multiple samples in the same plot.
}

\examples{
# Load example data
data(qPCRraw)
# Create the plots
plotCtHistogram(qPCRraw, card=2)
plotCtHistogram(qPCRraw, card=3, n=50, col="blue")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot }
