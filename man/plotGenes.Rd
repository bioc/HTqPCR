\name{plotCtLines}
\Rdversion{1.1}
\alias{plotCtLines}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting Ct values from qPCR across multiple samples.}
\description{This function is for displaying a set of features from a qPCRset across multiple samples, such as a timeseries or different treatments. Values for each feature are connected by lines, and the can be averaged across groups rather than shown for individual smaples.}

\usage{
plotCtLines(q, genes, groups, col = brewer.pal(10, "Spectral"), xlab = "Sample", ylab = "Ct", legend = TRUE, lwd = 2, lty, pch, xlim, ...)
}

\arguments{
 \item{q}{object of class qPCRset. }
 \item{genes}{numeric or character vector, selected genes to make the plot for.}
 \item{groups}{vector, the different groups that the samples in \code{q} belong to. See details.}
 \item{col}{vector, colours to use for the lines.}
 \item{xlab}{character string, label for the x-axis.}
 \item{ylab}{character string, label for the y-axis.}
 \item{legend}{logical, whether to include a colour legend or not.}
 \item{lwd}{numeric, the width of the lines.}
 \item{lty}{vector, line types to use. See \code{par} or \code{lines} for details.}
 \item{pch}{vector, if \code{groups} is set, the point types that will be used for each feature in \code{genes}.}
 \item{xlim}{vector of length two, the limits for the x-axis. Mainly used for adjusting the position of the legend.}
 \item{\dots}{any other arguments will be passed to the \code{matplot} function.}
}

\details{The default plot shows the Ct values across all samples in \code{q}, with lines connecting the samples. However, if \code{groups} is set the Ct values will be averaged within groups. Lines connect these averages, but the individual values are shown with different point types, as chosen in \code{pch}.
}

\value{A plot is created on the current graphics device.}

\author{Heidi Dvinge}

\seealso{\code{\link{matplot}} }

\examples{
# Load some example data
data(qPCRraw)
samples <- exFiles <- read.delim(file.path(system.file("exData", package="HTqPCR"), "files.txt"))
# Draw dfferent plots
plotCtLines(qPCRraw, genes=1:10)
plotCtLines(qPCRraw, genes=1:10, groups=samples$Treatment, xlim=c(0,3))
feat <- as.numeric(as.factor(featureType(qPCRraw)[1:10]))
plotCtLines(qPCRraw, genes=1:10, col=feat)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
