\name{plotCtSignificance}
\Rdversion{1.1}
\alias{plotCtSignificance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Barplot with Ct values between genes from qPCR.}
\description{Function for producing a barplot of the Ct values from high-throughput qPCR samples. A comparison is made between two groups which have been tested for differential expression, and all individual Ct values are shown, to identify potential outliers. }

\usage{
plotCtSignificance(qDE, q, comparison = 1, genes, p.val = 0.1, groups, calibrator, target, p.sig = 0.05, p.very.sig = 0.01, mark.sig = TRUE, col, un.col = "#D53E4F", point.col = "grey", legend = TRUE, mar, main, jitter = 0.5, ...)
}

\arguments{
 \item{qDE}{list or data.frame, the result from \code{ttestCtData} or \code{limmaCtData}.}
 \item{q}{the qPCRset data that was used for testing for differential expression.}
 \item{comparison}{integer or character string, indicating which component to use if \code{x} is a list.}
 \item{genes}{numeric or character vector, selected genes to make the plot for.}
 \item{p.val}{numeric between 0 and 1, if \code{genes} is not supplied all given with (adjusted) p-value below this threshold will be included.}
 \item{groups}{vector, the groups of all the samples in \code{q}.}
 \item{calibrator}{character string, which of the \code{groups} is the calibrator.}
 \item{target}{character string, which of the \code{groups} is the target.}
 \item{p.sig}{numeric, the cut-off for significant p-values that will be marked by *.}
 \item{p.very.sig}{numeric, the cut-off for very significant p-values that will be marked by ".}
 \item{mark.sig}{logical, should significant features be marked.}
 \item{col}{vector, colours to use for the two sets of bars, one per sample type.}
 \item{un.col}{integer or character string, the colour to use for all Ct values that are "Unreliable" or "Undetermined".}
 \item{point.col}{integer or character string, the colour to use for all other Ct values.}
 \item{legend}{logical, should a legend be included int eh barplot.}
 \item{mar}{vector with 4 values, the size of the margins. See \code{par} for more info.}
 \item{main}{character string, the image title. Default to the name of the chosen comparison.}
 \item{jitter}{numeric, between 0 and 1. If Ct values are very similar, the individual points might lie on top of each other in the bars. This adds a jittering factor along the x-axis. If 0 the points will all be aligned.}
 \item{\dots}{any other arguments will be passed to the \code{barplot} function.}
}

\details{This function will make a barplot with the average Ct values for the test and reference samples for the selected genes. All the individual Ct values are plotted on top of the bars though, and the "Unreliable" or "Undetermined" ones are marked, to do a visual assessment of the impact of non-valid measurements on the average.
	
It's up to the user to specify the correct \code{calibrator} and \code{target} for the given comparison; no checking is done.	
}

\value{A plot is created on the current graphics device.}

\author{Heidi Dvinge}

\seealso{\code{\link{barplot}} and \code{\link{plotCtRQ}} or \code{\link{plotCtOverview}} for a plot of the relative quantification between samples.
}

\examples{
# Load example data and calculate differential expression
data(qPCRpros)
grp <- factor(c("A", "B", "B", "A"))
qDE <- ttestCtData(qPCRpros[,1:4], groups=grp, calibrator="B")
# Plot
plotCtSignificance(qDE, q=qPCRpros, groups=grp, target="A", calibrator="B", genes=1:10, jitter=0.2)
plotCtSignificance(qDE, q=qPCRpros, groups=grp, target="A", calibrator="B", genes=featureNames(qPCRpros)[1:10], jitter=0, un.col="blue")
plotCtSignificance(qDE, q=qPCRpros, groups=grp, target="A", calibrator="B", p.val=0.05, jitter=0.1)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
