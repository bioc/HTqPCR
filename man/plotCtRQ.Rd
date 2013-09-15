\name{plotCtRQ}
\Rdversion{1.1}
\alias{plotCtRQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot the relative quantification of Ct values from qPCR experiments.}
\description{Function for plotting the relative quantification (RQ) between two groups of data, whose Ct values have been tested for significant differential expression.}

\usage{
plotCtRQ(qDE, comparison = 1, genes, transform = "log2", p.val = 0.1, mark.sig = TRUE, p.sig = 0.05, p.very.sig = 0.01, mark.un = TRUE, un.tar = "black", un.cal = "black", col, legend = TRUE, xlim, mar, main, ...)
}

\arguments{
  \item{qDE}{list or data.frame, the result from \code{ttestCtData} or \code{limmaCtData}.}
  \item{comparison}{integer or character string, indicating which component to use if \code{qDE} is a list.}
  \item{genes}{numeric or character vector, selected genes to make the plot for.}
  \item{transform}{character string, how should the data be displayed. Options are "none", "log2" or "log10". See details}
  \item{p.val}{numeric between 0 and 1, if \code{genes} is not supplied all given with (adjusted) p-value below this threshold will be included.}
  \item{mark.sig}{logical, should significant features be marked.}
  \item{p.sig}{numeric, the cut-off for significant p-values that will be marked by *.}
  \item{p.very.sig}{numeric, the cut-off for very significant p-values that will be marked by ".}
  \item{mark.un}{logical, should data with unreliable target or calibrator samples be marked. See details.}
  \item{un.tar}{colour to use for the undetermined targets. See details.}
  \item{un.cal}{colour to use for the undetermined calibrators. See details.}
  \item{col}{vector, colours to use for the bars.}
  \item{legend}{logical, should a legend be included in the barplot.}
  \item{xlim}{vector of length 2, the limits on the x-axis. Mainly used for moving the legend to the left of bars.}
  \item{mar}{vector with 4 values, the size of the margins. See \code{par} for more info.}
  \item{main}{character string, the image title. Default to the name of the chosen comparison.}
  \item{\dots}{any other arguments will be passed to the \code{barplot} function.}
}

\details{
The relative quantification is calculated as RQ=2^-ddCT, where ddCT is the deltadeltaCt value.	

If \code{mark.un=TRUE}, those bars where either the calibrator or target sample measurements were undetermined are marked using diagonal lines. Whether either of these are called undetermined (includes unreliable values) or not depends on all the input Ct values in \code{ttestCtData} or \code{limmaCtData}, and whether \code{stringent=TRUE} was used in these functions.
}

\value{A plot is created on the current graphics device.}

\author{Heidi Dvinge}

\seealso{
\code{\link{ttestCtData}} and \code{\link{limmaCtData}} for testing the Ct data for differential expression. 
}

\examples{
# Load example data and calculate differential expression
data(qPCRpros)
qDE <- ttestCtData(qPCRpros[,1:4], groups=factor(c("A", "B", "B", "A")), calibrator="B")
# Plotting the top 10 results or first 10 genes
plotCtRQ(qDE, genes=1:10)
plotCtRQ(qDE, genes=featureNames(qPCRpros)[1:10])
# Plot all results with p-value below 0.08
plotCtRQ(qDE, p.val=0.08, transform="none")
plotCtRQ(qDE, p.val=0.08, transform="log10")

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot }
