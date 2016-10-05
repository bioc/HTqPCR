\name{plotCVBoxes}
\Rdversion{1.1}
\alias{plotCVBoxes}
\title{Boxplots of CV for qPCR Ct values.}
\description{Function that will calculate the coefficients of variation across selected qPCR data, and plot the results in a boxplot.}

\usage{
plotCVBoxes(q, cards = TRUE, xlab = "", ylab = "CV", col = brewer.pal(5, "Spectral"), main = NULL, stratify, ...)
}

\arguments{
 \item{q}{object of class qPCRset.}
 \item{cards}{vector, the numbers of the cards to plot. Defaults to TRUE = all cards.}
 \item{xlab}{character string, label for the x-axis.}
 \item{ylab}{character string, label for the y-axis.}
 \item{col}{vector of colours to use.}
 \item{main}{character string, plot title.}
 \item{stratify}{character, specifying what to stratify the Ct values by. NULL, the default means no stratification, "type" is the feature types of the qPCRset, and "class" the feature class.}
 \item{\dots}{any other arguments will be passed to the \code{boxplot} function.}
}

\details{The CV is calculated across all the selected cards based on each well position, without taking possibly replicated genes on the cards into consideration. "type" and "class" are automatically extracted from the qPCRset using \code{featureType} and \code{featureClass}.}

\value{A plot is created on the current graphics device. The CV values are returned invisibly.}

\author{Heidi Dvinge}

\seealso{\code{\link{boxplot}}}

\examples{
# Load example data
data(qPCRraw)
# Make plot with all samples or just a few
plotCVBoxes(qPCRraw)
plotCVBoxes(qPCRraw, cards=c(1,4))
plotCVBoxes(qPCRraw, stratify="class")
x <- plotCVBoxes(qPCRraw, stratify="type")
x[1:10]
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
