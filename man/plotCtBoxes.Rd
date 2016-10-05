\name{plotCtBoxes}
\Rdversion{1.1}
\alias{plotCtBoxes}
\title{Boxplots for qPCR Ct values.}
\description{Function for making boxplots of Ct values from high-throughput qPCR data. The boxes can be made either using all values on each card, or stratified by different feature information.}

\usage{
plotCtBoxes(q, cards = TRUE, xlab = "", col, main = NULL, names, stratify = "type", mar = c(7, 4, 3, 1), ...)
}

\arguments{
 \item{q}{object of class qPCRset. }
 \item{cards}{vector, the numbers of the cards to plot. Defaults to TRUE = all cards.}
 \item{xlab}{character string, label for the x-axis.}
 \item{col}{vector of colours to use, defaults to different colour for each card.}
 \item{main}{character string, plot title.}
 \item{names}{vector, names to plot under the boxes. Defaults to sample names.}
 \item{stratify}{character, specifying what to stratify the Ct values by. NULL, the default means no stratification, "type" is the feature types of the qPCRset, and "class" the feature class.}
 \item{mar}{vector, the size of the margins. See \code{\link{par}} for details.}
 \item{\dots}{any other arguments will be passed to the \code{boxplot} or \code{par}function.}
}

\details{
% If \code{q} is a qPCRset it can contain either the raw or normalized Ct values. deltaCt or deltadeltaCt (abbreviated to dCt and ddCt) can be calculated and put into the required list format using the \code{comparativeCt} function.
	
For the stratified plots all boxes with Ct values from the same card are plotted in identical colours. "type" and "class" are automatically extracted from the qPCRset using \code{featureType} and \code{featureClass}.}

\value{A plot is created on the current graphics device.}

\author{Heidi Dvinge}

\seealso{\code{\link{boxplot}}}

\examples{
# Loading the data
data(qPCRraw)
# Make plot with all samples or just a few
plotCtBoxes(qPCRraw, stratify=NULL)
plotCtBoxes(qPCRraw, cards=c(1,4))
plotCtBoxes(qPCRraw, stratify="class")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
