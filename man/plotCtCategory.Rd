\name{plotCtCategory}
\Rdversion{1.1}
\alias{plotCtCategory}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarising the feature categories for Ct values.}
\description{This function will provide a summary of the \code{featureCategory} for a qPCRset. Focus can either be on categories across samples, or across features.}

\usage{
plotCtCategory(q, cards = TRUE, by.feature = FALSE, stratify, col, xlim, main, ...)
}

\arguments{
  \item{q}{object of class qPCRset.}
  \item{cards}{integers, the number of the cards (samples) to plot.}
  \item{by.feature}{logical, should the categories be summarised for features rather than samples. See details.}
  \item{stratify}{character string, either "type" or "class" indicating if the categories should be stratified by \code{featureType} or \code{featureClass} of \code{q}. Ignored if \code{by.features} is TRUE.}
  \item{col}{vector with the colours to use for the categories. Default is green for "OK", yellow for "Unreliable" and red for "Undetermined".}
  \item{xlim}{vector, the limits of the x-axis. If \code{by.feature} is FALSE, this can be used to adjust the size of the barplot to fit in the colour legend.}
  \item{main}{character string, the title of the plot.}
  \item{\dots}{further arguments passed to \code{barplot} or \code{heatmap}.}
}

\details{
This function is for generating two different types of plot. If \code{by.feature=FALSE} the number of each \code{featureCategory} will be counted for each card, and a barplot is made. If however \code{by.feature=TRUE}, then the categories for each feature across the selected cards will be clustered in a heatmap.

For \code{by.feature=TRUE} the plot can be modified extensively using calls to the underlying \code{heatmap} function, such as setting \code{cexRow} to adjust the size of row labels.
}

\value{A figure is produced on the current graphics device.}

\author{Heidi Dvinge}

\seealso{
\code{\link{setCategory}}, and \code{\link{heatmap}} for the underlying plotting function.
}

\examples{
# Load example preprocessed data
data(qPCRpros)
# Plot categories for samples
plotCtCategory(qPCRpros)
plotCtCategory(qPCRpros, cards=1:3, stratify="class")
# Categories for features
plotCtCategory(qPCRpros, by.feature=TRUE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot }
