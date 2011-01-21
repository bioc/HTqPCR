\name{plotCtCard}
\Rdversion{1.1}
\alias{plotCtCard}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Image plot of qPCR Ct values from a card format}
\description{Function for plotting high-throughput qPCR Ct values from a platform with a defined spatial layout, such as TaqMan Low Density Assay cards. The location of Ct values in the plot corresponds to the position of each well on the card.}

\usage{
plotCtCard(q, card = 1, plot = "Ct", main, nrow = 16, ncol = 24, col, col.range, na.col = "grey", na.value = 40, legend.cols, well.size = 3.1, zero.center = FALSE, unR = FALSE, unD = FALSE, ...)
}

\arguments{
  \item{q}{object of class qPCRset. }
  \item{card}{integer, the sample number to plot.}
  \item{plot}{character string among "Ct", "flag", "type", "class") indicating what type of plot to produce. See Details for a longer description.}
  \item{main}{character string, the title of the plot. Per deault this is the sample name corresponding to card.}
  \item{nrow}{integer, the numer of rows on the card (16 for a standard 384 well format).}
  \item{ncol}{integer, the numer of columns on the card (24 for a standard 384 well format).}
  \item{col}{vector of colors of the same length as the number of different groups for the categorical data, or the name of a colour scheme for the continuous data.}
  \item{col.range}{vector, the range of colours to use.}
  \item{na.col}{the colour used for well with NA (undetermined) Ct values.}
  \item{na.value}{numeric, if NA has been replaced by an (arbitrary) high Ct value in the data.}
  \item{legend.cols}{integer, how many columns should the legend text be split into (defaults to number of labels).}
  \item{well.size}{numeric, for adjusting the size of the wells on the card.}
  \item{zero.center}{logical, should the colours be shifted to be zero-centered.}
  \item{unR}{logical, should wells from the category "Unreliable" be crossed out.}
  \item{unD}{logical, should wells from the category "Undetermined" be crossed out.}
  \item{\dots}{any other arguments will be passed to the \code{plot} and \code{points} functions.}
}

\details{This function may be used to plot the values of any well-specific information, such as the raw or normalized Ct values, or categorical data such as flag, gene class etc. The image follows the layout of an actual HTqPCR card. 
	
If \code{unR=TRUE} these will wells will be crossed out using a diagonal cross (X), whereas \code{unD=TRUE} will be marked with a horisontal/vertical cross.}

\value{A plot is created on the current graphics device.}

\author{Heidi Dvinge}

\seealso{\code{\link{image}}, and \code{\link{plotCtArray}} for plotting data from other high-throughput qPCR platforms (e.g. Fluidigm arrays).}

\examples{
# Load some example data
data(qPCRraw)
# Plot Ct values from first card
plotCtCard(qPCRraw)
plotCtCard(qPCRraw, card=2, col.range=c(10,35))
plotCtCard(qPCRraw, unR=TRUE, unD=TRUE)
# Other examples
plotCtCard(qPCRraw, plot="class")
plotCtCard(qPCRraw, plot="type")
plotCtCard(qPCRraw, plot="flag")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
