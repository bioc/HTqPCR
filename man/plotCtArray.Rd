\name{plotCtArray}
\alias{plotCtArray}

\title{Image plot of qPCR Ct values from an array format}
\description{Function for plotting high-throughput qPCR Ct values from a platform with a defined spatial layout, such as Fluidigm Dynamic Arrays (BioMark) or OpenArray from Applied Biosystems. The location of Ct values in the plot corresponds to the position of each well on the array.}

\usage{
plotCtArray(q, plot = "Ct", main, col, col.range, na.col = "grey", na.value = 40, chamber.size, ...)
}

\arguments{
  \item{q}{object of class qPCRset. }
  \item{plot}{character string indicating what type of plot to produce. Currently only "Ct" is implemented.}
  \item{main}{character string, the title of the plot. Per default "Ct values".}
  \item{col}{the name of a colour scheme.}
  \item{col.range}{vector, the range of colours to use.}
  \item{na.col}{the colour used for well with NA (undetermined) Ct values.}
  \item{na.value}{numeric, if NA has been replaced by an (arbitrary) high Ct value in the data.}
  \item{chamber.size}{numeric, for adjusting the size of the reaction chamber on the card.}
  \item{\dots}{any other arguments will be passed to the \code{plot} function.}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{A plot is created on the current graphics device.}

\author{Heidi Dvinge}

\seealso{\code{\link{plotCtCard}} for plotting data from other high-throughput qPCR platforms.}

\examples{
	
# Locate example data
exPath <- system.file("exData", package="HTqPCR")
exFiles <- "BioMark_sample.csv"
# Create qPCRset object
raw	<- readCtData(exFiles, path=exPath, n.features=48, n.data=48, format="BioMark")
# Plot
plotCtArray(raw)
# Change colour and range
plotCtArray(raw, col=brewer.pal(11, "Spectral"), col.range=c(10,35))

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
