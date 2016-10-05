\name{plotCtReps}
\Rdversion{1.1}
\alias{plotCtReps}
\title{Scatter plot of features analysed twice during each qPCR experiment.}
\description{In high-throughput qPCR data some features may be present twice on each card (sample). This function will make a scatter plot of one replicate versus the other for each sample individually, as well as mark genes with very deviating replicate values.
}

\usage{
plotCtReps(q, card = 1, percent = 20, verbose = TRUE, col = 1, ...)
}

\arguments{
 \item{q}{object of class qPCRset.}
 \item{card}{integer, the sample number to plot.}
 \item{percent}{numeric, features with replicate values differ more than this percentage from their average will be marked on the plot.}
 \item{verbose}{logical, should the deviating genes and their Ct values be printed to the terminal.}
 \item{col}{integer or character; the colour of the points in the scatter plot.}
 \item{\dots}{any other arguments are passed to \code{\link{plot}}.}
}

\details{This function will look through the data in the qPCRset, find all genes with are presented twice on the array, and plot the Ct values of these replicated genes against each other. Whether a genes goes to the x or y-axis depends on the first occurrence of the gene names.
	
All genes where abs(rep1-rep2) > percent/100*replicate mean will be marked by an open circle, and the gene names written in red letters.}

\value{
An plot is created on the current graphics device. Also, a data.frame with the names and values of deviating genes is returned invisibly.}

\author{Heidi Dvinge}

\seealso{\code{\link{plot}}, and \code{\link{par}} for the plotting parameters.
}

\examples{
# Load example data
data(qPCRraw)
# Plot replicates
plotCtReps(qPCRraw, card=1, percent=30)
plotCtReps(qPCRraw, card=2, percent=10)
reps <- plotCtReps(qPCRraw, card=2, percent=20)
reps
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
