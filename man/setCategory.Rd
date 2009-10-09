\name{setCategory}
\Rdversion{1.1}
\alias{setCategory}
\title{Assign categories to Ct values from qPCR data.}
\description{Data in qPCRset objects will have feature categories ("Unreliable", "Undtermined") assigned to them based on different Ct criteria. }

\usage{
setCategory(q, Ct.max = 35, Ct.min = 10, replicates = TRUE, quantile = 0.9, groups, flag = TRUE, flag.out = "Failed", verbose = TRUE, plot = FALSE, ...)
}

\arguments{
  \item{q}{qPCRset object.}
  \item{Ct.max}{numeric, the maximum tolerated Ct value. Everything above this will be "Undetermined".}
  \item{Ct.min}{numeric, the minimum tolerated Ct value. Everything below this will be "Unreliable".}
  \item{replicates}{logical, should Ct values from genes replicated within each sample be collapsed for the standard deviation.}
  \item{quantile}{numeric from 0 to 1, the quantile interval accepted for standard deviations. \code{NULL} means that variation between replicates is not used for setting the categories. See details.}
  \item{groups}{vector, grouping of cards, for example biological or technical replicates}
  \item{flag}{logical, should categories also be set to "Unreliable" according to the content of \code{flag(q)}.}
  \item{flag.out}{character vector, if \code{flag=TRUE}, what are the flag(s) to be set as "Unreliable".}
  \item{verbose}{logical, should a summary about category counts per sample be printed to the prompt.}
  \item{plot}{logical, should some plots of the standard deviations be created.}
  \item{\dots}{any other arguments are passed to \code{plot}.}
}

\details{Categories can be assigned to the \code{featureCategory} of the qPCRset using either just simple criteria (max/min of Ct values or \code{flag} of \code{q}) or by looking at the standard deviation of Ct values across biological and technical replicates for each gene.

When looking at replicates, the standard deviation and mean are calculated and a normal distribution following these parameters is generated. Individual Ct values that are outside the interval set by \code{quantile} are set as "Unreliable". %This interval means whether they eviate from the e.g. 90\% quantile estimeated from a normal distribution of the given mena and standard deviation.
So if e.g. \code{quantile=90} the values outside the top 5\% and lower 5\% of the normal distribution with the given mean and standard deviation are removed. 

"Undetermined" has priority over "Unreliable", so if a value is outside \code{quantile} but also above \code{Ct.max} it will be "Undetermined".
	
NB: When setting categories based on replicates, the Ct values are assumed to follow a normal distribution. This might not be the case if the number of samples within each group is small, and there are no replicates on the genes within each sample.

If the number of replicates vary significantly between biological groups, this will influence the thresholds used for determining the range of "OK" Ct values.}

\value{If \code{plot=TRUE} one figure per sample grous is returned to the current graphics device.
A qPCRset with the new feature categories is returned invisibly.}

%\references{
%% ~put references to the literature/web site here ~
%}

\note{
It's adviced to try several different values for \code{quantile}, depending on the input data set. Using the function \code{PlotCtCategory(..., by.feature=FALSE)} or \code{plotCtCategory(..., by.feature=TRUE)} might help assess the result of different \code{quantile} choices.
}

\author{Heidi Dvinge}

\seealso{
\code{\link{filterCategory}}, \code{\link{plotCtCategory}}
}

\examples{
# Load example data
data(qPCRraw)
exFiles <- read.delim(file.path(system.file("exData", package="HTqPCR"), "files.txt"))
 # Set categories in various ways
setCategory(qPCRraw, flag=FALSE, quantile=NULL)
setCategory(qPCRraw[,1:4], groups=exFiles$Treatment[1:4], plot=TRUE)
setCategory(qPCRraw[,1:4], groups=exFiles$Treatment[1:4], plot=TRUE, quantile=0.80)
x <- setCategory(qPCRraw, groups=exFiles$Treatment, verbose=FALSE, quantile=0.80)
# Plot the categories
plotCtCategory(x)
}

\keyword{htest }
