\name{cbind}
\alias{cbind}
\alias{rbind}
\alias{cbind.qPCRset}
\alias{rbind.qPCRset}
\title{Combine qPCRset objects}
\description{Functions for combining multiple \code{qPCRset} objects into one, by either adding columns (samples) or rows (features).}

\usage{
\method{cbind}{qPCRset}(..., deparse.level = 1)
\method{rbind}{qPCRset}(..., deparse.level = 1)
}

\arguments{
  \item{\dots}{\code{qPCRset} objects that are to be combined.}
  \item{deparse.level}{not implemented currently. See \code{\link[base]{cbind}}.}
}

\details{
In some cases it might be desirable to merge multiple \code{qPCRset} objects, that have been read into R or processed individually. This can be done for either identical samples across multiple different cards (such as a 384 well plate), or if more samples have been run on cards with the same layout. 
		
\code{cbind} combines data assuming that all experiments have been carried out on identical cards, i.e. that \code{featureNames}, \code{featureType}, \code{featurePos} and \code{featureClass} is identical across all the \code{qPCRset} objects. \code{rbind} combines data assuming that the same samples have been analysed using different qPCR cards.

For both functions, the \code{getCtHistory} of all the individual objects will be added to the combined \code{qPCRset}.
}

\value{A combined qPCRset object.}

\author{Heidi Dvinge}

\seealso{
\code{\link[base]{cbind}}}


\examples{
# Load some example data and split into multiple qPCRset objects
data(qPCRraw)
q1 <- qPCRraw[,1:2]
q2 <- qPCRraw[,2:4]
q3 <- qPCRraw[,5:6]
# Combine together by samples
q.samples <- cbind(q1,q3,q2)
n.wells(q.samples)
n.samples(q.samples)
# Combine as if the same samples had been run on multiple different cards
sampleNames(q3)	<- sampleNames(q1)
q.features <- rbind(q1,q3)
n.wells(q.features)
n.samples(q.features)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
