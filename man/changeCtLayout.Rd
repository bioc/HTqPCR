\name{changeCtLayout}
\alias{changeCtLayout}
\title{Changing the dimensions (rows x columns) of qPCRset objects}
\description{
A function for splitting up the individual qPCR cards, in case there are multiple samples present on each card. I.e. for cases where the layout isn't 1 sample x 384 features, but for example 4 samples x 96 features on each 384 well card.
}

\usage{
changeCtLayout(q, sample.order)
}

\arguments{
 \item{q}{a qPCRset object.}
 \item{sample.order}{vector, same length as number of features on each card (e.g. 384). See details.}
}

\details{
The result from each qPCR run of a given card typically gets presented together, such as in a file with 384 lines, one per feature, for 384 well plates. However, some cards may contain multiple samples, such as commercial cards that are designed to be loaded with two separate samples and then include 192 individual features.

Per default, each card is read into the \code{qPCRset} object as consisting of a single sample, and hence one column in the Ct data matrix. When this is not the case, the data can subsequently be split into the correct features x samples (rows x columns) dimensions using this function. The parameter \code{sample.order} is a vector, that for each feature in the \code{qPCRset} indicates what sample it actually belongs to.

In the new \code{qPCRset} the samples (Ct columns) are ordered first by \code{sample.order} then by the original \code{sampleNames}, as shown in the examples below.
}

\value{A qPCRset object like the input, but with the dimensions changed according to the new layout.}

\author{Heidi Dvinge}

\note{
Since the actual biological samples are likely to differ on each card, after applying \code{changeCtLayout} renaming of the samples in \code{qPCRset} using \code{sampleNames} is advisable.

The features are assumed to be identical for all samples on a given card! I.e. if for example \code{sample.order=rep(c("A", "B"), each=192)}, then feature number 1 (the first for sample A) should be the same as feature number 193 (the first for sample B). The new \code{featureNames} are taken for those features listed as belonging to the first sample in \code{sample.order}.
}

\examples{
# Example data
data(qPCRraw)
# With e.g. 2 or 4 samples per 384 well card.
sample2.order <- rep(c("subSampleA", "subSampleB"), each=192)
sample4.order <- rep(c("subA", "subB", "subC", "subD"), each=96)
# Splitting the data into all individual samples
qPCRnew2 <- changeCtLayout(qPCRraw, sample.order=sample2.order)
show(qPCRnew2)
qPCRnew4 <- changeCtLayout(qPCRraw, sample.order=sample4.order)
show(qPCRnew4)
sampleNames(qPCRnew4)
}

\keyword{ htest }
