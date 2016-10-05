\name{filterCtData}
\Rdversion{1.1}
\alias{filterCtData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Filter out features (genes) from qPCR data.}

\description{This function is for filtering Ct data from high-throughput qPCR platforms like the TaqMan Low Density Arrays. This can for example be done prior to analysing the statistical significance of the data, to remove genes where the results are of low quality, or that are not of interest to the analysis in question.}

\usage{
filterCtData(q, remove.type, remove.name, remove.class, remove.category, n.category = 3, remove.IQR, verbose = TRUE)
}

\arguments{
 \item{q}{object of class qPCRset.}
 \item{remove.type}{character vector, the feature type(s) to be removed from the data object.}
 \item{remove.name}{character vector, the feature name(s) to be removed from the data object.}
 \item{remove.class}{character vector, the feature class(es) to be removed from the data object.}
 \item{remove.category}{character vector, the features categories(s) to be assessed across samples.}
 \item{n.category}{numeric, all features with more than this number of \code{remove.category} across samples are removed.}
 \item{remove.IQR}{numeric, all features with an interquartile range (IQR) below this limit across samples will be removed.}
 \item{verbose}{boolean, should some information be printed to the prompt.}
}

\details{This function may be used to exclude individual or small groups of features that are irrelevant to a given analysis. However, it can also be used on a more general basis, to for example split the data into separate \code{qPCRset} objects based on features with different characteristics, such as groups of markers or other gene classes present in \code{featureClass}.
	
\code{remove.IQR} can be used to exclude features that show only little variation across the samples. These are unlikely to be differentially expressed, so including them in downstream analysis such as \code{limmaCtData} or \code{ttestCtData} would result in a slight loss of power caused by the adjustment of p-values required due to multiple testing across all features.}

\value{An object of class qPCRset like the input, but with the required features removed.}

\author{Heidi Dvinge}

\note{After removing features the function \code{\link{plotCtCard}} will no longer work, since the number of features is now smaller than the card dimensions.

When using \code{remove.category} or \code{remove.IQR} and there are replicated features present on the array, it might no longer be possible to use the \code{ndups} parameter of \code{\link{limmaCtData}}, since the number of replicates isn't identical for each feature.

Filtering can be performed either before or after normalization, but in some cases normalization might be affected by this, for example if many features are removed, making it difficult to identify rank-invariant genes.
}

\examples{
# Load some example data
data(qPCRpros)
show(qPCRpros)
# Filter based on different feature type
qFilt <- filterCtData(qPCRpros, remove.type=c("Endogenous Control"))
# Filter based on feature type and name
qFilt <- filterCtData(qPCRpros, remove.type=c("Endogenous Control"), remove.name=c("Gene1", "Gene20", "Gene30"))
# Filter based on feature class
qFilt <- filterCtData(qPCRpros, remove.class="Kinase")
# Filter based on feature categories, using two different cut-offs
qFilt <- filterCtData(qPCRpros, remove.category="Undetermined")
qFilt <- filterCtData(qPCRpros, remove.category="Undetermined", n.category=5)
# Remove features without much variation across samples
iqr <- apply(exprs(qPCRpros), 1, IQR, na.rm=TRUE)
hist(iqr, n=20)
qFilt <- filterCtData(qPCRpros, remove.IQR=2)
}

\keyword{ htest }
