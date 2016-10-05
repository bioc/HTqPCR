\name{mannwhitneyCtData}
\alias{mannwhitneyCtData}

\title{Differentially expressed features with qPCR: Mann-Whitney}
\description{Function for calculating p-values across two groups for the features present in high-throughput qPCR data, such as from TaqMan Low Density Arrays. Also known as two sample Wilcoxon test.}

\usage{
mannwhitneyCtData(q, groups = NULL, calibrator, alternative = "two.sided", paired = FALSE, replicates = TRUE, sort = TRUE, stringent = TRUE, p.adjust = "BH", ...)
}

\arguments{
 \item{q}{qPCRset object.}
 \item{groups}{factor, assigning each sample to one of two groups.}
 \item{calibrator}{which of the two groups is to be considered as the reference and not the test? Defaults to the first group in \code{groups}.}
 \item{alternative}{character string (first letter is enough), specifying the alternative hypothesis, "two.sided" (default), "greater" or "less".}
 \item{paired}{logical, should a paired t-test be used.}
 \item{replicates}{logical, if replicated genes are present on the array, the statistics will be calculated for all the replicates combined, rather than the individual wells.}
 \item{sort}{boolean, should the output be sorted by p-values.}
 \item{stringent}{boolean, for flagging results as "Undetermined". See details.}
 \item{p.adjust}{character string, which method to use for p-value adjustment for multiple testing. See details.}
 \item{\dots}{any other arguments will be passed to the \code{\link{wilcox.test}} function.}
}

\details{Once the Ct values have been normalised, differential expression can be calculated. This function deals with just the simple case, where there are two types of samples to compare. For a parametric test see \code{ttestCtData} and \code{limmaCtData} for more complex studies.
	
The underlying statistics is calculated by \code{wilcox.test}. Due to the high possibility of ties for each feature between samples, the test is run with \code{exact=FALSE}.

All results are assigned to a category, either "OK" or "Undetermined" depending on the input Ct values. If \code{stringent=TRUE} any unreliable or undetermined measurements among technical and biological replicates will result in the final result being "Undetermined". For \code{stringent=FALSE} the result will be "OK" unless at least half of the Ct values for a given gene are unreliable/undetermined.

The argument \code{p.adjust} is passed on to the \code{\link{p.adjust}} function. Options include e.g. "BH" (Benjamini & Hochberg, the default), "fdr" and "bonferroni". See \code{\link{p.adjust}} for more information on the individual methods.
}


\value{A data.frame containing the following information:
	\item{genes}{The names of the features on the card.}
	\item{feature.pos}{The \code{featurePos} of the genes. If replicated genes are used, the feature positions will be concatenated together.}
	\item{MB.test}{The name and value of the test statistic.}
	\item{p.value}{The corresponding p-value.}
	\item{ddCt}{The delta delta Ct values.}
	\item{FC}{The fold change; 2^(-ddCt).}
	\item{meanCalibrator}{The average expression level of each gene in the calibrator sample(s).}
	\item{meanTarget}{The average expression level of each gene in the target sample(s).}
	\item{categoryCalibrator}{The category of the Ct values ("OK", "Undetermined") across the calibrator.}
	\item{categoryTarget}{Ditto for the target.}
}

\author{Heidi Dvinge}

\seealso{\code{\link{wilcox.test}}, \code{\link{ttestCtData}}, \code{\link{limmaCtData}}. \code{\link{plotCtRQ}} and \code{\link{plotCtSignificance}} can be used for visualising the results.}


\examples{
# Load example preprocessed data
data(qPCRpros)
# Test between two groups, collapsing replicated features
diff.exp <- mannwhitneyCtData(qPCRpros[,1:4], groups=factor(c("A", "B", "B", "A")), calibrator="B")
diff.exp[1:10,]
# The same test, taking replicated features individually
diff.exp <- mannwhitneyCtData(qPCRpros[,1:4], groups=factor(c("A", "B", "B", "A")), calibrator="B", replicates=FALSE)
# Using another method for p-value adjustment
diff.exp <- mannwhitneyCtData(qPCRpros[,1:4], groups=factor(c("A", "B", "B", "A")), calibrator="B", p.adjust="holm")
}

\keyword{htest }
