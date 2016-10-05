\name{limmaCtData}
\Rdversion{1.1}
\alias{limmaCtData}
\title{Differentially expressed features with qPCR: limma}
\description{Function for detecting differentially expressed genes from high-throughput qPCR Ct values, based on the framework from the \code{limma} package. Multiple comparisons can be performed, and across more than two groups of samples.}

\usage{
limmaCtData(q, design = NULL, contrasts, sort = TRUE, stringent = TRUE, ndups = 1, spacing = NULL, dupcor, ...)
}

\arguments{
 \item{q}{object of class qPCRset.}
 \item{design}{matrix, design of the experiment rows corresponding to cards and columns to coefficients to be estimated. See details.}
 \item{contrasts}{matrix, with columns containing contrasts. See details}
 \item{sort}{boolean, should the output be sorted by adjusted p-values.}
 \item{stringent}{boolean, for flagging results as "Undetermined". See details.}
 \item{ndups}{integer, the number of times each feature is present on the card.}
 \item{spacing}{integer, the spacing between duplicate spots, spacing=1 for consecutive spots}
 \item{dupcor}{list, the output from \code{duplicateCorrelation}. See details.}
 \item{\dots}{any other arguments are passed to \code{lmFit}, \code{contrasts.fit}, \code{eBayes} or \code{decideTests}.}
}

\details{This function is a wrapper for the functions \code{lmFit}, \code{contrasts.fit} (if a contrast matrix is supplied) and \code{eBayes} from the limma package. See the help pages for these functions for more information about setting up the design and contrast matrices.

All results are assigned to a category, either "OK" or "Unreliable" depending on the input Ct values. If \code{stringent=TRUE} any unreliable or undetermined measurements among technical and biological replicates will result in the final result being "Undetermined". For \code{stringent=FALSE} the result will be "OK" unless at least half of the Ct values for a given gene are unreliable/undetermined.

Note that when there are replicated features in the samples, each feature is assumed to be present the same number of times, and with regular spacing between replicates. Reordering the sample by \code{featureNames} and setting \code{spacing=1} is recommendable.

If technical sample replicates are available, \code{dupcor} can be used. It is a list containing the estimated correlation between replicates. \code{limmaCtData} will then take this correlation into account when fitting a model for each gene. It can be calculate using the function \code{\link[limma]{duplicateCorrelation}}. Technical replicates and duplicated spots can't be assessed at the same time though, so if \code{dupcor} is used, \code{ndups} should be 1.
}

\value{A list of data.frames, one for each column in \code{design}, or for each comparison in \code{contrasts} if this matrix is supplied. Each component of the list contains the result of the given comparisons, with one row per gene and has the columns:
	\item{genes }{Feature IDs.}
	\item{feature.pos}{The unique feature IDs from \code{featurePos} of the \code{q} object. Useful if replicates are not collapsed, in which case there might be several features with identical names.}
	\item{t.test }{The result of the t-test.}
	\item{p.value }{The corresponding p.values.}
	\item{adj.p.value }{P-values after correcting for multiple testing using the Benjamini-Holm method.}
	\item{ddCt }{The deltadeltaCt values.}
	\item{FC}{The fold change; 2^(-ddCt).}
	\item{meanTest }{The average Ct across the test samples for the given comparison.}
	\item{meanReference }{The average Ct across the reference samples for the given comparison.}
	\item{categoryTest }{The category of the Ct values ("OK", "Undetermined") across the test samples for the given comparison.}
	\item{categoryReference }{The category of the Ct values ("OK", "Undetermined") across the reference samples for the given comparison.}
Also, the last item in the list is called "Summary", and it's the result of calling \code{decideTests} from limma on the fitted data. This is a data frame with one row per feature and one column per comparison, with down-regulation, no change and up-regulation marked by -1, 0 and 1.
}

\references{
 Smyth, G. K. (2005). Limma: linear models for microarray data. In: Bioinformatics and Computational Biology Solutions using R and Bioconductor.
R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber (eds), Springer, New York, pages 397--420.
}

\author{Heidi Dvinge}


\seealso{
\code{\link[limma]{lmFit}}, \code{\link[limma]{contrasts.fit}} and \code{\link[limma]{ebayes}} for more information about the underlying limma functions. \code{\link{mannwhitneyCtData}} and \code{\link{ttestCtData}} for other functions calculating differential expression of Ct data. \code{\link{plotCtRQ}}, \code{\link{heatmapSig}} and \code{\link{plotCtSignificance}} can be used for visualising the results.
}

\examples{
# Load example preprocessed data
data(qPCRpros)
samples <- read.delim(file.path(system.file("exData", package="HTqPCR"), "files.txt"))
# Define design and contrasts
design <- model.matrix(~0+samples$Treatment)
colnames(design) <- c("Control", "LongStarve", "Starve")
contrasts <- makeContrasts(LongStarve-Control, LongStarve-Starve, Starve-Control, levels=design)
# The actual test
diff.exp <- limmaCtData(qPCRpros, design=design, contrasts=contrasts)
# Some of the results
diff.exp[["LongStarve - Control"]][1:10,]
# Example with duplicate genes on card.
# Reorder data to get the genes in consecutive rows
temp <- qPCRpros[order(featureNames(qPCRpros)),]
diff.exp <- limmaCtData(temp, design=design, contrasts=contrasts, ndups=2, spacing=1)
# Some of the results
names(diff.exp)
diff.exp[["LongStarve - Control"]][1:10,]
diff.exp[["Summary"]][1:10,]
}

\keyword{ htest }
