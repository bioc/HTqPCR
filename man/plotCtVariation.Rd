\name{plotCtVariation}
\alias{plotCtVariation}
\title{Plot variation in Ct values across replicates}
\description{Examine the variation in Ct values, either across features present multiple times on each card, or for within different groups of samples. The function supports both a summarised and a more detailed output.}

\usage{
plotCtVariation(q, cards = TRUE, variation = "var", type = "summary", sample.reps, feature.reps, log = FALSE, add.featurenames = FALSE, ylab, n.col, ...)
}

\arguments{
 \item{q}{object of class qPCRset. }
 \item{cards}{vector, the numbers of the cards to plot. Defaults to TRUE = all cards.}
 \item{variation}{character string indication whether to calculate the variation, "var", or standard deviation, "sd".}
 \item{type}{character string indicating whether to output the results in a summarised boxplot, "summary" or as a more detailed scatter plot, "detail". See Details and the examples.}
 \item{sample.reps}{a vector grouping the samples (see Details). Overrides \code{feature.reps}.}
 \item{feature.reps}{a vector grouping the features according to which are replicates. Per default \code{featureNames(q)} are used.}
 \item{log}{logical, should the results be converted into log10 values.}
 \item{add.featurenames}{logical, if \code{type="detail"} should the names of each feature be added to the scatter plot.}
 \item{ylab}{character, the label of the y-axis.}
 \item{n.col}{integer, if \code{type="detail"} how many columns should the scatterplots be presented in. Defaults to 3, or \code{n.samples(q)} if <3.}
 \item{\dots}{further arguments passed to \code{boxplot} or \code{plot}.}
}

\details{It is often useful to examine the data to determine if some samples are inherently more variable than other, or if the concordance between replicates on each qPCR card is acceptable. Using \code{type="summary"} generates a boxplot with all the variation values, either across genes (if \code{sample.reps} is set) or with each samples (default, or if \code{feature.reps} is set). That way the general distribution of variation or standard deviation values can be compared quickly.
	
If it looks like there's an unacceptable (or interesting) difference in the variation, this can be further investigated using \code{type="detail"}. This will generate multiple sub-plots, containing a single scatterplot of variation versus mean for each gene (if \code{sample.reps} is set) or each sample (default, or if \code{feature.reps} is set). Including the mean in the plot can be used to assess heteroskedasticity in the data. }

\value{A plot is created on the current graphics device. The variation and mean across each type of replicate is returned invisibly in a list with "Var" and "Mean" slots.}

\author{Heidi Dvinge}

\seealso{
\code{\link{plotCtReps}} for cases where the qPCR card only contains two replicates of each feature. \code{\link{plotCVBoxes}} for other ways of plotting variation within different groups.
}

\examples{

# Load some example data
data(qPCRraw)

# Get a summary of the standard deviation across replicated features
plotCtVariation(qPCRraw, variation="sd", log=TRUE)
# Summary of the first 40 genes, assuming there are 3 groups of samples
plotCtVariation(qPCRraw[1:40,], sample.reps=rep(1:2,3))

# Detailed summary of variation versus mean Ct value for replicated features within each sample
plotCtVariation(qPCRraw, type="detail", log=TRUE)
plotCtVariation(qPCRraw, type="detail")
# Add feature names to see which the highly varying replicates are.
plotCtVariation(qPCRraw, type="detail", add.featurenames=TRUE, pch=" ", cex=0.8)
# Use different information to indicate which features are replicates
plotCtVariation(qPCRraw, type="detail", feature.reps=paste("test", rep(1:96, each=4)))
# Examine variation across samples for the first 9 features
plotCtVariation(qPCRraw[1:9,], type="detail", sample.reps=paste("mutant", rep(1:3,2)), add.featurenames=TRUE)

# Examine the output
test <- plotCtVariation(qPCRraw, variation="sd")
names(test)
head(test[["Var"]])

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
