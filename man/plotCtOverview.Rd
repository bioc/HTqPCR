\name{plotCtOverview}
\Rdversion{1.1}
\alias{plotCtOverview}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Overview plot of qPCR Ct values across multiple conditions.}

\description{Function for high-throughput qPCR data, for showing the average Ct values for features in a barplot, either for individual samples or averaged across biological or technical groups. If Ct values are shown, error bars can be included, or the Ct values can be displayed relative to a calibrator sample.}

\usage{
plotCtOverview(q, cards = TRUE, genes, groups, calibrator, replicates = TRUE, col, conf.int = FALSE, legend = TRUE, ...)
}

\arguments{
  \item{q}{object of class qPCRset.}
  \item{cards}{integer, the cards (samples) to use. Defaults to all.}
  \item{genes}{vector selecting the features to show. See Details.}
  \item{groups}{vector with groups to average the samples across. If missing all the samples are displayed individually.}
  \item{calibrator}{the value in \code{groups} to use as calibrator sample. See Details.}
  \item{replicates}{logical, if should values from replicated features in each sample be collapsed or kept separate.}
  \item{col}{colours to use for each sample or group.}
  \item{conf.int}{logical, should the 95 percent confidence interval be shown. See Details.}
  \item{legend}{logical, should a legend be included in the plot.}
  \item{\dots}{further arguments passed to \code{barplot}.}
}

\details{
If a calibrator is chosen all values will be displayed relative to this, i.e. as Ct(sample)-Ct(calibrator). If there is no calibrator, the full Ct values are shown, including 95\% confidence interval if selected. For confidence intervals when there is a calibrator, it's the variation across Ct(sample)-average(Ct(calibrator)) that is shown.	

When setting \code{replicates=TRUE} it is often better to specify \code{genes} by name rather than selecting for example the first 10 features using 1:10. This literally only takes the first 10 rows of the data, although some of these features might be replicated elsewhere in the data.
}

\value{A figure is produced in the current graphics device.}

\author{Heidi Dvinge}

\examples{
# Load example data
data(qPCRraw)
exPath <- system.file("exData", package="HTqPCR")
samples	<- read.delim(file.path(exPath, "files.txt"))
# Show all samples for the first 10 genes
g	<- featureNames(qPCRraw)[1:10]
plotCtOverview(qPCRraw, genes=g, xlim=c(0,90))
plotCtOverview(qPCRraw, genes=g, xlim=c(0,50), groups=samples$Treatment)
plotCtOverview(qPCRraw, genes=g, xlim=c(0,60), groups=samples$Treatment, conf.int=TRUE, ylim=c(0,55))
# Relative to a calibrator sample
plotCtOverview(qPCRraw, genes=g, groups=samples$Treatment, calibrator="Control")
plotCtOverview(qPCRraw, genes=g, groups=samples$Treatment, calibrator="Control", conf.int=TRUE, ylim=c(-0.5,0.5))
plotCtOverview(qPCRraw, genes=g, groups=samples$Treatment, calibrator="LongStarve")
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot }
