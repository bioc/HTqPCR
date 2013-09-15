\name{heatmapSig}
\Rdversion{1.1}
\alias{heatmapSig}
\title{Heatmap of deltadeltaCt values from qPCR data.}
\description{Heatmap and clustering of deltadeltaCt values from different sample comparisons using qPCR data.}

\usage{
heatmapSig(qDE, comparison = "all", col, zero.center = TRUE, mar, dist = "pearson", ...)
}

\arguments{
  \item{qDE}{data.frame or list, as created by \code{ttestCtData} or \code{limmaCtData}.}
  \item{comparison}{integers or the names of the comparisons to include in the plot. Defaults to all results in the \code{qDE} data, but a minimum of two is required.}
  \item{col}{colour scheme to use for the plot.}
  \item{zero.center}{logical, should the colour scale be centered around 0.}
  \item{mar}{vector of length two, the bottom and right side margins mof the heatmap.}
  \item{dist}{character string, either "pearson" (default) or "euclidean" indicating what type of distance is used for the clustering.}
  \item{\dots}{further arguments passed to \code{heatmap.2}.}
}

\details{
This function can be useful if multiple conditions are compared, for detecting features with similar behaviour in comparisons, and look at the general level of up and down regulation. }

\value{A plot if produced in the current graphics device.}

\author{Heidi Dvinge}

\seealso{\code{\link[gplots]{heatmap.2}} for modifying the plot, and \code{\link{ttestCtData}} or \code{\link{limmaCtData}} for generating the data used for the plotting.
}

\examples{
# Load example preprocessed data 
data(qPCRpros)
samples <- read.delim(file.path(system.file("exData", package="HTqPCR"), "files.txt"))
# Define design and contrasts
design <- model.matrix(~0+samples$Treatment)
colnames(design) <- c("Control", "LongStarve", "Starve")
contrasts	<- makeContrasts(LongStarve-Control, LongStarve-Starve, Starve-Control, levels=design)
# Reorder data to get the genes in consecutive rows
temp	<- qPCRpros[order(featureNames(qPCRpros)),]
# The actual test 
qDE	<- limmaCtData(temp, design=design, contrasts=contrasts, ndups=2, spacing=1)
# Plotting the heatmap
heatmapSig(qDE)
heatmapSig(qDE, dist="euclidean")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot }
