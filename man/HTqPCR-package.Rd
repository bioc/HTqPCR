\name{HTqPCR-package}
\Rdversion{1.1}
\alias{HTqPCR-package}
\alias{HTqPCR}
\docType{package}
\title{Analysis of High-Throughput qPCR data (HTqPCR)}

\description{This package is for analysing high-throughput qPCR data. Focus is on data from Taqman Low Density Arrays, but any kind of qPCR performed across several samples is applicable. Cycle threshold (Ct) data from different cards (samples) is read in, normalised, processed and the genes are tested for differential expression across different samples. Results are visualised in various ways.
}

\details{
\tabular{ll}{
Package: \tab HTqPCR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2009-07-03\cr
License: \tab Artistic\cr
LazyLoad: \tab yes\cr
Depends: \tab methods\cr
}
}

\author{
Maintainer: Heidi Dvinge <heidi@ebi.ac.uk> 
Maintainer: Paul Bertone <bertone@ebi.ac.uk>
}

%\references{
%~~ Literature or other references for background information ~~
%}

\keyword{ package }

%\seealso{
%\code{\link[limma:limma-package]{limma}}
%}

\examples{
# Locate example data and create qPCRset object
exPath <- system.file("exData", package="HTqPCR")
exFiles <- read.delim(file.path(exPath, "files.txt"))
raw <- readCtData(files=exFiles$File, path=exPath)
# Preprocess
raw.cats	<- setCategory(raw, groups=exFiles$Treatment, plot=FALSE)
norm	<- normalizeCtData(raw.cats, norm="scale.rank")
# Various plots
plotCtDensity(norm)
plotCtBoxes(norm)
plotCtOverview(norm, groups=exFiles$Treatment, genes=featureNames(raw)[1:10], calibrator="Control")
plotCtCor(norm)
plotCtScatter(norm, cards=c(1,4), col="type")
# Define design and contrasts for testing differential expression
design <- model.matrix(~0+exFiles$Treatment)
colnames(design) <- c("Control", "LongStarve", "Starve")
contrasts	<- makeContrasts(LongStarve-Control, LongStarve-Starve, Starve-Control, levels=design)
# Reorder by featureNames (2 replicates of each feature) and the actual test
norm2	<- norm[order(featureNames(norm)),]
diff.exp <- limmaCtData(norm2, design=design, contrasts=contrasts, ndups=2, spacing=1)
# Some of the results
names(diff.exp)
diff.exp[["LongStarve - Control"]][1:10,]
diff.exp[["Summary"]][1:10,]
# Some plots of results
plotCtRQ(diff.exp, genes=1:10)
plotCtSignificance(qDE=diff.exp, q=norm2, groups=exFiles$Treatment, calibrator="Control", target="LongStarve", genes=1:10, jitter=0.2)
plotCtSignificance(qDE=diff.exp, q=norm2, comparison="LongStarve - Starve", groups=exFiles$Treatment, calibrator="Starve", target="LongStarve", genes=1:10, jitter=0.2)
}
