\name{readCtData}
\Rdversion{1.1}
\alias{readCtData}
\alias{.readCtPlain}
\alias{.readCtSDS}
\alias{.readCtBioMark}
\alias{.readCtLightCycler}
\alias{.readCtCFX}
\alias{.readCtOpenArray}
\title{Reading Ct values from qPCR experiments data into a qPCRset}
\description{
This function will read tab separated text files with Ct values and feature meta-data from high-throughput qPCR experiments into a qPCRset containing all the relevant information.
}

\usage{
readCtData(files, path = NULL, n.features = 384, format="plain", column.info, flag, feature, type, position, Ct, header = FALSE, SDS = FALSE, n.data = 1, samples, na.value = 40, sample.info, ...)
}

\arguments{
  \item{files}{ character vector with the names of the files to be read.}
  \item{path}{ character string with the path to the folder containing the data files.}
  \item{n.features}{integer, number of features present on each array (e.g. 384). See details.}
  \item{format}{character, the format of the input file. Options are "plain", "SDS", "LightCycler", "CFX", "OpenArray" and "BioMark". See Details. }
  \item{column.info}{list, indicating which column number or name the information of interest is in. It is set automatically by \code{format}, but this can be overridden manually. The names list slots can be 'flag', 'feature', 'position', 'type and 'Ct'. See Details. Note than when indicating column names, these are sometimes changed by R to be syntactically valid, so replacing e.g. brackets by dots. }
  \item{flag,feature,
  Ct,type,
  position}{deprecated, use \code{column.info} instead.}
  \item{header}{logical, does the file contain a header row or not. Only used for \code{format="plain"}.}
  \item{SDS}{deprecated, use \code{format="SDS"} instead.}
  \item{n.data}{integer vector, same length as \code{files}. Indicates the number of samples that are present in each file. For each file in \code{files}, n.data*n.features lines will be read.}
 \item{samples}{ character vector with names for each sample. Per default the file names are used.}
  \item{na.value}{integer, a Ct value that will be assigned to all undetermined/NA wells.}
  \item{sample.info}{object of class \code{AnnotatedDataFrame}, given the phenoData of the object. Can be added later.}
  \item{\dots}{ any other arguments are passed to \code{\link{read.table}} or \code{\link{read.csv}}.}
}

\details{
This is the main data input function for the HTqPCR package for analysing qPCR data. It extracts the threshold cycle, Ct value, of each well on the card, as well as information about the quality (e.g.~passed/failed) of the wells. The function is tuned for data from TaqMan Low Density Array cards, but can be used for any kind of qPCR data.

The information to be extracted is:
\itemize{
	\item{\code{flag}}{ integer indicating the number of column containing information about the flags.}
	\item{\code{feature}}{ integer indicating the number of column containing information about the individual features (typically gene names). }
	\item{\code{type}}{ integer indicating the number of column containing information about the type of each feature.}
	\item{\code{position}}{ integer indicating the number of column containing information about the position of features on the card.}
	\item{\code{Ct}}{ integer indicating the number of column containing information about the Ct values.}
Per default, this information is assumed to be in certain columns depending on the input format.  
}

\code{featureNames}, \code{featureType} and \code{featurePos} will be extracted from the first file. If \code{flag}, \code{type} or \code{position} are not included into \code{column.info}, this means that this information is not available in the file. \code{flag} will then be set to "Passed", \code{type} to "Target" and \code{position} to "feature1", "feature2", ... etc until the end of the file. Especially \code{position} might not be available in case the data does not come from a card format, but it is required in subsequent functions in order to disambiguate in case some features are present multiple times.

\code{format} indicates the format of the input file. The options currently implemented are:
\itemize{
	\item{\code{plain}}{ A tab-separated text file, containing no header unless \code{header=TRUE}. The information extracted defaults to \code{column.info=list(flag=4, feature=6, type=7, position=3, Ct=8)}.}
	\item{\code{SDS}}{ An output file from the SDS Software program. This is often used for the TaqMan Low Density Arrays from Applied Biosystems, but can also be used for assays from other vendors, such as Exiqon. \code{column.info} is the same as for "plain".}
	\item{\code{OpenArray}}{ The TaqMan OpenArray Real-Time PCR Plates. The information extracted defaults to \code{column.info=list(flag="ThroughHole.Outlier", feature="Assay.Assay.ID", type="Assay.Assay.Type", position="ThroughHole.Address", Ct="ThroughHole.Ct")}.}
	\item{\code{BioMark}}{ The BioMark HD System from Fluidigm, currently including the 48.48 and 96.96 assays. The information extracted defaults to \code{column.info=list(flag="Call", feature="Name.1", position="ID", Ct="Value")}.}
	\item{\code{CFX}}{ The CFX Automation System from Bio-Rad. The information extracted defaults to \code{column.info=list(feature="Content", position="Well", Ct="Cq.Mean")}.}
	\item{\code{LightCycler}}{ The LightCycler System from Roche Applied Science. The information extracted defaults to \code{column.info=list(feature="Name", position="Pos", Ct="Cp")}.}
}
The BioMark and OpenArray assays always contain information multiple samples on each assay, such as 48 features for 48 samples for the BioMark 48.48. The results across these samples are always present in a single file, e.g. with 48x48=2304 rows. Setting n.features=2304 will read in all the information and create a \code{qPCRset} object with dimensions 2304x1. Setting n.data=48 and n.features=48 will however automatically convert this into a 48x48 \code{qPCRset}. See \code{openVignette(package="HTqPCR")} for examples. The samples are being read in the order in which they're present in the file, i.e. from row 1 onwards, regardless of how they're loaded onto the particular platform.

If the data was analysed using for example SDS Software it may contain a variable length header specifying parameters for files that were analysed at the same time. If \code{format="SDS"} then \code{readCtData} will scan through the first 100 lines of each file, and skip all lines until (and including) the line beginning with "#", which is the header. The end of the file might also contain some plate ID information, but only the number of lines specified in \code{n.features} will be read.

\code{n.features} indicates the number of features present on each array. For example, for a 384 well plate with just 1 sample, the number would be 382. For a plate with 2 individual samples loaded onto it, \code{n.features=196} and \code{n.data=2}. For 1 file with 5 plates and 2 samples per plate, the numbers are \code{n.features=196} and \code{n.data=10}. n.features*n.data must correspond to the total number of lines to be read from each file.
}                               

\section{Warnings}{
The files are all assumed to belong to the same design, i.e.~have the same features (genes) in them and in identical order.}

\value{
A \code{"\link[=qPCRset-class]{qPCRset}"} object. }

\author{ Heidi Dvinge}

\seealso{\code{\link{read.delim}} for further information about reading in data, and \code{"\link[=qPCRset-class]{qPCRset}"} for a definition of the resulting object.}

\examples{
# Locate example data and create qPCRset object
exPath <- system.file("exData", package="HTqPCR")
exFiles <- read.delim(file.path(exPath, "files.txt"))
raw <- readCtData(files=exFiles$File, path=exPath)
# Example of adding missing information (random data in this case)
featureClass(raw)	<- factor(rep(c("A", "B", "C"), each=384/3))
pData(raw)[,"rep"] 	<- c(1,1,2,2,3,3)

## See the package vignette for more examples, including different input formats.

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}
