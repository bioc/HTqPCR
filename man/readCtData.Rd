\name{readCtData}
\Rdversion{1.1}
\alias{readCtData}
\title{Reading Ct values from qPCR experiments data into a qPCRset}
\description{
This function will read tab separated text files with Ct values and feature meta-data from high-throughput qPCR experiments into a qPCRset containing all the relevant information.
}

\usage{
readCtData(files, path = NULL, n.features = 384, flag = 4, feature = 6, type = 7, position = 3, Ct = 8, header = FALSE, SDS = FALSE, n.data = 1, samples, na.value = 40, ...)
}

\arguments{
  \item{files}{ character vector witht he names of the files to be read.}
  \item{path}{ character string with the path to the folder containing the data files.}
  \item{n.features}{integer, number of features present in each file.}
  \item{flag}{ integer indicating the number of column containing information about the flags. See Details.}
  \item{feature}{ integer indicating the number of column containing information about the individual features (typically gene names). }
  \item{type}{ integer indicating the number of column containing information about the type of each feature. See Details.}
  \item{position}{integer indicating the number of column containing information about the position of features on the card. See Details.}
  \item{Ct}{integer indicating the number of column containing information about the Ct values.}
  \item{header}{logical, does the file contain a header row or not.}
  \item{SDS}{logical, is the data in the output format from Sequence Detection Systems (SDS) Software. See Details.}
  \item{n.data}{integer vector, same length as \code{files}. Indicates the number of samples that are present in each file. For each file n.data*n.features lines will be read.}
 \item{samples}{ character vector with names for each sample. Per default the file names are used.}
  \item{na.value}{integer, a Ct value that will be assigned to all undetermined/NA wells.}
  \item{\dots}{ any other arguments are passed to \code{\link{read.table}}.}
}

\details{
This is the main data input function for the HTqPCR package for analysing qPCR data. It extracts the threshold cycle, Ct value, of each well on the card, as well as information about the quality (e.g.~passed/failed) of the wells. The function is tuned for data from TaqMan Low Density Array cards, but can be used for any kind of qPCR data.

\code{featureNames}, \code{featureType} and \code{featurePos} will be extracted from the first file. If \code{flag}, \code{type} or \code{position} er set to NULL, this means that this information is not available in the file. \code{flag} will then be set to "Passed", \code{type} to "Target" and \code{position} to "feature1", "feature2", ... etc until the end of the file. Especially \code{position} might not be available in case the data does not come from a card format, but it is required in subsequent functions in order to disambiguate between features in case some features are present multiple times.

If the data was analysed using SDS Software it may contain a variable length header specifying parameters for files that were analysed at the same time. If \code{SDS=TRUE} then \code{readCtData} will scan through the first 100 lines of each file, and skip all lines until (and including) the line beginning with "#", which is the header. The end of the file might also contain some plate ID information, but only the number of lines specified in \code{n.features} will be read.
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
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}
