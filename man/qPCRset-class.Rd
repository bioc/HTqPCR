\name{qPCRset-class}
\Rdversion{1.1}
\docType{class}
\alias{qPCRset-class}
\alias{[,qPCRset-method}
\alias{exprs,qPCRset-method}
\alias{exprs<-,qPCRset,matrix-method}
\alias{featureNames,qPCRset-method}
\alias{featureNames<-,qPCRset,character-method}
\alias{sampleNames,qPCRset-method}
\alias{sampleNames<-,qPCRset,character-method}
\alias{featureClass}
\alias{featureClass<-}
\alias{featureCategory}
\alias{featureCategory<-}
\alias{featurePos}
\alias{featurePos<-}
\alias{featureType}
\alias{featureType<-}
\alias{flag}
\alias{flag<-}
\alias{n.wells}
\alias{n.samples}
\alias{normalized}
\alias{show,qPCRset-method}
\alias{summary,qPCRset-method}

\title{Class "qPCRset"}

\description{This is a class for containing the raw or normalized data, such as Ct values and some related quality information, from TaqMan Low Density Arrays or any other type of (high-throughput) qPCR data. It is similar to \code{\link[Biobase:class.eSet]{eSet}} for microarray data.}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("qPCRset", ...)} or using \code{readCtData}.
}

\section{Slots}{
	 \describe{
    \item{\code{featureNames}:}{Object of class \code{"character"} giving the names of the features, such as genes or miRNAs, in the samples.}
    \item{\code{sampleNames}:}{Object of class \code{"character"} containing the sample names.}
    \item{\code{exprs}:}{Object of class \code{"matrix"} containing the Ct values.}
    \item{\code{flag}:}{Object of class \code{"data.frame"} containing the flag for each Ct value, as supplied by the input files.}
    \item{\code{featureType}:}{Object of class \code{"factor"} representing the different types of features on the card, such as controls and target genes.}
    \item{\code{featurePos}:}{Object of class \code{"character"} representing the location "well" of a gene on the card. If data does not come from a card format, the positions will be given consecutive names.}
    \item{\code{featureClass}:}{Object of class \code{"factor"} with some meta-data about the genes, for example if it is a marker, transcription factor or similar.}
    \item{\code{featureCategory}:}{Object of class \code{"data.frame"} representing the quality of the measurement for each Ct value, such as "OK", or "Unreliable" if the Ct value is considered too migh.}
    \item{\code{normalized}:}{Object of class \code{"character"} indicating if the data has been normalized, and if so then what method was used.}
  }
}

\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "qPCRset")}: Subsets by genes or samples. }
    \item{exprs}{\code{signature(object = "qPCRset")}: Extracts the Ct matrix. }
    \item{exprs<-}{\code{signature(object = "qPCRset", value = "matrix")}: Replaces the Ct matrix. }
    \item{featureNames}{\code{signature(object = "qPCRset")}: Extracts the features (gene names) on the card. }
    \item{featureNames<-}{\code{signature(object = "qPCRset", value = "character")}: Replaces the features (gene names) on the card. }
    \item{sampleNames}{\code{signature(object = "qPCRset")}: Extracts the sample names.}
    \item{sampleNames<-}{\code{signature(object = "qPCRset", value = "character")}: Replaces the sample names. }
    \item{featureType}{\code{signature(object = "qPCRset")}: Extracts the feature type for each gene. }
    \item{featureType<-}{\code{signature(object = "qPCRset", value = "factor")}: Replaces the feature type for each gene. }
    \item{featurePos}{\code{signature(object = "qPCRset")}: Extracts the position of each feature (gene) on the card. }
    \item{featurePos<-}{\code{signature(object = "qPCRset", value = "character")}: Replaces the position of each feature (gene) on the card. }
    \item{featureClass}{\code{signature(object = "qPCRset")}: Extracts the feature class for each gene. }
    \item{featureClass<-}{\code{signature(object = "qPCRset", value = "factor")}: Replaces the feature class for each gene. }
    \item{featureCategory}{\code{signature(object = "qPCRset")}: Extracts the category of each Ct value. }
    \item{featureCategory<-}{\code{signature(object = "qPCRset", value = "data.frame")}: Replaces the category of each Ct value. }
     \item{n.wells}{\code{signature(object = "qPCRset")}: Extracts information about the number of wells on the card. }
     \item{n.samples}{\code{signature(object = "qPCRset")}: Extracts information about the number of samples in the set. }
     \item{normalized}{\code{signature(object = "qPCRset")}: Extracts information about normalization. }
     \item{normalized}{\code{signature(object = "qPCRset")}: Extracts information about normalization. }
    \item{show}{\code{signature(object = "qPCRset")}: Displays som abbreviated information about the data object.}
    \item{summary}{\code{signature(object = "qPCRset")}: Displays a summary of the Ct values from each sample. }
	 }
}

\author{Heidi Dvinge}

\examples{
data(qPCRraw)
show(qPCRraw)
showClass("qPCRset")
str(qPCRraw)
}

\keyword{classes}
