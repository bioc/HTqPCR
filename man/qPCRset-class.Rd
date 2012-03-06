\name{qPCRset-class}
\Rdversion{1.1}
\docType{class}
\alias{qPCRset-class}
\alias{[,qPCRset-method}
\alias{exprs,qPCRset-method}
\alias{exprs<-,qPCRset,ANY-method}
\alias{featureCategory}
\alias{featureCategory,qPCRset-method}
\alias{featureCategory<-}
\alias{featureCategory<-,qPCRset-method}
\alias{featureClass}
\alias{featureClass,qPCRset-method}
\alias{featureClass<-}
\alias{featureClass<-,qPCRset-method}
\alias{featureNames,qPCRset-method}
\alias{featureNames<-,qPCRset,character-method}
\alias{featurePos}
\alias{featurePos,qPCRset-method}
\alias{featurePos<-}
\alias{featurePos<-,qPCRset-method}
\alias{featureType}
\alias{featureType,qPCRset-method}
\alias{featureType<-}
\alias{featureType<-,qPCRset-method}
\alias{flag}
\alias{flag,qPCRset-method}
\alias{flag<-}
\alias{flag<-,qPCRset-method}
\alias{sampleNames,qPCRset-method}
\alias{sampleNames<-,qPCRset,character-method}
\alias{show,qPCRset-method}
\alias{summary,qPCRset-method}
\alias{getCt}
\alias{setCt<-}
\alias{n.wells}
\alias{n.samples}
\alias{getCtHistory}
\alias{setCtHistory<-}

\title{Class "qPCRset"}

\description{This is a class for containing the raw or normalized cycle threshold (Ct) values and some related quality information. It is suitable for TaqMan Low Density Arrays or any other type of (high-throughput) qPCR data, where gene expression is measured for any number of genes, across several samples/conditions. It inherits from \code{\link[Biobase:class.eSet]{eSet}} for microarray data.}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("qPCRset", assayData, phenoData, featureData, experimentData, annotation, protocolData, ...)} or using \code{readCtData}.
}

\section{Slots}{
  \describe{
    \item{\code{CtHistory}:}{Object of class \code{"data.frame"} indicating how the data has been read in, normalized, filtered etc. Gives the exact commands used during these operations.}
    \item{\code{assayData}:}{Object of class \code{"AssayData"}, containing the Ct values.}
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}, where information about samples can be added.}
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"}, where information about features can be added. If the object is from \code{readCtData}, the featureData will contain the columns 'featureName', 'featurePos' and 'featureType'. }
    \item{\code{experimentData}:}{Object of class \code{"MIAxE"}, where details about the experiment can be stored. }
    \item{\code{annotation}:}{Object of class \code{"character"}, where the name of the qPCR assay can be stored.}
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"}, where details of the protocol can be stored. }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}.}
  }
  Furthermore, the following information is contained within the object.
  \describe{
     \item{\code{flag}:}{Object of class \code{"data.frame"} containing the flag for each Ct value, as supplied by the input files.}
    \item{\code{featureCategory}:}{Object of class \code{"data.frame"} representing the quality of the measurement for each Ct value, such as "OK", "Undetermined" or "Unreliable" if the Ct value is considered too high.}
    }
}

\section{Extends}{
Class \code{"\linkS4class{eSet}"}, directly.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "eSet", distance 2.
Class \code{"\linkS4class{Versioned}"}, by class "eSet", distance 3.
}

\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "qPCRset")}: Subsets by genes or samples. }
    \item{exprs}{\code{signature(object = "qPCRset")}: Extracts the Ct matrix. Is identical to \code{getCt}}
    \item{exprs<-}{\code{signature(object = "qPCRset", value = "matrix")}: Replaces the Ct matrix. Is identical to \code{setCt<-}}
    \item{getCt}{\code{signature(object = "qPCRset")}: Extracts the Ct matrix. Is identical to \code{exprs}.}
    \item{setCt<-}{\code{signature(object = "qPCRset", value = "matrix")}: Replaces the Ct matrix.  Is identical to \code{exprs<-}.}
    \item{featureNames}{\code{signature(object = "qPCRset")}: Extracts the features (gene names) on the card. }
    \item{featureNames<-}{\code{signature(object = "qPCRset", value = "character")}: Replaces the features (gene names) on the card. }
    \item{sampleNames}{\code{signature(object = "qPCRset")}: Extracts the sample names.}
    \item{sampleNames<-}{\code{signature(object = "qPCRset", value = "character")}: Replaces the sample names. }
    \item{featureType}{\code{signature(object = "qPCRset")}: Extracts the different types of features on the card, such as controls and target genes.}
    \item{featureType<-}{\code{signature(object = "qPCRset", value = "factor")}: Replaces the feature type for each gene. }
    \item{featurePos}{\code{signature(object = "qPCRset")}: Extracts the position of each feature (gene) on the assay, representing the location "well" (such as well A1, A2, ...). If data does not come from a card format, the positions will be given consecutive names. }
    \item{featurePos<-}{\code{signature(object = "qPCRset", value = "character")}: Replaces the position of each feature (gene) on the card. }
    \item{featureClass}{\code{signature(object = "qPCRset")}: Extracts the feature class for each gene. }
    \item{featureClass<-}{\code{signature(object = "qPCRset", value = "factor")}: Replaces the feature class for each gene, for example if it is a marker, transcription factor or similar. }
    \item{featureCategory}{\code{signature(object = "qPCRset")}: Extracts the category of each Ct value.}
    \item{featureCategory<-}{\code{signature(object = "qPCRset", value = "data.frame")}: Replaces the category of each Ct value. }
    \item{flag}{\code{signature(object = "qPCRset")}: Extracts the flag of each Ct value. }
    \item{flag<-}{\code{signature(object = "qPCRset")}: Replaces the flag of each Ct value. }
     \item{n.wells}{\code{signature(object = "qPCRset")}: Extracts information about the number of wells on the card. }
     \item{n.samples}{\code{signature(object = "qPCRset")}: Extracts information about the number of samples in the set. }
     \item{getCtHistory}{\code{signature(object = "qPCRset")}: Extracts data frame containing information about the history of the object (which operations have been performed on it). }
     \item{setCtHistory<-}{\code{signature(object = "qPCRset")}: Add information about the history of the object. }
    \item{show}{\code{signature(object = "qPCRset")}: Displays some abbreviated information about the data object.}
    \item{summary}{\code{signature(object = "qPCRset")}: Displays a summary of the Ct values from each sample. }
	 }
}

\author{Heidi Dvinge}

\examples{
# The data format
data(qPCRraw)
show(qPCRraw)
getCtHistory(qPCRraw)
showClass("qPCRset")
str(qPCRraw)
# Information about samples
phenoData(qPCRraw)
pData(qPCRraw)
pData(qPCRraw)[,"Rep"] <- c(1,1,2,2,3,3)
# Information about features
featureData(qPCRraw)
head(fData(qPCRraw))
}

\keyword{classes}
