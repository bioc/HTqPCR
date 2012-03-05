\name{qPCRraw}
\Rdversion{1.1}
\docType{data}
\alias{qPCRraw}
\title{Example raw qPCR data.}
\description{Six qPCR samples, performed on the TaqMan Low Density Arrays from Applied Biosystem. Each sample contains 384 PCR reactions, and there are 3 different samples with 2 replicates each. To be used as example data in the HTqPCR package.}

\usage{data(qPCRraw)}

\format{
  The format is:
Formal class 'qPCRset' [package ".GlobalEnv"] with 8 slots
  ..@ CtHistory        :'data.frame':	1 obs. of  1 variable:
  .. ..$ history: chr "readCtData(files = exFiles$File, path = exPath)"
  ..@ assayData        :<environment: 0x118094e30> 
  ..@ phenoData        :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	1 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr "Sample numbering"
  .. .. ..@ data             :'data.frame':	6 obs. of  1 variable:
  .. .. .. ..$ sample: int [1:6] 1 2 3 4 5 6
  .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ featureData      :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	4 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr [1:4] NA NA NA NA
  .. .. ..@ data             :'data.frame':	384 obs. of  4 variables:
  .. .. .. ..$ featureNames: Factor w/ 191 levels "Gene1","Gene10",..: 1 104 115 126 137 148 159 170 181 2 ...
  .. .. .. ..$ featureType : Factor w/ 2 levels "Endogenous Control",..: 1 2 2 2 2 2 2 2 2 2 ...
  .. .. .. ..$ featurePos  : Factor w/ 384 levels "A1","A10","A11",..: 1 12 18 19 20 21 22 23 24 2 ...
  .. .. .. ..$ featureClass: Factor w/ 3 levels "Kinase","Marker",..: 1 2 1 3 2 2 2 3 1 2 ...
  .. .. ..@ dimLabels        : chr [1:2] "featureNames" "featureColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ experimentData   :Formal class 'MIAME' [package "Biobase"] with 13 slots
  .. .. ..@ name             : chr ""
  .. .. ..@ lab              : chr ""
  .. .. ..@ contact          : chr ""
  .. .. ..@ title            : chr ""
  .. .. ..@ abstract         : chr ""
  .. .. ..@ url              : chr ""
  .. .. ..@ pubMedIds        : chr ""
  .. .. ..@ samples          : list()
  .. .. ..@ hybridizations   : list()
  .. .. ..@ normControls     : list()
  .. .. ..@ preprocessing    : list()
  .. .. ..@ other            : list()
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 2
  .. .. .. .. .. ..$ : int [1:3] 1 0 0
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ annotation       : chr(0) 
  ..@ protocolData     :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	0 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr(0) 
  .. .. ..@ data             :'data.frame':	6 obs. of  0 variables
  .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. ..@ .Data:List of 3
  .. .. .. ..$ : int [1:3] 2 14 0
  .. .. .. ..$ : int [1:3] 2 14 0
  .. .. .. ..$ : int [1:3] 1 3 0
}

\examples{
data(qPCRraw)
}

\keyword{datasets}
