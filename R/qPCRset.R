#------------------------------------------------------------------
# Define the class
#------------------------------------------------------------------

setClass("qPCRset", contains = "eSet", representation(CtHistory="data.frame"))


#------------------------------------------------------------------
# Functions also in Biobase
# Accessors
#------------------------------------------------------------------

setMethod("exprs", signature(object="qPCRset"), definition =
 function (object) {x <- assayDataElement(object, "exprs"); rownames(x) <- featureNames(object); x}
)

setMethod("featureNames",
  signature(object="qPCRset"),
  function(object) as.character(fData(object)$featureNames))

setMethod("sampleNames",
  signature(object="qPCRset"),
  function(object) colnames(exprs(object)))


#------------------------------------------------------------------
# Functions also in Biobase
# Replacement methods
#------------------------------------------------------------------

setReplaceMethod("exprs", signature(object="qPCRset"), definition =
 function (object, value) assayDataElementReplace(object, "exprs", value, validate=TRUE)
)

setReplaceMethod("featureNames",
			signature(object="qPCRset", value="character"),
			function(object, value) {fData(object)[,"featureNames"]<- value; object})

setReplaceMethod("sampleNames",
			signature(object="qPCRset", value="character"),
			function(object, value) {colnames(exprs(object))<-value; object})


#------------------------------------------------------------------
# Subsetting
#------------------------------------------------------------------

setMethod("[", "qPCRset",
function(x, i=TRUE, j=TRUE, drop=FALSE) {
	# if(!missing(i)) {
		# exprs(x) <- exprs(x)[i, , drop=FALSE]
		# if (nrow(featureCategory(x))>0)
			# featureCategory(x) <- featureCategory(x)[i, , drop=FALSE]
		# if (nrow(flag(x))>0)
			# flag(x) <- flag(x)[i, , drop=FALSE]
		# featureData(x) <- featureData(x)[i,,drop=FALSE]
	# }
	# if(!missing(j)) {
		# exprs(x) <- exprs(x)[, j, drop=FALSE]
		# if (ncol(featureCategory(x))>0)
			# featureCategory(x) <- featureCategory(x)[, j, drop=FALSE]
		# if (ncol(flag(x))>0)
			# flag(x) <- flag(x)[, j, drop=FALSE]
		# phenoData(x) <- phenoData(x)[j,,drop=FALSE]
	# }
	# return(x)
	out <- new("qPCRset", 
				exprs=exprs(x)[i, j, drop=FALSE], 
				featureCategory=featureCategory(x)[i, j, drop=FALSE], 
				flag=flag(x)[i, j, drop=FALSE])
	phenoData(out) <- phenoData(x)[j,,drop=FALSE]
	featureNames(out) <- featureNames(x)[i]
	out
})


#------------------------------------------------------------------
# Various displays and similar things
#------------------------------------------------------------------

setMethod("show","qPCRset",
function(object) {
	cat("An object of class \"", class(object), "\"\n",sep="")
	cat("Size: ", length(featureNames(object)), "features,", ncol(exprs(object)), "samples\n")
	cat("Feature types:\t\t", paste(levels(featureType(object)), collapse=", "), "\n")
	cat("Feature names:\t\t", featureNames(object)[1:3], "...\n")
	cat("Feature classes:\t\t", paste(levels(featureClass(object)), collapse=", "), "\n")
	cat("Feature categories:\t", paste(unique(unlist(featureCategory(object))), collapse=", "), "\n")
	n <- ifelse(length(sampleNames(object))>2, 3, length(sampleNames(object)))
	cat("Sample names:\t\t", sampleNames(object)[1:n], "...\n")
})

setMethod("summary","qPCRset",
function(object) {
	s <- summary(exprs(object))
	if (ncol(s) > 1) {
		rows <- gsub("(.+): *[0-9\\.]+.+", "\\1", s[,1])
		s <- apply(s[, 1:ncol(s), drop=FALSE], 2, function(x) gsub(".+:( *[0-9\\.]+).+", "\\1", x))
		rownames(s) <- rows
	}
	colnames(s) <- sampleNames(object)
	s	
})


#------------------------------------------------------------------
# Ct values (NB: Same as exprs() functions)
#------------------------------------------------------------------

getCt <-
function(object) {exprs(object)}

`setCt<-` <-
function (object, value) {exprs(object)<-value; object}


#------------------------------------------------------------------
# Feature position
#------------------------------------------------------------------

setMethod("featurePos", signature = "qPCRset", definition =
 function (object) as.character(fData(object)$featurePos)
)

setReplaceMethod("featurePos", signature = "qPCRset", definition =
 function (object, value) {fData(object)[,"featurePos"] <- value; object}
)


#------------------------------------------------------------------
# Feature types
#------------------------------------------------------------------

setMethod("featureType", signature = "qPCRset", definition =
 function (object) as.character(fData(object)$featureType)
)

setReplaceMethod("featureType", signature = "qPCRset", definition =
 function (object, value) {fData(object)[,"featureType"] <- value; object}
)


#------------------------------------------------------------------
# Feature categories
#------------------------------------------------------------------

setMethod("featureCategory", signature = "qPCRset", definition =
 function (object) {x <- assayDataElement(object, "featureCategory"); rownames(x) <- make.unique(featureNames(object)); x}
# function (object) {x <- assayDataElement(object, "featureCategory"); x}
)

setReplaceMethod("featureCategory", signature = "qPCRset", definition =
 function (object, value) assayDataElementReplace(object, "featureCategory", value, validate=FALSE)
)


#------------------------------------------------------------------
# Feature classes
#------------------------------------------------------------------

setMethod("featureClass", signature = "qPCRset", definition =
 function (object) as.character(fData(object)$featureClass)
)

setReplaceMethod("featureClass", signature = "qPCRset", definition =
 function (object, value) {fData(object)[,"featureClass"] <- value; object}
)


#------------------------------------------------------------------
# Flags
#------------------------------------------------------------------

setMethod("flag", signature = "qPCRset", definition =
 function (object) {x <-assayDataElement(object, "flag"); rownames(x) <- make.unique(featureNames(object)); x}
# function (object) {x <-assayDataElement(object, "flag"); x}
)

setReplaceMethod("flag", signature = "qPCRset", definition =
 function (object, value) assayDataElementReplace(object, "flag", value, validate=FALSE)
)


#------------------------------------------------------------------
# Size of qPCRset
#------------------------------------------------------------------

n.samples <-
function(object) ncol(exprs(object))

n.wells <-
function(object) nrow(exprs(object))


#------------------------------------------------------------------
# History
#------------------------------------------------------------------

getCtHistory <-
function(object) object@CtHistory

`setCtHistory<-` <-
function(object, value) {object@CtHistory <- value; object}
