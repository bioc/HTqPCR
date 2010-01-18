#------------------------------------------------------------------
# Define the class
#------------------------------------------------------------------

setClass("qPCRset",
	representation(featureNames	= "character",
					sampleNames	= "character",
					exprs	= "matrix",
					flag	= "data.frame",
					featureType	= "factor",
					featurePos	= "character",
					featureClass	= "factor",
					featureCategory	= "data.frame",
					history	= "data.frame"),
	prototype = list(featureNames=character(), sampleNames=character(), exprs=matrix(), flag=data.frame(), featureType=factor(), featurePos=character(), featureClass=factor(), featureCategory=data.frame(), history=data.frame(stringsAsFactors=FALSE))
)


#------------------------------------------------------------------
# Functions also in Biobase
# Accessors 
#------------------------------------------------------------------

setMethod("exprs",
          signature(object="qPCRset"),
          function(object) {x<-object@exprs; 
          					rownames(x)<-featureNames(object); 
          					colnames(x)<-sampleNames(object); 
          					x})

setMethod("featureNames",
          signature(object="qPCRset"),
          function(object) object@featureNames)

setMethod("sampleNames",
          signature(object="qPCRset"),
          function(object) object@sampleNames)


#------------------------------------------------------------------
# Functions also in Biobase
# Replacement methods
#------------------------------------------------------------------

setReplaceMethod("exprs", 
			signature(object="qPCRset", value="matrix"),
			function(object, value) {object@exprs <- value; object})

setReplaceMethod("featureNames",
			signature(object="qPCRset", value="character"),
			function(object, value) {object@featureNames <- value; object})

setReplaceMethod("sampleNames",
			signature(object="qPCRset", value="character"),
			function(object, value) {object@sampleNames <- value; object})


#------------------------------------------------------------------
# Subsetting
#------------------------------------------------------------------

setMethod("[", "qPCRset", 
function(x, i, j, drop=FALSE) {
	if(!missing(i)) {
		exprs(x) <- exprs(x)[i, , drop=FALSE]
		featureNames(x) <- featureNames(x)[i]
		featureType(x) <- featureType(x)[i]
		featurePos(x) <- featurePos(x)[i]
		featureClass(x) <- featureClass(x)[i]
		featureCategory(x) <- featureCategory(x)[i, , drop=FALSE]
		flag(x)	<- flag(x)[i, , drop=FALSE]
	}
	if(!missing(j)) {
		exprs(x) <- exprs(x)[, j, drop=FALSE]
		sampleNames(x)	<- sampleNames(x)[j]
		featureCategory(x)	<- featureCategory(x)[, j, drop=FALSE]
		flag(x)	<- flag(x)[, j, drop=FALSE]
	}
	return(x)
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
	cat("Feature classes:\t", paste(levels(featureClass(object)), collapse=", "), "\n")
	cat("Feature categories:\t", paste(unique(unlist(featureCategory(object))), collapse=", "), "\n")
	cat("Sample names:\t\t", sampleNames(object)[1:3], "...\n")
})

setMethod("summary","qPCRset",
function(object) {
	s	<- summary(exprs(object))
	if (ncol(s) > 1) {
		rows	<- gsub("(.+): *[0-9\\.]+.+", "\\1", s[,1])
		s	<- apply(s[, 1:ncol(s), drop=FALSE], 2, function(x) gsub(".+:( *[0-9\\.]+).+", "\\1", x))
		rownames(s)	<- rows
	}
	colnames(s)	<- sampleNames(object)
	s	
})


#------------------------------------------------------------------
# Ct values (NB: Same as exprs() functions)
#------------------------------------------------------------------

getCt <- 
function(object) {x<-object@exprs; 
          	rownames(x)<-featureNames(object); 
          	colnames(x)<-sampleNames(object); 
          	x}

`setCt<-` <- 
function(object, value) {
			object@exprs <- value; object}


#------------------------------------------------------------------
# Feature position
#------------------------------------------------------------------

featurePos <-
function(object) object@featurePos

`featurePos<-` <-
function(object, value) {
			object@featurePos <- value; object}


#------------------------------------------------------------------
# Feature types
#------------------------------------------------------------------

featureType <-
function(object) object@featureType

`featureType<-` <-
function(object, value) {
			object@featureType <- value; object}


#------------------------------------------------------------------
# Feature categories
#------------------------------------------------------------------

featureCategory <-
function(object) object@featureCategory

`featureCategory<-` <-
function(object, value) {
			object@featureCategory <- value; object}


#------------------------------------------------------------------
# Feature classes
#------------------------------------------------------------------

featureClass <-
function(object) object@featureClass

`featureClass<-` <-
function(object, value) {
			object@featureClass <- value; object}


#------------------------------------------------------------------
# Flags
#------------------------------------------------------------------

flag <-
function(object) object@flag

`flag<-` <-
function(object, value) {
			object@flag <- value; object}


#------------------------------------------------------------------
# Size of qPCRset
#------------------------------------------------------------------

n.samples <-
function(object) ncol(object@exprs)

n.wells <-
function(object) nrow(object@exprs)


#------------------------------------------------------------------
# History
#------------------------------------------------------------------

history <-
function(object) object@history

