readCtData	<- 
function(files,
	path	= NULL,
	n.features	= 384,
	format	= "plain",
	column.info,
	flag,
	feature,
	type,
	position,
	Ct,
	header	= FALSE,
	SDS	= FALSE,
	n.data	= 1,
	samples,
	na.value	= 40,
	sample.info,
	...)
{
	# Initial checks
	if(missing(files))
		stop("No input files specified")
	# Various parameters
	if (length(files)!=length(n.data) & length(n.data)!=1)
		stop("n.data must either be a single integer, or same length as number of files")
	if (length(n.data)==1) 
		n.data	<- rep(n.data, length(files))
	nsamples	<- sum(n.data)
	ncum	<- cumsum(n.data)
	s.names	<- NULL
	nspots	<- n.features
	# Warn if some 'old style' parameters are set
	if (SDS) {
		warning("Please use format='SDS'. The SDS' parameter is retained for backward compatibility only.")
		format	<- "SDS"
	}
	if (!missing(flag) | !missing(feature) | !missing(type) | !missing(position) | !missing(Ct)) {
		warning("Please use 'column.info' for providing a list of column numbers containing particular information. The use of 'flag', 'feature', 'type', 'position' and 'Ct' is deprecated and will be removed in future versions.")
	}
	# Which values to set for information contained within samples
	if (missing(column.info)) {
		column.info	<- switch(format,
			plain	= list(flag=4, feature=6, type=7, position=3, Ct=8),
			SDS	= list(flag=4, feature=6, type=7, position=3, Ct=8),
			LightCycler	= list(feature="Name", position="Pos", Ct="Cp"),
			CFX	= list(feature="Content", position="Well", Ct="C.t..Mean"),
			OpenArray	= list(flag="ThroughHole.Outlier", feature="Assay.Assay.ID", type="Assay.Assay.Type", position="ThroughHole.Address", Ct="ThroughHole.Ct"),
			BioMark	= list(flag="Call", feature="Name.1", position="ID", Ct="Value"))
	}
	# Initializing data required for qPCRset object.
	X <- matrix(0,nspots,nsamples)
	X.flags	<- as.data.frame(X)
	X.cat <- data.frame(matrix("OK", ncol=nsamples, nrow=nspots), stringsAsFactors=FALSE)
	# Read through all files (assume features are sorted based on position)
	for (i in seq_along(files)) {
		# Which columns of the qPCRset object are to be filled out
		if (i==1) {
			cols	<- 1:ncum[i]
		} else {
			cols	<- (ncum[i-1]+1):ncum[i]
		}
		# The file to read
		readfile	<- ifelse(is.null(path), files[i], file.path(path, files[i]))
		# Format dependent reading of files
		sample	<- switch(format,
			plain	= .readCtPlain(readfile=readfile, header=header, n.features=n.features, n.data=n.data, i=i, ...),
			SDS	= .readCtSDS(readfile=readfile, n.data=n.data, i=i, nspots=nspots, ...),
			LightCycler	= .readCtLightCycler(readfile=readfile, n.data=n.data, i=i, nspots=nspots, ...),
			CFX	= .readCtCFX(readfile=readfile, n.data=n.data, i=i, nspots=nspots, ...),
			OpenArray = .readCtOpenArray(readfile=readfile, n.data=n.data, i=i, nspots=nspots, ...),
			BioMark	= .readCtBioMark(readfile=readfile, n.data=n.data, i=i, nspots=nspots, ...))
		# Assign Ct values to object
		data	<- matrix(sample[,column.info[["Ct"]]], ncol=n.data[i])
		undeter	<- apply(data, 2, function(x) x %in% c("Undetermined", "No Ct"))
  		X.cat[,cols][undeter] <- "Undetermined"
		if (is.null(na.value)) {
 			data[data=="Undetermined"]	<- NA
 		} else {
 			data[data %in% c("Undetermined", "No Ct", "999")] <- na.value
 		}
		X[,cols]	<- apply(data, 2, function(x) as.numeric(as.character(x)))
		if ("flag" %in% names(column.info)) {
			flags	<-  matrix(sample[,column.info[["flag"]]], ncol=n.data[i])
			flags[flags=="-"]	<- "Failed"
			flags[flags=="+"]	<- "Passed"
			X.flags[,cols]	<- flags
		} else {
			X.flags[,cols]	<- "Passed"
		}
		# Get sample names if present
		if (format=="OpenArray") {
			s.names	<- c(s.names, unique(sample$SampleInfo.SampleID))
		} else if (format %in% c("plain", "SDS")) {
			s.names	<- c(s.names, unique(sample[,2]))
		} else {
			s.names	<- s.names
		}
		# Assemble the feature data info to the qPCRset
		if (i==1) {
			# Add default values if any are missing from file
			featPos	<- paste("feature", 1:nspots, sep="")
			if ("position" %in% names(column.info))
				featPos	<- as.character(sample[1:nspots,column.info[["position"]]])
			featType	<- factor(rep("Target", nspots))
			if ("type" %in% names(column.info))		
				featType	<- sample[1:nspots,column.info[["type"]]]
			featName	<- paste("feature", 1:nspots, sep="")
			if ("feature" %in% names(column.info)) 
				featName	<- as.character(sample[1:nspots,column.info[["feature"]]])
			# Create feature data object
			df <- data.frame(featureNames=featName, featureType=as.factor(featType), featurePos=featPos)
			metaData <- data.frame(labelDescription=c(
               	"Name of the qPCR feature (gene)",
               	"Type pf feature",
               	"Position on assay"))
			featData	<- AnnotatedDataFrame(data=df, varMetadata=metaData)
		}
	}
	# Provide sample names
	if (!missing(samples)) {
		if (length(samples) < nsamples) {
			warning("Not enough sample names provided; using Sample1, Sample2, ... instead\n")
			samples	<- paste("Sample", 1:nsamples, sep="")
		} else if (length(samples) == nsamples) {
			samples	<- samples	
		}
	} else if (missing(samples)) {
		if (length(files)==nsamples) {
			samples	<- gsub("(.+)\\..+", "\\1", files)
		} else if (length(s.names)==nsamples) {
			samples	<- s.names
		} else {
			samples	<- paste("Sample", 1:nsamples, sep="")
		}
	}
	# In case there are some duplicated ones
	samples	<- make.unique(samples)
	# Add warnings if there are any NAs
	if (any(is.na(X)))
		warning("One or more samples contain NAs. Consider replacing these with e.g. Ct=40 now.")
	# Add some sort of phenoData
	if (missing(sample.info)) {
		pdata <- data.frame(sample=1:length(samples), row.names=samples)
		sample.info <- new("AnnotatedDataFrame", data = pdata,
            varMetadata = data.frame(labelDescription = "Sample numbering", row.names = "Sample names"))
	}
	# Define the 'history'
	X.hist	<- data.frame(history=capture.output(match.call(readCtData)), stringsAsFactors=FALSE)
	# Create the qPCRset object
    out	<- new("qPCRset", exprs=X, phenoData=sample.info, featureData=featData, featureCategory=X.cat, flag=X.flags, CtHistory=X.hist)
	# Return the object	
	out
}

.readCtPlain	<- 
function(readfile=readfile, header=header, n.features=n.features, n.data=n.data, i=i, ...)
{
	# A check for file dimensions. Read a single file.
	sample	<- read.delim(file=readfile, header=header, ...)
	nspots	<- nrow(sample)
	if (nspots != n.features*n.data[1])
		warning(paste(n.features, "gene names (rows) expected, got", nspots))
	# Read in the required file
	out	<- read.delim(file=readfile, header=header, colClasses="character", nrows=nspots*n.data[i], ...)
	# Return
	out
} # .readCtPlain

.readCtSDS	<- 
function(readfile=readfile, n.data=n.data, i=i, nspots=nspots, ...)
{
	# Scan through beginning of file, max 100 lines
	file.header <- readLines(con=readfile, n=100)
	n.header	<- grep("^#", file.header)
	if (length(n.header)==0) 
		n.header	<- 0
	# Read data, skip the required lines
	out	<- read.delim(file=readfile, header=FALSE, colClasses="character", nrows=nspots*n.data[i], skip=n.header, strip.white=TRUE, ...)
	# Return
	out
} # .readCtSDS

.readCtLightCycler	<- 
function(readfile=readfile, n.data=n.data, i=i, nspots=nspots, ...)
{
	# Read data, skip the required lines
	out	<- read.delim(file=readfile, header=TRUE, as.is=TRUE, nrows=nspots*n.data[i], skip=1, strip.white=TRUE, ...)
	# Return
	out
} # .readCtLightCycler

.readCtCFX	<- 
function(readfile=readfile, n.data=n.data, i=i, nspots=nspots, ...)
{
	# Read data, skip the required lines
	out	<- read.delim(file=readfile, header=TRUE, as.is=TRUE, nrows=nspots*n.data[i], strip.white=TRUE, dec = ",", ...)
	# Return
	out
} # .readCtCFX

.readCtOpenArray	<- 
function(readfile=readfile, n.data=n.data, i=i, nspots=nspots, ...)
{
	# Read data
	out	<- read.csv(file=readfile, header=TRUE, as.is=TRUE, nrows=nspots*n.data[i], strip.white=TRUE, ...)
	# Regard those marked as outliers as "Unreliable"
	out$ThroughHole.Outlier[out$ThroughHole.Outlier=="False"]	<- "OK"
	out$ThroughHole.Outlier[out$ThroughHole.Outlier=="True"]	<- "Unreliable"	
	# Return
	out
} # .readCtOpenArray

.readCtBioMark	<- 
function(readfile=readfile, n.data=n.data, i=i, nspots=nspots, ...)
{
	# Scan through beginning of file, max 100 lines
	file.header <- readLines(con=readfile, n=100)
	n.header	<- grep("^ID", file.header)-1
	if (length(n.header)==0) 
		n.header	<- 0
	# Read data, skip the required lines
	out	<- read.csv(file=readfile, header=TRUE, as.is=TRUE, nrows=nspots*n.data[i], skip=n.header, strip.white=TRUE, ...)
	# Convert the calls into flags
	out$Call[out$Call=="Pass"] <- "OK"
	out$Call[out$Call=="Fail"] <- "Undetermined"	
	# Return
	out
} # .readCtBioMark
