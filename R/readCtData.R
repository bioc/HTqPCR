readCtData <-
function(files,
	path	= NULL,
	n.features	= 384,
#	format	= "none",
	flag	= 4,
	feature	= 6,
	type	= 7,
	position	= 3,
	Ct	= 8,
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
	if (is.null(Ct) & !SDS) 
		stop("Either input must be specified as SDS, or a column containing CT values given")
	# Various parameters
	if (length(files)!=length(n.data) & length(n.data)!=1)
		stop("n.data must either be a single integer, or same length as number of files")
	if (length(n.data)==1) 
		n.data	<- rep(n.data, length(files))
	nsamples	<- sum(n.data)
	ncum	<- cumsum(n.data)
	s.names	<- NULL
	# If all files are "easy" format
	if (!SDS) {
		# Read a single file, to get dimensions etc.
		readfile	<- ifelse(is.null(path), files[1], file.path(path, files[1]))
		sample	<- read.delim(file=readfile, header=header, ...)
		nspots	<- nrow(sample)
		if (nspots != n.features*n.data[1])
			warning(paste(n.features, "gene names (rows) expected, got", nspots))
	}
	nspots	<- n.features
	# Initializing data required for qPCRset object.
	X <- matrix(0,nspots,nsamples)
	X.flags	<- as.data.frame(X)
	X.cat <- data.frame(matrix("OK", ncol=nsamples, nrow=nspots), stringsAsFactors=FALSE)
	# Read through all files (assume features are sorted based on position)
	for (i in 1:length(files)) {
		# Which columns of the qPCRset object are to be filled out
		if (i==1) {
			cols	<- 1:ncum[i]
		} else {
			cols	<- (ncum[i-1]+1):ncum[i]
		}
		# The file to read
		readfile	<- ifelse(is.null(path), files[i], file.path(path, files[i]))
		# Format dependent reading of files
		if (!SDS) {
			# Simple - just read data
#			if (format=="stratagene") {
#				sample	<- read.delim(file=readfile, header=TRUE, colClasses="character", nrows=nspots*n.data[i], strip.white=TRUE, check.names=FALSE, ...)
#			} else {
				sample	<- read.delim(file=readfile, header=header, colClasses="character", nrows=nspots*n.data[i], ...)
#			}	
		} else if (SDS) {
			# Scan through beginning of file, max 100 lines
			file.header <- readLines(con=readfile, n=100)
			n.header	<- grep("^#", file.header)
			if (length(n.header)==0) 
				n.header	<- 0
			# Read data, skip the required lines
			sample	<- read.delim(file=readfile, header=FALSE, colClasses="character", nrows=nspots*n.data[i], skip=n.header, strip.white=TRUE, ...)
		}
		# Assign Ct values to object
		data	<- matrix(sample[,Ct], ncol=n.data[i])
		undeter	<- apply(data, 2, function(x) x %in% c("Undetermined", "No Ct"))
  		X.cat[,cols][undeter] <- "Undetermined"
		if (is.null(na.value)) {
 			data[data=="Undetermined"]	<- NA
 		} else {
 			data[data %in% c("Undetermined", "No Ct", "999")] <- na.value
 		}
		X[,cols]	<- apply(data, 2, function(x) as.numeric(as.character(x)))
		if (!is.null(flag)) {
			flags	<-  matrix(sample[,flag], ncol=n.data[i])
			flags[flags=="-"]	<- "Failed"
			flags[flags=="+"]	<- "Passed"
			X.flags[,cols]	<- flags
		} else {
			X.flags[,cols]	<- "Passed"
		}
		# Get sample names if present
		s.names	<- c(s.names, unique(sample[,2]))
		# Assemble the feature data info to the qPCRset
		if (i==1) {
			# Add default values if any are missing from file
			featPos	<- paste("feature", 1:nspots, sep="")
			if (!is.null(position))
				featPos	<- as.character(sample[1:nspots,position])
			featType	<- factor(rep("Target", nspots))
			if (!is.null(type))		
				featType	<- sample[1:nspots,type]
			featName	<- paste("feature", 1:nspots, sep="")
			if (!is.null(feature)) 
				featName	<- as.character(sample[1:nspots,feature])
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

