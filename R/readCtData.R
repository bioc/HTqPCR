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
#	# If the files are from a particular manufacturer, set the corresponding parameters
#	if (format=="tlda") {
#		flag	= 4
#		feature	= 6
#		type	= 7
#		position= 3
#		Ct		= 8
#	} else if (format=="stratagene") {
#		flag	= "Final Call (dR)"
#		feature	= "Well Name"
#		type	= "Well Type"
#		position= "Well"
#		Ct		= "Ct (dR)"
#		header	= TRUE
#	}
	# If all files are "easy" format
	if (!SDS) {
		# Read a single file, to get dimensions etc.
		readfile	<- ifelse(is.null(path), files[1], file.path(path, files[1]))
		sample	<- read.delim(file=readfile, header=header)
		nspots	<- nrow(sample)
		if (nspots != n.features*n.data[1])
			stop(paste(n.features, "gene names (rows) expected, got", nspots))
	}
	nspots	<- n.features
	# Initializing qPCRset object.
	X <- matrix(0,nspots,nsamples)
	cat <- data.frame(matrix("OK", ncol=nsamples, nrow=nspots), stringsAsFactors=FALSE)
	out <- new("qPCRset", exprs=X, flag=as.data.frame(X), featureCategory=cat)
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
  		featureCategory(out)[,cols][undeter] <- "Undetermined"
		if (is.null(na.value)) {
 			data[data=="Undetermined"]	<- NA
 		} else {
 			data[data %in% c("Undetermined", "No Ct")] <- na.value
 		}
		out@exprs[,cols]	<- apply(data, 2, function(x) as.numeric(as.character(x)))
		if (!is.null(flag)) {
			flags	<-  matrix(sample[,flag], ncol=n.data[i])
			flags[flags=="-"]	<- "Failed"
			flags[flags=="+"]	<- "Passed"
			out@flag[,cols]	<- flags
		} else {
			out@flag[,cols]	<- "Passed"
		}
		# Get sample names if present
		s.names	<- c(s.names, unique(sample[,2]))
		# Add the additional info to the qPCRset
		if (i==1) {
			# Add default values if any are missing from file
			featPos	<- paste("feature", 1:nspots, sep="")
			if (!is.null(position))
				featPos	<- as.character(sample[,position])
			featType	<- factor(rep("Target", nspots))
			if (!is.null(type))		
				featType	<- sample[,type]
			featName	<- paste("feature", 1:nspots, sep="")
			if (!is.null(feature)) 
				featName	<- as.character(sample[1:nspots,feature])
			# Add the data
			featureNames(out)	<- featName
			featureType(out)	<- as.factor(featType)
			featurePos(out)	<- featPos
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
	sampleNames(out)	<- samples
	colnames(flag(out))	<- samples
	colnames(featureCategory(out))	<- samples
	# Return the object	
	out
}

