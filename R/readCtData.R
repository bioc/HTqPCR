readCtData <-
function(files,
	path	= NULL,
	n.features	= 384,
	flag	= 4,
	feature	= 6,
	type	= 7,
	position	= 3,
	Ct	= 8,
	header	= FALSE,
	SDS	= FALSE,
	samples,
	na.value	= 40,
	...)
{
	# Initial checks
	if(missing(files))
		stop("No input files specified")
	# Various parameters
	nsamples	<- length(files)
	# If all files are "easy" format
	if (!SDS) {
		# Read a single file, to get dimensions etc.
		sample	<- read.delim(file=paste(path, files[1], sep="/"), header=header)
		nspots	<- nrow(sample)
		if (nspots != n.features)
			stop(paste(n.features, "gene names (rows) expected, got", nspots))
	}
	nspots	<- n.features
	# Initializing qPCRset object.
	X <- matrix(0,nspots,nsamples)
	cat <- data.frame(matrix("OK", ncol=nsamples, nrow=nspots), stringsAsFactors=FALSE)
	out <- new("qPCRset", exprs=X, flag=as.data.frame(X), featureCategory=cat)
	# Read through all files (assume features are sorted based on position)
	for (i in 1:nsamples) {
		# Format dependent reading of files
		if (!SDS) {
			# Simple - just read data
			sample	<- read.delim(file=paste(path, files[i], sep="/"), header=header, colClasses="character", nrows=nspots)
		} else if (SDS) {
			# Scan through beginning of file, max 100 lines
			file.header <- readLines(con=paste(path, files[i], sep="/"), n=100)
			n.header	<- grep("^#", file.header)
			if (length(n.header)==0) 
				n.header	<- 0
			# Read data, skip the required lines
			sample	<- read.delim(file=paste(path, files[i], sep="/"), header=FALSE, colClasses="character", nrows=nspots, skip=n.header)
		}
		# Assign Ct values to object
		data	<- sample[,Ct]
  		featureCategory(out)[data=="Undetermined",i] <- "Undetermined"
		if (is.null(na.value)) {
 			data[data=="Undetermined"]	<- NA
 		} else {
 			data[data=="Undetermined"]	<- na.value
 		}
		out@exprs[,i]	<- as.numeric(as.character(data))
		if (!is.null(flag)) {
			out@flag[,i]	<- sample[,flag]
		} else {
			out@flag[,i]	<- rep("Passed", nspots)
		}
		# Add the additional info to the qPCRset
		if (i==1) {
			# Add default values if any are missing from file
			featPos	<- paste("feature", 1:nspots, sep="")
			if (!is.null(position))
				featPos	<- as.character(sample[,position])
			featType	<- factor(rep("Target", nspots))
			if (!is.null(type))		
				featType	<- sample[,type]
			# Add the data
			featureNames(out)	<- as.character(sample[,feature])
			featureType(out)	<- as.factor(featType)
			featurePos(out)	<- featPos
		}
	}
	# Provide sample names
	if (missing(samples))
		samples	<- gsub("(.+)\\..+", "\\1", files)
	sampleNames(out)	<- samples
	# Return the object	
	out
}

