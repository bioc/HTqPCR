rbind.qPCRset <-
function(..., deparse.level=1) {
	# Make list of all the qPCRsets to be combined
	qsets	<- list(...)
	n.qsets	<- length(qsets)
	# Get the names of the qPCRsets to be joined
	qset.names	<- as.list(sys.call(which=1))[-1]
	# If this seems to have not worked, name them consecutively
	if (length(qset.names)!=n.qsets)
		qset.names	<- paste("qPCRset", 1:n.qsets, sep="")
	# Prepare output object
	out <- qsets[[1]]
	if (n.qsets==1) 
		return(out)
	# Run through all qPCRsets
	for (i in 2:n.qsets) {
		# Skip if dimensions don't match
		if (n.samples(qsets[[1]])!=n.samples(qsets[[i]])) {
			warning(paste(qset.names[i], "has a different number of samples than", qset.names[1], "and is not included in the combined object."))
			next	
		}
		# Warn if not all sample names are identical
		if (!all(sampleNames(qsets[[i]])==sampleNames(qsets[[1]])))
			warning(paste("Sample names in", qset.names[i], "are not identical to", qset.names[1]))
#		# Bind together the remaining info - factors
#		featureType(out)	<- factor(c(as.vector(featureType(out)), as.vector(featureType(qsets[[i]]))))
#		featureClass(out)	<- factor(c(as.vector(featureClass(out)), as.vector(featureClass(qsets[[i]]))))
		# Bind together the remaining info - matrix
		exprs(out)	<- rbind(exprs(out), exprs(qsets[[i]]))
		# Bind together the remaining info - data frames
		flag(out)	<- as.data.frame(rbind(as.matrix(flag(out)), as.matrix(flag(qsets[[i]]))), stringsAsFactors=FALSE)
		featureCategory(out)	<- as.data.frame(rbind(as.matrix(featureCategory(out)), as.matrix(featureCategory(qsets[[i]]))), stringsAsFactors=FALSE)
		# Bind together the remaining info - vectors
		#featureNames(out)	<- c(featureNames(out), featureNames(qsets[[i]]))
		#featurePos(out)	<- c(featurePos(out), featurePos(qsets[[i]]))
		featureData(out)	<- AnnotatedDataFrame(rbind(fData(out), fData(qsets[[i]])))

	} # for i in qsets
	# Update the history slot
	for (q in seq_along(qsets)) {
		if (nrow(getCtHistory(qsets[[q]]))==0)
			setCtHistory(qsets[[q]])	<- data.frame(history="Manually created qPCRset object.", stringsAsFactors=FALSE)
	}
	all.hist	<- sapply(qsets, getCtHistory)
	new.hist	<- paste(rep(qset.names, times=sapply(all.hist, length)), unlist(all.hist), sep=": ")
	new.hist	<- data.frame(history=new.hist, stringsAsFactors=FALSE)
	setCtHistory(out)	<- rbind(new.hist, capture.output(match.call(rbind)))

	# Return object
	out
}
