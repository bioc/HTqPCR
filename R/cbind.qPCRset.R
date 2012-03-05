cbind.qPCRset <-
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
		if (n.wells(qsets[[1]])!=n.wells(qsets[[i]])) {
			warning(paste(qset.names[i], "has a different number of features than", qset.names[1], "and is not included in the combined object."))
			next	
		}
		# Warn if not all "feature markers" are identical
		if (!all(featureNames(qsets[[i]])==featureNames(qsets[[1]])))
			warning(paste("Feature names in", qset.names[i], "are not identical to", qset.names[1]))
		if (!all(featureType(qsets[[i]])==featureType(qsets[[1]])))
			warning(paste("Feature types in", qset.names[i], "are not identical to", qset.names[1]))
		if (!all(featurePos(qsets[[i]])==featurePos(qsets[[1]])))
			warning(paste("Feature positions in", qset.names[i], "are not identical to", qset.names[1]))
		if (!all(featureClass(qsets[[i]])==featureClass(qsets[[1]])))
			warning(paste("Feature classes in", qset.names[i], "are not identical to", qset.names[1]))
		# Bind together the remaining info
		exprs(out)	<- cbind(exprs(out), exprs(qsets[[i]]))
		flag(out)	<- cbind(flag(out), flag(qsets[[i]]))
		featureCategory(out)	<- cbind(featureCategory(out), featureCategory(qsets[[i]]))
		phenoData(out)	<- AnnotatedDataFrame(rbind(pData(out), pData(qsets[[i]])))
		#sampleNames(out)	<- c(sampleNames(out), sampleNames(qsets[[i]]))
	} # for i in qsets
	# Update the history slot
	for (q in seq_along(qsets)) {
		if (nrow(getCtHistory(qsets[[q]]))==0)
			setCtHistory(qsets[[q]])	<- data.frame(history="Manually created qPCRset object.", stringsAsFactors=FALSE)
	}
	all.hist	<- sapply(qsets, getCtHistory)
	new.hist	<- paste(rep(qset.names, times=sapply(all.hist, length)), unlist(all.hist), sep=": ")
	new.hist	<- data.frame(history=new.hist, stringsAsFactors=FALSE)
	setCtHistory(out)	<- rbind(new.hist, capture.output(match.call(cbind)))

	# Return object
	out
}

