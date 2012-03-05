changeCtLayout <-
function(q, sample.order) {
	# Check if the sample order is a factor
	if (!is.factor(sample.order))
		sample.order	<- as.factor(sample.order)

	# Split Ct values
	X	<- split(as.data.frame(exprs(q)), sample.order)
	Xnew <- do.call("cbind", X)
	colnames(Xnew)	<- paste(rep(levels(sample.order), each=n.samples(q)), sampleNames(q), sep=":")
	# Split categories
	cats	<- split(as.data.frame(featureCategory(q)), sample.order)
	catsnew <- do.call("cbind", cats)
	# Split flags
	flags <- split(as.data.frame(flag(q)), sample.order)
	flagsnew <- do.call("cbind", flags)

	# Make new qPCRset
	out <- new("qPCRset", exprs=as.matrix(Xnew), flag=flagsnew, featureCategory=catsnew)
	# Add some extra info 
#	sampleNames(out)	<- paste(rep(levels(sample.order), each=n.samples(q)), sampleNames(q), sep=":")
	featureData(out)	<- featureData(q)[as.numeric(sample.order)==1]
	#featureNames(out)	<- featureNames(q)[as.numeric(sample.order)==1]
	#featureType(out)	<- featureType(q)[as.numeric(sample.order)==1]
	#featureClass(out)	<- featureClass(q)[as.numeric(sample.order)==1]
	phenoData(out)	<- new("AnnotatedDataFrame", data=data.frame(sample=1:ncol(out), row.names=sampleNames(out)), varMetadata=data.frame(labelDescription = "Sample numbering", row.names = "Sample names"))
	
	# Opdate "history" slot
	if (nrow(getCtHistory(out))==0)
		setCtHistory(out)	<- data.frame(history="Manually created qPCRset object.", stringsAsFactors=FALSE)
	setCtHistory(out)	<- rbind(getCtHistory(out), capture.output(match.call(changeCtLayout)))

	# Return the new object
	out
}

