changeCtLayout <-
function(q, sample.order) {
	# Check if the sample order is a factor
	if (!is.factor(sample.order))
		sample.order	<- as.factor(sample.order)

	# Split Ct values
	X	<- split(as.data.frame(exprs(q)), sample.order)
	Xnew <- do.call("cbind", X)
	# Split categories
	cats	<- split(as.data.frame(featureCategory(q)), sample.order)
	catsnew <- do.call("cbind", cats)
	# Split flags
	flags <- split(as.data.frame(flag(q)), sample.order)
	flagsnew <- do.call("cbind", flags)

	# Make new qPCRset
	out <- new("qPCRset", exprs=as.matrix(Xnew), flag=flagsnew, featureCategory=catsnew)
	# Add some extra info (not in matrix or data.frame form)
	sampleNames(out)	<- paste(rep(levels(sample.order), each=n.samples(q)), sampleNames(q), sep=":")
	featureNames(out)	<- featureNames(q)[as.numeric(sample.order)==1]
	featureType(out)	<- featureType(q)[as.numeric(sample.order)==1]
	featureClass(out)	<- featureClass(q)[as.numeric(sample.order)==1]
	
	# Opdate "history" slot
	if (nrow(getCtHistory(q))==0)
		q@history	<- data.frame(history="Manually created qPCRset object.", stringsAsFactors=FALSE)
	out@history	<- rbind(q@history, capture.output(match.call(changeCtLayout)))

	# Return the new object
	out
}

