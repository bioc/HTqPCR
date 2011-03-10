filterCtData <-
function(q,
	remove.type,
	remove.name,
	remove.class, 
	remove.category,
	n.category	= 3,
	remove.IQR,
	verbose	= TRUE)
{
	# Check that the necessary parameters are there
	if (missing(remove.type) & missing(remove.name) & missing(remove.class) & missing(remove.category) & missing(remove.IQR))
		stop("At least one \'remove\' parameter must be specified.\n")
	# Filter based on feature type
	if (!missing(remove.type)) {
		# Identifying + removing required features
		index	<- featureType(q) %in% remove.type
		q	<- q[!index,]
		# Print info if required
		if (verbose) 
			cat(paste("Removed ", sum(index), " \'", paste(remove.type, collapse="/"), "\' features based on featureType(q).\n", sep=""))
	}
	# Filter based on feature name
	if (!missing(remove.name)) {
		# Identifying + removing required features
		index	<- featureNames(q) %in% remove.name
		q	<- q[!index,]
		# Print info if required
		if (verbose) 
			cat(paste("Removed ", sum(index), " \'", paste(remove.name, collapse="/"), "\' features based on featureNames(q).\n", sep=""))
	}
	# Filter based on feature class
	if (!missing(remove.class)) {
		# Identifying + removing required features
		index	<- featureClass(q) %in% remove.class
		q	<- q[!index,]
		# Print info if required
		if (verbose)
			cat(paste("Removed ", sum(index), " \'", paste(remove.class, collapse="/"), "\' features based on featureClass(q).\n", sep=""))
	}
	# Filter based on feature category
	if (!missing(remove.category)) {
		# Identifying + removing required features
		cat.count	<- apply(featureCategory(q), 1, function(x) sum(x %in% remove.category))
		index	<- cat.count > n.category
		q	<- q[!index,]
		# Print info if required
		if (verbose)
			cat(paste("Removed ", sum(index), " features with >", n.category, " \'", paste(remove.category, collapse="/"), "\' based on featureCategory(q).\n", sep=""))
	}
	# Filter based on un-variable features
	if (!missing(remove.IQR)) {
		# Identifying + removing required features
		data	<- exprs(q)
#		data[data>na.IQR]	<- NA
		iqr	<- apply(data, 1, IQR, na.rm=TRUE)
		index	<- iqr < remove.IQR
		q	<- q[!index,]
		# Print info if required
		if (verbose)
			cat(paste("Removed ", sum(index, na.rm=TRUE), " features with IQR <", remove.IQR, " based on exprs(q).\n", sep=""))
	}	
	# Add to the history of the object
	if (nrow(getCtHistory(q))==0)
		q@history	<- data.frame(history="Manually created qPCRset object.", stringsAsFactors=FALSE)
	q@history	<- rbind(q@history, capture.output(match.call(filterCtData)))
	# Return the filtered object
	q	
}

