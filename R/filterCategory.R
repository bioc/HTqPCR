filterCategory <-
function(q,
	na.categories	= c("Unreliable", "Undetermined"))
{
	index	<- apply(featureCategory(q), 2, "%in%", na.categories)
	exprs(q)[index] <- NA
	# Add to the history of the object
	if (nrow(getCtHistory(q))==0)
		q@history	<- data.frame(history="Manually created qPCRset object.", stringsAsFactors=FALSE)
	q@history	<- rbind(q@history, capture.output(match.call(filterCategory)))
	q		
}

