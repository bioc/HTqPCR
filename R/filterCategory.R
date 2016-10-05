filterCategory <-
function(q,
	na.categories	= c("Unreliable", "Undetermined"))
{
	index <- apply(featureCategory(q), 2, "%in%", na.categories)
	exprs(q)[index] <- NA
	# Add to the history of the object
	if (nrow(getCtHistory(q))==0)
		setCtHistory(q) <- data.frame(history="Manually created qPCRset object.", stringsAsFactors=FALSE)
	setCtHistory(q) <- rbind(getCtHistory(q), capture.output(match.call(filterCategory)))
	q		
}

