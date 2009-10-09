filterCategory <-
function(q,
	na.categories	= c("Unreliable", "Undetermined"))
{
	index	<- apply(featureCategory(q), 2, "%in%", na.categories)
	exprs(q)[index] <- NA
	q		
}

