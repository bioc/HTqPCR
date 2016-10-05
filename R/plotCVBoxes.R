plotCVBoxes <-
function(q,
	cards	= TRUE,
	xlab	= "",
	ylab	= "CV",
	col	= brewer.pal(5, "Spectral"),
	main	= NULL,
	stratify,
	...)
{
	# Extract the relevant data
	data <- exprs(q)[,cards]
	# Calculate coefficient of variation
	cv <- apply(data, 1, function(x) {sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE)})
 # Plotting stritified by spot type or otherwise
	if (missing(stratify)) {
 	boxplot(cv, xlab=xlab, ylab=ylab, main=main, col=col, ...)
	} else if (stratify=="type") {
 	groups <- list(cv=cv, type=featureType(q))
	 boxplot(cv~type, data=groups, xlab=xlab, ylab=ylab, main=main, names=levels(featureType(q)), col=col, ...)
	} else if (stratify=="class") {
 	groups <- list(cv=cv, type=featureClass(q))
	 boxplot(cv~type, data=groups, xlab=xlab, ylab=ylab, main=main, names=levels(featureClass(q)), col=col, ...)
 } else {
 	# Other stratifications?
 }
 # Return CVs if required
 invisible(cv)	
}

