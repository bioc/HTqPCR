plotCtVariation <-
function(q,
	cards	= TRUE,
	variation	= "var",
	type	= "summary",
	sample.reps,
	feature.reps,
	log	= FALSE,
	add.featurenames	= FALSE,
	ylab,
	n.col,
	...)
{
	# If plotting parameters are changed
	old.par <- par(no.readonly=TRUE)
	on.exit(par(old.par))
	# Get the data
	if (class(q)=="matrix") {
		data <- q[,cards]	
	} else if (class(q)=="qPCRset") {
		data <- exprs(q)[,cards,drop=FALSE]
	} else {
		stop("Data is of wrong format, only qPCRset and matrices are supported.\n")
	}
	# Which type of variation to calculate
	switch(variation,
		var 	= {	var.fun = var;
				if (missing(ylab)) ylab="Variation"},
		sd	= {	var.fun = sd;
				if (missing(ylab)) ylab="Standard deviation" },
		{stop("Only 'var' and 'sd' are supported types of variation.\n")}
	) # switch
	# Replicate based on samples or features
	if (missing(sample.reps)) {
		if (missing(feature.reps)) {
			replicates <- featureNames(q)
		} else {
			replicates <- feature.reps
		}
	} else {
		replicates <- sample.reps
		data <- t(data)
	}
	# Calculate the variation
	ct.var <- aggregate(data, list(replicates), var.fun)
	ct.mean <- aggregate(data, list(replicates), mean)
	# Using the log10 instead?
	if (log) {
		ct.var[,2:ncol(ct.var)] <- log10(ct.var[,2:ncol(ct.var)])
		if (ylab %in% c("Variation", "Standard deviation"))
			ylab <- paste("log10", ylab)
	}
	# Plots the selected type of variation
	if (type=="summary") {
		boxplot(ct.var[,-1], ylab=ylab, ...)
	} else if (type=="detail") {
		# Set plotting parameters
 		if (missing(n.col))
			n.col <- min(ncol(data), 3)
		n.row <- ceiling(ncol(data)/n.col)	
		par(mfrow=c(n.row, n.col), mar=c(3,3,2,1), mgp=c(1.75, 0.75,0))
		# Plot the selected variation
		for (i in 1:(ncol(ct.var)-1)) {
			plot(ct.mean[,1+i], ct.var[,i+1], main=colnames(data)[i], ylab=ylab, xlab="Ct value", ...)
			# Add feature names if required
			if (add.featurenames)
				text(ct.mean[,i+1], ct.var[,i+1], ct.mean[,1], pch=NULL, ...)
		} # for i	
	} # if type
	# Return the results if desired
	invisible(list(Var=ct.var, Mean=ct.mean))	
}

