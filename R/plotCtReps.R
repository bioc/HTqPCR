plotCtReps <-
function(q,
	card	= 1,
	percent	= 20,
	verbose	= TRUE,
	col	= 1,
	...)
{
	# Get the data
	if (class(q)=="qPCRset") {
		data <- exprs(q)
	} else if (class(q)=="matrix") {
		data <- q
	} else {
		stop("Data is of wrong format, only qPCRsets and matrices are supported.\n")
	}
	genes <- featureNames(q)
	data <- data[,card]
	# Split the data by replicated genes
	split.data <- split(data, genes)
	split.data <- split.data[listLen(split.data)==2]
	# Prepare the plot
	r <- range(split.data, na.rm=TRUE)
	plot(x=r, y=r, type="n", xlab="", ylab="", ...)
	abline(0, 1, col="darkgrey", lty=2, lwd=2)
	# Plot the rep points
	out <- data.frame(rep1=vector(), rep2=vector())
	for (s in 1:length(split.data)) {
		# The "default" points
		dat <- split.data[[s]]
		names(split.data)[s]
		points(dat[1], dat[2], pch=19, col=col, ...)
		# Those with difference larger than specifyied by percent
		if (abs(dat[1]-dat[2]) > percent/100*abs(mean(dat))) {
			# Plot them
			gene <- names(split.data)[s]
			points(dat[1], dat[2], pch=1, cex=2)
			text(dat[1], dat[2], gene, pos=1, col=2)
			# Add to (potential) output
			out[gene,] <- dat
		}	
	}
	# Print details if required
	if (verbose) {
		cat(paste("Replicates differing > ", percent, "% on card ", card, ":\n", sep=""))
		print(out)
	}
	# Return differing reps if required
	invisible(out)
}

