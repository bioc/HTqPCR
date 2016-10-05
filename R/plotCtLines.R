plotCtLines <-
function(q,
	genes,
	groups,
	col	= brewer.pal(10, "Spectral"),
	xlab	= "Sample",
	ylab	= "Ct",
	legend	= TRUE,
	lwd	= 2,
	lty,
	pch,
	xlim,
	...)
{
	# Extract the genes of interest
	if (missing(genes))
		stop("Please select the genes to plot.")
	data <- exprs(q)[genes,]
	g.names <- rownames(data)
	# If groups are supplied, take the average of each feature within groups
	if (!missing(groups)) {
		temp <- aggregate(t(data), by=list(groups), FUN=mean)
		data2 <- t(temp[,-1])
		colnames(data2) <- temp[,1]
	} else {
		data2 <- data
	}	 	
	# Setting some plotting parameters
	if (missing(lty))
		lty <- 1:ncol(data2)	
	if (missing(pch))
		pch <- 1:nrow(data2)
	if (missing(xlim))
		xlim <- c(ifelse(legend,-1,1), ncol(data2))
	# The actual plotting
	matplot(x=1:ncol(data2), y=t(data2), type="l", col=col, lwd=lwd, lty=lty, xlab=xlab, ylab=ylab, ylim=c(min(data), max(data)), xlim=xlim, xaxt="n", ...)
	axis(1, at=1:ncol(data2), labels=colnames(data2))
	# Plot the individaul points if averaged
	if (!missing(groups)) {
		if (!is.factor(groups))
			groups <- as.factor(groups)
		matpoints(x=as.numeric(groups), y=t(data), col=col, pch=1:ncol(data))
	}
	# Plot legend, if required
	if (legend) {
		legend(x=xlim[1], y=max(data, na.rm=TRUE), legend=g.names, col=col, bty="n", lwd=lwd, lty=lty)
	}
}

