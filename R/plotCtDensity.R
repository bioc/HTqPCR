plotCtDensity <-
function(q,
	cards	= TRUE,
	xlab	= "Ct",
	ylab	= "Density",
	col,
	main	= NULL,
	legend	= TRUE,
	lwd	= 2,
	...) 
{
	# Get the data
	if (class(q)=="matrix") {
		data	<- q[,cards]	
	} else if (class(q)=="qPCRset") {
		data	<- exprs(q)[,cards]
	} else {
		stop("Data is of wrong format, only qPCRset and matrices are supported.\n")
	}
	# Calculate density; collect x and y values
	dens	<- apply(data, 2, density, na.rm=TRUE)
	x.pos <- do.call(cbind, lapply(dens, function(d) d$x))
    y.pos <- do.call(cbind, lapply(dens, function(d) d$y))
    # Plotting
    if (missing(col))
    	col	<- colorRampPalette(brewer.pal(11, "Spectral"))(ncol(data))
    matplot(x.pos, y.pos, xlab=xlab, ylab=ylab, col=col, main=main, type="l", lwd=lwd, ...)
    # Plot legend if required
    if (legend) 
    	legend(min(x.pos), max(y.pos), legend=sampleNames(q)[cards], col=col, lty=1:ncol(data), lwd=lwd, bty="n", ...)
}

