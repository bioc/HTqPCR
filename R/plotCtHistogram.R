plotCtHistogram <-
function(q,
	card	= 1,
	xlab	= "Ct",
	col,
	main,
	n	= 30,
	...) 
{
	# Get the data
	if (class(q)=="matrix") {
		data	<- q[,card]	
		if (missing(main))
			main	<- colnames(q)[card]
	} else if (class(q)=="qPCRset") {
		data	<- exprs(q)[,card]
		if (missing(main))
			main	<- sampleNames(q)[card]
	} else {
		stop("Data is of wrong format, only qPCRset and matrices are supported.\n")
	}
	# Setting some plotting parameters
    if (missing(col))
    	col	<- "#66C2A5"
	# Plot histogram
	hist(data, breaks=n, col=col, main=main, xlab=xlab, ...)
	
}

