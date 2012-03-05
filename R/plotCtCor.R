plotCtCor <-
function(q,
	col,
	col.range	= c(0,1),
	main,
	mar,
	...)
{
	# Get the data
	if (class(q)=="qPCRset") {
		data	<- exprs(q)
	} else if (class(q)=="matrix") {
		data	<- q
	} else {
		stop("Data is of wrong format, only qPCRsets and matrices are supported.\n")
	}
	# Calculate 1 minus the correlation
	x	<- 1-cor(data)
	# Set the colours
	if (missing(col)) {
		col	<- colorRampPalette(brewer.pal(11, "Spectral"))(20)
	}
	# Set the breaks
	b	<- seq(col.range[1], col.range[2], length.out=length(col)+1)
	# Setting different plotting parameters
	if (missing(main)) {
		main	<- "Correlation between Ct values"
	}
	if (missing(mar)) {
		max	<- max(nchar(colnames(data)))
		mar	<- 0.4*max+5
		mar	<- c(mar, mar)
	}
	# Plot
	heatmap.2(x, col=col, breaks=b, scale="none", dendrogram="row", trace="none", main=main, density.info="none", mar=mar, ...)
}
