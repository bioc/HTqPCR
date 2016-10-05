heatmapSig <-
function(qDE,
	comparison	= "all",
	col,
	zero.center	= TRUE,
	mar,
	dist	= "pearson",
	...)
{
	# Get the data. Remove the test summary
	if (is.numeric(comparison))
		qDE <- qDE[comparison]
	qDE <- qDE[names(qDE)!="Summary"]	
 	data <- sapply(qDE, function(x) x$ddCt)
 	colnames(data) <- names(qDE)
 	# Prepare colour
 	if (missing(col))
 		col <- colorRampPalette(brewer.pal(9, "RdYlGn"))(21)
 	if (zero.center) {
 		b <- seq(-max(abs(data)), max(abs(data)), length.out=length(col)+1)
 	} else {
 		b <- seq(min(data), max(data), length.out=length(col)+1)
 	}
	# Set margin sizes
	if (missing(mar)) {
		mar <- c(0.4*max(nchar(colnames(data)))+5, 0.4*max(nchar(as.character(qDE[[1]][,"genes"])))+1)
	}
	# Set clustering method
	if (dist=="pearson") {
		d <- function(x) as.dist(1-cor(t(x)))
	} else if (dist=="euclidean") {
		d <- function(x) dist(x, method="euclidean")
	} else {
		stop(paste("Distance method \'", dist, "\' is not implemented\n", sep=""))
	}
	# The actual plotting
	heatmap.2(data, trace="none", density.info="none", col=col, distfun=d, breaks=b, mar=mar, ...)
}

