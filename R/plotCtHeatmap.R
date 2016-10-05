plotCtHeatmap <-
function(q,
	main	= NULL,
	col,
	col.range,
	dist	= "pearson",
	zero.center,
	mar,
	gene.names,
	sample.names,
	...)
{
	# Extract the data
	if (class(q)=="matrix") {
		data <- q	
	} else if (class(q)=="qPCRset") {
		data <- exprs(q)
	} else {
		stop("Data is of wrong format, only qPCRsets and matrices are supported.\n")
	}
	# Replace row/column names if required
	if (!missing(sample.names)) {
		if (length(sample.names)==1)
			sample.names <- rep(sample.names, ncol(data))
		colnames(data) <- sample.names
	}
	if (!missing(gene.names)) {
		if (length(gene.names)==1)
			gene.names <- rep(gene.names, nrow(data))
		rownames(data) <- gene.names
	}
	# Remove NAs, if there are any left
	na.index <- apply(data, 1, function(x) any(is.na(x)))
	data <- data[!na.index,]
	# Remove identical rows (e.g. all Ct=40)
	index <- apply(data, 1, function(x) all(x==x[1]))
	data <- data[!index,]
# TO DO: Adjust the colours - only for Ct?
	# Set colours
	if (missing(col)) {
		if (any(data<0)) {
			col <- colorRampPalette(brewer.pal(11, "RdYlGn"))(21)
		} else {
			col <- colorRampPalette(brewer.pal(11, "Spectral"))(20)
		}
	}
	# Set colour range
	if (missing(zero.center)) {
		if (any(data<0)) {
			zero.center <- TRUE	
		} else {
			zero.center <- FALSE	
		}	
	}
	if (missing(col.range)) {
		if (zero.center) {
			m <- max(abs(data))
			breaks <- seq(-m, m, length.out=length(col)+1)
		} else {
			breaks <- seq(min(data), max(data), length.out=length(col)+1)
		}
	} else {
		breaks <- seq(col.range[1], col.range[2], length.out=length(col)+1)
	}
	# Set margin sizes
	if (missing(mar)) {
		mar <- c(0.4*max(nchar(colnames(data)))+3, 0.4*max(nchar(rownames(data)))+1)
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
	heatmap.2(data, trace="none", density.info="none", main=main, col=col, distfun=d, breaks=breaks, mar=mar, ...)
}

