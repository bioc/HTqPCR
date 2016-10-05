clusterCt <-
function(q,
	main	= NULL,
	type	= "genes",
	dist	= "pearson",
	xlab	= "Cluster dendrogram",
	n.cluster,
	h.cluster,
	select.cluster	= FALSE,
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
	# Transpose if samples rather than genes are chosen
	if (type=="samples")
		data <- t(data)
	# Remove NAs (there should be any though...)
	na.index <- apply(data, 1, function(x) any(is.na(x)))
	data <- data[!na.index,]
	# Remove identical rows (e.g. all Ct=40)
	index <- apply(data, 1, function(x) all(x==x[1]))
	data <- data[!index,]
	# Calculate the distance
	if (dist=="pearson") {
		d <- as.dist(1-cor(t(data)))
	} else if (dist=="euclidean") {
		d <- dist(data, method="euclidean")
	} else {
		stop(paste("Distance method \'", dist, "\' is not implemented\n", sep=""))
	}
	# Prepare the tree
	tree <- hclust(d)
	# The actual plotting
	plot(tree, main=main, xlab=xlab, sub="", ...)
	# Preselected clusters
	blocks <- s.blocks <- NULL
	c.col <- brewer.pal(11, "Spectral")[c(1,3,5,9:11)]
	if (!missing(n.cluster)) {
		blocks <- rect.hclust(tree, k=n.cluster, border=c.col)
	} else if (!missing(h.cluster)) {
		blocks <- rect.hclust(tree, h=h.cluster, border=c.col)
	}
	# Manually selecting some blocks
	if (select.cluster) {
		s.blocks <- identify(tree)
	}
	invisible(c(blocks, s.blocks))
}

