plotCtPCA <-
function(q,
	s.names,
	g.names) 
{
	# Extract the data
	if (class(q)=="qPCRset") {
		data	<- exprs(q)
	} else if (class(q)=="matrix") {
		data	<- q
	} else {
		stop("Data is of wrong format, only qPCRsets and matrices are supported.\n")
	}
	# The actual PCA
	pca	<- prcomp(data[!apply(data, 1, function(x) any(is.na(x))),]) 
	# Change dimension names
	if (missing(g.names)) {
		rownames(pca$x)	<- rep("*", nrow(pca$x))
	} else {
		rownames(pca$x)	<- g.names
	}
	if (!missing(s.names)) {
		rownames(pca$rotation) <- s.names
	}
	# Plot
	biplot(pca)
}

