plotCtPCA <-
function(q,
	s.names,
	f.names,
	scale	= TRUE,
	features	= TRUE,
	col,
	cex	= c(1,1)) 
{
	# Extract the data
	if (class(q)=="qPCRset") {
		data	<- exprs(q)
	} else if (class(q)=="matrix") {
		data	<- q
	} else {
		stop("Data is of wrong format, only qPCRsets and matrices are supported.\n")
	}
	# Setting colours
	if (missing(col)) {
		col	<- colorRampPalette(brewer.pal(11, "Spectral"))(ncol(data))
	}
	# Remove columns that are all NAs
	data	<- data[,!apply(data, 2, function(x) all(is.na(x)))]
	# The actual PCA
	pca	<- prcomp(data[!apply(data, 1, function(x) any(is.na(x))),], scale.=scale) 
	# Change dimension names
	if (missing(f.names)) {
		rownames(pca$x)	<- rep("*", nrow(pca$x))
	} else {
		rownames(pca$x)	<- f.names
	}
	if (!missing(s.names)) {
		rownames(pca$rotation) <- s.names
	}
	# Plot - either including features, or just for samples
	if (features) {
		biplot(pca, cex=cex) 
	} else {
		x	<- pca$rotation[,1]
		y	<- pca$rotation[,2]
		plot(x, y, pch=19, col=col, 
				xlim=range(x)+c(-abs(min(x))*.1, abs(max(x))*.1),  
				ylim=range(y)+c(-abs(min(y))*.1, abs(max(y))*.1))
		text(pca$rotation[,1], pca$rotation[,2], colnames(data), pos=1)
	}
}

