plotCtCategory <-
function(q,
	cards	= TRUE,
	by.feature	= FALSE,
	stratify,
	col,
	xlim,
	main,
	...) 
{
	# Get the data
	data	<- featureCategory(q)[,cards,drop=FALSE]
	colnames(data)	<- sampleNames(q)[cards]
	# Count categories. NB: Not every category always present.
	cats	<- sort(unique(unlist(data)))
	temp	<- array(0, c(length(cats), ncol(data)), list(cats, colnames(data)))
	if (length(temp)==1) {
		count <- matrix(nrow(data), dimnames=list(cats, colnames(data)))
	} else {
		count <- sapply(1:ncol(data), function(x) {xx <- table(data[,x]); temp[names(xx),x] <- xx; temp[,x]})
	}
	colnames(count)	<- colnames(data)
	# Various plotting parameters
	if (missing(col))
		col	<- c("#66C2A5", "#9E0142", "#FEE08B", "grey", "lightblue", "orange")
	if (missing(main))
		main	<- "Feature categories"
	# Adjust colours so they're consistent even if not all categories are represented
	all.cat.levels	<- sort(unique(unlist(featureCategory(q))))
	col	<- col[1:length(all.cat.levels)][all.cat.levels %in% cats] 
	# If features rather than samples is to be plotted
	if (by.feature) {
		data.num	<- apply(data, 2, function(x) as.numeric(as.factor(x)))
		heatmap(data.num, col=col, scale="none", labCol=colnames(data), labRow=featureNames(q), ...)
	} else if (!missing(stratify)) {
		# Stratify if required
		if (stratify=="type") {
			strat	<- featureType(q)
		} else if (stratify=="class") {
			strat	<- featureClass(q)
		} else {
			stop(paste("Plot type \'", stratify, "\'isn't implemented\n"))
		}
		# Divide data into stratified parts
		l.strat	<- length(unique(strat))
		s.data	<- list()
		for (i in 1:ncol(data)) {
			s.list	<- split(data[,i], strat)
			temp	<- array(0, c(length(cats), l.strat), list(cats, levels(strat)))
			count <- sapply(1:l.strat, function(x) {xx <- table(s.list[[x]]); temp[names(xx),x] <- xx; temp[,x]})
			colnames(count)	<- names(s.list)
			s.data[[i]]	<- count
		}
		all.data	<- do.call("cbind", s.data)
		# Misc parameters
		if (missing(xlim)) 
			xlim	<- c(0,ncol(all.data)*1.6)
		space	<- rep(c(0.3,rep(0.1, l.strat-1)), ncol(data))
		# Stratified plot
		barplot(all.data, col=col, las=2, main=main, legend.text=cats, space=space, xlim=xlim, ...)
		mtext(text=colnames(data), side=1, line=0, at=seq(0.3+1.1*l.strat/2, ncol(data)*(0.2+l.strat*1.1), l.strat*1.1+0.2))
	} else {
		if (missing(xlim))
			xlim	<- c(0, ncol(data)*1.6)
		# Standard plot
		barplot(count, col=col, las=2, legend.text=cats, xlim=xlim, main=main, ...)
	}
}

