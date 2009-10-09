plotCtScatter <-
function(q,
	cards	= c(1,2),
	col	= "class",
	pch	= 20,
	diag	= FALSE,
	cor	= TRUE,
	Ct.max	= 35,
	legend	= TRUE,
	...)
{
	# Get the data
	data	<- exprs(q)[,cards]
	# Prepare colours
	col.selection	<- brewer.pal(10, "Spectral")[c(1,8,9,5,10,3,7,2,4,6)]
	if (col=="class") {
		colours	<- col.selection[as.numeric(featureClass(q))]
		col.names	<- levels(featureClass(q))
	} else if (col=="type") {
		colours	<- col.selection[as.numeric(featureType(q))]
		col.names	<- levels(featureType(q)) 
	} else {
		colours	<- col	
	}
	# Do the plotting
	plot(data, pch=pch, col=colours, ...)
	# Add diagonal if required
	if (diag)
		abline(0, 1, col="grey", lwd=2)
	# Add some information about the correlation if required
	if (cor) {
		index	<- apply(data, 1, function(x) any(x>Ct.max))
		corr	<- cor(data)[1,2]	
		corr2	<- cor(data[!index,])[1,2]
		text(min(data), c(max(data), max(data)-2)-2, labels=paste(c("R^2", paste("R^2 (Ct<", Ct.max, ")", sep="")), format(c(corr, corr2), digits=3), sep=": "), pos=4)
	}
	# Add legends if required
	if (legend & col %in% c("class", "type")) {
		legend(mean(range(data)), min(data)+5+length(col.names), legend=col.names, col=col.selection[1:length(col.names)], pch=pch, bty="n")
	}
}

