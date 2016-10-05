plotCtPairs <-
function(q,
	cards	= TRUE,
	lower.panel	= panel.Ct.cor,
	upper.panel	= panel.Ct.scatter,
	Ct.max	= 35,
	col	= "type",
	pch	= 20,
	cex.cor	= 2,
	cex.pch	= 1,
	diag	= TRUE,
	...
)
{
	# Get the data
	data <- exprs(q)[,cards]
	# Prepare colours
	col.selection <- brewer.pal(10, "Spectral")[c(1,8,9,5,10,3,7,2,4,6)]
	if (length(col)==1 & col=="class") {
		feat <- as.factor(featureClass(q))
		if (length(feat)==0)
			stop("No featureClass available for object 'q'")
		colours <- col.selection[as.numeric(feat)]
		col.names <- levels(feat)
	} else if (length(col)==1 & col=="type") {
		feat <- as.factor(featureType(q))
		colours <- col.selection[as.numeric(feat)]
		col.names <- levels(feat)
	} else {
		colours <- col	
	}
	# The actual plotting
	suppressWarnings(pairs(data, lower.panel=lower.panel, upper.panel=upper.panel, point.col=colours, point.pch=pch, draw.diag=diag, Ct.max=Ct.max, text.cex=cex.cor, point.cex=cex.pch, ...))
}

