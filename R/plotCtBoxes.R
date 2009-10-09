plotCtBoxes <-
function(q,
	cards	= TRUE,
	xlab	= "",
	col,
	main	= NULL,
	names,
	stratify	= "type",
	mar	= c(7,4,3,1),
	...) 
{
	# Get the data
	if (class(q)=="matrix") {
		data	<- q[,cards]	
	} else if (class(q)=="qPCRset") {
		data	<- exprs(q)[,cards]
	} else {
		stop("Data is of wrong format, only qPCRset and matrices are supported.\n")
	}
    # Plotting parameters
    if (missing(col))
    	col	<- colorRampPalette(brewer.pal(11, "Spectral"))(ncol(data))
    if (missing(names))
    	names	<- sampleNames(q)[cards]
    par(mar=mar, ...)
    # Plotting stritified by spot type, gene class or simple version
    if (!is.null(stratify)) {
    	groups	<- list(dat=c(data), sample=rep(sampleNames(q)[cards], each=nrow(data)))
    	# Divide the data
 	    if (stratify=="type") {
 	    	groups[["strat"]]	<- rep(featureType(q), ncol(data))
  			n	<- length(levels(featureType(q)))
  			lev	<- levels(featureType(q))
  		} else if (stratify=="class") {
 	    	groups[["strat"]]	<- rep(featureClass(q), ncol(data))
  			n	<- length(levels(featureClass(q)))
  			lev	<- levels(featureClass(q))
		} else {
			stop(paste("Plot type \'", stratify, "\'isn't implemented\n"))
		}
		# The actual plotting
	    boxplot(dat~strat+sample, data=groups, col=rep(col, each=n), xlab=xlab, main=main, names=NA, ...)
	    mtext(1, at=seq(length(lev)-(length(lev)-1)/2, n*ncol(data), n), text=sampleNames(q)[cards], line=0.5)
	    mtext(1, at=1:(n*ncol(data)), text=gsub(" ", "\n", rep(lev), ncol(data)), line=2, las=2)
	} else {
	    boxplot(data, col=col, xlab=xlab, main=main, names=names, ...)
	}
}

