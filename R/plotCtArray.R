plotCtArray <-
function(q,
	plot	= "Ct",
	main,
	col,
	col.range,
	na.col	= "grey",
	na.value	= 40,
	chamber.size,
	...)
{
	# Define the plot layout
	old.par <- par(no.readonly=TRUE)
	on.exit(par(old.par))
	layout(rbind(2,1), widths=5, heights=c(5,1), respect=TRUE)
	# Some general "size" parameters
	ncol <- n.samples(q)
	nrow <- n.wells(q)
	if (missing(chamber.size)) {
		chamber.size <- 70/nrow
	}	
	# Ct or other values on a continuous scale (might add others later)
	if (plot %in% c("Ct")) {
		# Extract the Ct values
		data <- exprs(q)
		data[data==na.value] <- NA
		# Define the colours for the card
		if (missing(col)) {
			col <- rev(colorRampPalette(brewer.pal(9, "YlGnBu"))(20))
		} 		
		if (missing(col.range)) {
			breaks <- cut(data, breaks=length(col))
		} else {
			breaks <- cut(data, seq(col.range[1], col.range[2], length.out=length(col)+1), labels=FALSE)
			breaks[data<col.range[1]] <- 1
			breaks[data>col.range[2]] <- length(col)
		}
		values <- col[breaks]
		values[is.na(values)] <- na.col
		# Plot the legend
		par(mar=c(0,1,0,1))
		l <- length(col)
		data <- data[!is.na(data)]
		at <- seq(0,1,length.out=5)
		plot(0:1,0:1, type="n", xaxt="n", yaxt="n", ylab="", bty="n", xlab="")
		rect(seq(0,1,length.out=l+1)[-(l+1)],0.6,seq(0,1,length.out=l+1)[-1], 0.75, col=col)
		if (!missing(col.range)) {
			lab <- format(quantile(col.range, at), digits=2)
			lab <- paste(c("<", rep("", length(at)-2), ">"), lab, sep="")
		} else {
			lab <- format(quantile(range(data), at), digits=2)
		}
		text(x=at, y=0.6, labels=lab, cex=1, pos=1)
	} # if (plot %in% c("Ct"))
	# The actual plotting of the card
	if (missing(main))
		main <- "Ct values"
	par(mar=c(0,1,2,1), mgp=c(1,0.6,0))
	x <- rep(1:ncol,nrow)/ncol
	y <- rep(nrow:1, each=ncol)/nrow
	plot(x=x, y=y, cex=chamber.size, xaxt="n", yaxt="n", ylab="", xlab="", main=main, bg=values, pch=22, ...)
}

