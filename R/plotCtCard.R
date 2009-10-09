plotCtCard <-
function(q,
	card	= 1,
	plot	= "Ct",
	main,
	nrow	= 16,
	ncol	= 24,
	col,
	col.range,
	na.col	= "grey",
	na.value	= 40,
	legend.cols,
	well.size	= 3.1,
	zero.center	= FALSE,
	unR	= FALSE,
	unD	= FALSE,
	...)
{
	# Define the plot layout
	layout(rbind(2,1), widths=5, heights=c(3.5,1), respect=TRUE)
	# Extract the gene categories
	if (class(q)=="qPCRset") {
		category	<- featureCategory(q)[,card]
	} 
	
	# Ct or other values on a continuous scale (might add other later)
	if (plot %in% c("Ct")) {
		# Extract the values of interest
		data	<- exprs(q)[,card]
		data[data==na.value]	<- NA
		# Define the colours for the card
		if (missing(col)) {
			col	<- rev(colorRampPalette(brewer.pal(9, "YlGnBu"))(20))
		} 		
		if (missing(col.range)) {
			if (zero.center) {
				limit	<- max(abs(data), na.rm=TRUE)
				breaks	<- cut(data, breaks=seq(-limit, limit, length.out=length(col)))
			} else {
				breaks	<- cut(data, breaks=length(col))
			}
		} else {
			breaks	<- cut(data, seq(col.range[1], col.range[2], length.out=length(col)+1), labels=FALSE)
			breaks[data<col.range[1]]	<- 1
			breaks[data>col.range[2]]	<- length(col)
		}
		values	<- col[breaks]
		values[is.na(values)]	<- na.col
		# Plot the legend
		par(mar=c(0,4,0,1))
		l	<- length(col)
		data	<- data[!is.na(data)]
		at	<- seq(0,1,length.out=5)
		plot(0:1,0:1, type="n", xaxt="n", yaxt="n", ylab="", bty="n", xlab="")
		rect(seq(0,1,length.out=l+1)[-(l+1)],0.6,seq(0,1,length.out=l+1)[-1], 0.75, col=col)
		if (!missing(col.range)) {
			lab	<- format(quantile(col.range, at), digits=2)
			lab	<- paste(c("<", rep("", length(at)-2), ">"), lab, sep="")
		} else if (zero.center) {
			lab	<- format(quantile(c(-limit, limit), at), digits=2)
		} else {
			lab	<- format(quantile(range(data), at), digits=2)
		}
		text(x=at, y=0.6, labels=lab, cex=1, pos=1)
	}
	# Quality check, i.e. passed/failed or similar
	else {
		# Extract the values of interest
		if (plot=="flag") {
			data	<- as.factor(flag(q)[,card]) 
		} else if (plot=="type") {
			data	<- featureType(q) 
		} else if (plot=="class") {
			data	<- featureClass(q)
		} else if (plot=="category") {
			data	<- as.factor(featureCategory(q)[,card])
		} else {
			stop(paste("Plot type \'", plot, "\' isn't implemented\n", sep=""))
		}
		l	<- levels(data)
		# Define the colours for the card
		if (missing(col)) {
			if (length(l) < 3) {
				col	<- brewer.pal(10, "Spectral")[c(5,7)]
			} else {
				col	<- rev(brewer.pal(length(l), "Spectral"))
			}
		}
		values	<- rep("white", length(data))
		for (i in 1:length(l)) {
			values[data==l[i]]	<- col[i]
		}
		# Plot the legend
		if (missing(legend.cols)) 
			legend.cols	<- length(l)
		par(mar=c(0,2.5,0,1))
		plot(0:1,0:1, type="n", xaxt="n", yaxt="n", ylab="", bty="n", xlab="")
		legend(x=0, y=1, legend=l, pch=19, col=col, bty="n", pt.cex=well.size, ncol=legend.cols)
	}
	# User supplied factor?
	## TO DO 
	# The actual plotting of the card
	if (missing(main))
		main	<- sampleNames(q)[card]
	par(mar=c(2,4,2,1), mgp=c(1,0.6,0))
	x	<- rep(1:ncol,nrow)/ncol
	y	<- rep(nrow:1, each=ncol)/nrow
	plot(x=x, y=y, cex=well.size, xaxt="n", yaxt="n", ylab="", xlab="", main=main, bg=values, pch=21, ...)
	if (unR) {
		index	<- category=="Unreliable"
		points(x=x[index], y=y[index], cex=well.size, pch=4, ...)
	}
	if (unD) {
		index	<- category=="Undetermined"
		points(x=x[index], y=y[index], cex=well.size, pch=3, ...)
	}
	axis(1, at=1:ncol/ncol, labels=1:ncol, cex.axis=0.8)
	axis(2, at=1:nrow/nrow, labels=LETTERS[nrow:1], las=2, cex.axis=0.8)
	axis(2, at=(1:(nrow*2)/(nrow*2))[seq(3,nrow*2,4)], labels=paste("Port", (nrow/2):1), las=2, cex.axis=1, line=1, tick=FALSE)
#	segments(-0.04, (1:(nrow*2)/(nrow*2))[seq(2,nrow*2,4)], -0.04, (1:(nrow*2)/(nrow*2))[seq(4,nrow*2,4)], xpd=TRUE, lwd=2)
}

