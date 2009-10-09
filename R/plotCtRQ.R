plotCtRQ <-
function(qDE,
	comparison	= 1,
	genes,
	transform	= "log2",
	p.val	= 0.1,
	mark.sig	= TRUE,
	p.sig	= 0.05,
	p.very.sig	= 0.01,
	mark.un	= TRUE,
	un.tar	= "black",
	un.cal	= "black",
	col,
	legend	= TRUE,
	xlim,
	mar,
	main,
	...)
{
	# Get the data for the significance test
	if (class(qDE)=="list") {
		data	<- qDE[[comparison]]
	} else if (class(qDE)=="data.frame") {
		data	<- qDE	
	} else {
		stop("Unexpected format of qDE data.")	
	}
	ddCt	<- data$ddCt
	# Perform the chosen transformation
	if (transform=="none") {
		FC	<- 2^-ddCt	
	} else if (transform=="log10") {
		FC	<- log10(FC	<- 2^-ddCt)
	} else if (transform=="log2") {
		FC	<- log2(FC	<- 2^-ddCt)
	} else {
		stop("Transform \'", transform, "\' not valid.", sep="")
	}
	# Use adjusted p-values if present
	if ("adj.p.value" %in% colnames(data)) {
		p.values	<- data$adj.p.value
	} else {
		p.values	<- data$p.value
	}
	# Select which genes to show
	if (!missing(genes)) {
		if (is.numeric(genes)) {
			index	<- 1:nrow(data) %in% genes
		} else if (is.character(genes)) {
			index	<- data$genes %in% genes
		} else {
			index	<- genes
		}
	} else {
		index	<- p.values < p.val
	}
	if (!sum(index))
		stop("No genes fulfill the selected criteria")
	# Prepare colour schemes
	if (missing(col))
		col	<- "#3288BD"
	# Adjust various plotting parameters
	old.par <- par(no.readonly = TRUE)
	on.exit(par(old.par))
	if (missing(mar)) {
		max	<- max(nchar(data$genes[index]))
		mar	<- c(0.5*max+3,3,2,1)
	}
	if (missing(main))
		if (class(qDE)=="list") {
			main	<- paste("Comparison:", names(qDE)[comparison])
		} else {
			main	<- "Relative quantification"
		}
	if (missing(xlim))
		xlim	<- c(-1, sum(index)*1.2)
	# Make the barplot
	par(mar=mar)
	y.max	<- ifelse(max(FC[index])>0, max(FC[index])*9/8, 0)
	y.min	<- min(FC[index])*9/8
	barplot(FC[index], col=col, las=2, names=data$genes[index], main=main, ylim=c(y.min, y.max), xlim=xlim, ...)
	# Mark significant tests
	if (mark.sig) {
		sig	<- p.values[index] < p.sig & p.values[index] > p.very.sig
		very.sig	<- p.values[index] < p.very.sig
		at	<- seq(0.75, sum(index)*1.2, 1.2)
		if (any(sig))
			mtext("*", side=1, at=at[sig])	
		if (any(very.sig))
			mtext("\"", side=1, at=at[very.sig])
		if (legend)
			legend(xlim[1], y.max*0.95, legend=paste("p-value <", c(p.sig, p.very.sig)), pch=c("*", "\""), bty="n")
	}
	# Mark "reliability" of columns
	if (mark.un) {
		# Define locations of boxes
		if (sum(index)>1) {
			x.left	<- seq(0.2, sum(index)*1.2-1, 1.2)
		} else {
			x.left	<- 0.2
		}
		y.top	<- FC[index]
		# Mark unreliable targets
		tar	<- data$categoryTarget[index] != "OK"
		if (sum(tar)>0)
			rect(x.left[tar], 0, x.left[tar]+1, y.top[tar], angle=45, density=10, col=un.tar)
		# Mark unreliable calibrators
		cal	<- data$categoryCalibrator[index] != "OK"
		if (sum(cal)>0)
			rect(x.left[cal], 0, x.left[cal]+1, y.top[cal], angle=-45, density=10, col=un.cal)
		# Add legend
		if (legend)
			legend(xlim[1], y.min*0.75, legend=c("Target undetermined", "Calibrator undetermined"), bty="n", angle=c(45, -45), density=20)
	}
}

