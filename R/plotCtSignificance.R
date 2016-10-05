plotCtSignificance <-
function(qDE,
	q,
	comparison	= 1,
	genes,
	p.val	= 0.1,
	groups,
	calibrator,
	target,
	p.sig	= 0.05,
	p.very.sig	= 0.01,
	mark.sig	= TRUE,
	col,
	un.col	= "#D53E4F",
	point.col	= "grey",
	legend	= TRUE,
	mar,
	main,
	jitter	= 0.5,
	...)
{
	# Get the data for the significance test
	if (class(qDE)=="list") {
		data <- qDE[[comparison]]
	} else if (class(qDE)=="data.frame") {
		data <- qDE	
	} else {
		stop("Unexpected format of qDE data.")	
	}
#	means <- subset(data, select=grep("mean", colnames(data)))
	means <- data[,c("meanTarget", "meanCalibrator")]
	# Use adjusted p-values if present
	if ("adj.p.value" %in% colnames(data)) {
		p.values <- data$adj.p.value
	} else {
		p.values <- data$p.value
	}
	# Select which genes to show
	# Select which genes to show
	if (!missing(genes)) {
		if (is.numeric(genes)) {
			index <- 1:nrow(data) %in% genes
		} else if (is.character(genes)) {
			if (!any(grepl("Gene", data$genes)))
			{
				data$genes <- paste0("Gene", data$genes)
			}
			index <- data$genes %in% genes
		} else {
			index <- genes
		}
	} else {
		index <- p.values < p.val
	}
	if (!sum(index))
		stop("No genes fulfill the selected criteria")
	# Prepare colour schemes
	if (missing(col))
		col <- brewer.pal(7, "Spectral")[c(3,6)]
	col <- rep(col, sum(index))
	# Adjust various plotting parameters
	old.par <- par(no.readonly = TRUE)
	on.exit(par(old.par))
	if (missing(mar)) {
		max <- max(nchar(colnames(exprs(q))[index]), na.rm=TRUE)
		mar <- c(0.4*max+2,3,2,1)
	}
	if (missing(main))
		main <- paste("Comparison:", names(qDE)[comparison])
	# Make the barplot
	par(mar=mar)
	y.max <- max(means[index,])+max(means[index,])/8
	barplot(t(as.matrix(means)[index,]), beside=TRUE, col=col, las=2, names=data$genes[index], main=main, ylim=c(0, y.max), ...)
	if (legend)
		legend(0, y.max, legend=colnames(means), fill=col, bty="n")
	# Mark significant tests
	if (mark.sig) {
		sig <- p.values[index] < p.sig & p.values[index] > p.very.sig
		very.sig <- p.values[index] < p.very.sig
		at <- seq(2, sum(index)*3, 3)
		if (any(sig))
			mtext("*", side=1, at=at[sig])	
		if (any(very.sig))
			mtext("\"", side=1, at=at[very.sig])
	}
	# Add the actual values (might want to add more than two levels at some point)
	if (missing(target) | missing(calibrator)) {
		stop("Target and calibrator samples need to be supplied.")
	}
	if (!is.factor(groups))
		groups <- as.factor(groups)
	l.groups <- length(levels(groups))
	# If the reference present in the groups
	cal <- groups==calibrator
	tar <- groups==target
	# Run through all the genes
	sub.pos <- strsplit(as.character(data$feature.pos[index]), split=";", fixed=TRUE)
	for (i in 1:sum(index)) {
		# The x position in plot
		x.pos <- c(3*i-1.5)
		# The feature(s) to extract data for
		feat.index <- featurePos(q) %in% sub.pos[[i]]
		# Run through all levels.
		for (j in c("tar", "cal")) {
			# Get the data
			sample.index <- get(j)
			point.data <- exprs(q)[feat.index,sample.index]
			point.cats <- featureCategory(q)[feat.index,sample.index]
			# Set the point colour
			bg.col <- rep(point.col, length(point.data))
			bg.col[unlist(point.cats) %in% c("Unreliable", "Undetermined")] <- un.col
			# Set the x positions
			x.pos2 <- ifelse(j=="tar", x.pos, x.pos+1)
			if (jitter==0)
				jitter <- 0.00001
			x.pos2 <- jitter(rep(c(x.pos2), length(point.data)), amount=jitter)
			# The actual plotting
			points(x.pos2, c(point.data), pch=21, bg=bg.col, xpd=TRUE)
		}
	}
}

