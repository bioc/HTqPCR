plotCtOverview <-
function(q,
	cards	= TRUE,
	genes,
	groups,
	calibrator,
	replicates	= TRUE,
	col,
	conf.int	= FALSE,
	legend	= TRUE,
	...) 
{
	# Get the data
	data	<- exprs(q)[,cards]
	if (missing(genes)) {
		index	<- TRUE	
	} else if (is.numeric(genes) | is.logical(genes)) {
		index	<- genes
	} else if (is.character(genes)) {
		index	<- rownames(data) %in% genes
	} else if (is.factor(genes)) {
		genes	<- as.character(genes)
		index	<- rownames(data) %in% genes
	} else {
		stop("Unknown format of 'genes'")	
	}
	data	<- data[index,]
	# If genes are replicated, average across them
	if (replicates) {
		feature.split	<- rownames(data)	
	} else {
		feature.split	<- paste(featureNames(q), featurePos(q), sep="_")[index]
	}
	data.feat	<- split(as.data.frame(data), feature.split)
	# Split into groups if provided
	if (!missing(groups)) {
		sample.split	<- groups	
	} else {
		sample.split	<- colnames(data)	
	}
	data.samp	<- lapply(data.feat, function(x) split(t(x), sample.split))
	# Calculate mean and sd of each
	data.all	<- sapply(data.samp, 
		function(xx) { 
			means <- sapply(xx, mean)
			stdev <- sqrt(sapply(xx, var))
  			n     <- sapply(xx,length)
  			ciw   <- qt(0.975, n) * stdev / sqrt(n)
  			c(M=means, SD=ciw)
		})
	M	<- data.all[grep("M", rownames(data.all)),]
	# If calibrator, take the ratio compared to that
	if (!missing(calibrator)) {
		cali.mean	<- M[grep(calibrator, rownames(M)),]
		M	<- t(log2(t(M)/cali.mean))
		ylab	<- paste("Log2 ratio compared to group", calibrator)
	} else {
		ylab	<- "Ct values for samples"	
	}
	# Set some plotting parameters
	if (missing(col))
		col	<- brewer.pal(10, "Spectral")[c(1,8,9,5,10,2,7,3,4,6)][1:nrow(M)]
	if (legend) {
		bar.names	<- gsub("M\\.(.+)", "\\1", rownames(M))
	} else {
		bar.names	<- NULL
	}	
	# The actual plotting
	barplot(M, beside=TRUE, las=2, col=col, legend.text=bar.names, ylab=ylab, xpd=TRUE, ...)
	# Add error bars if required - ABSOLUTE Ct VALUES
	if	(missing(calibrator) & conf.int) {
		SD	<- data.all[grep("SD", rownames(data.all)),]
		# Get the right X positions
		x.pos	<- seq(1.5, ncol(SD)*(nrow(SD)+1), nrow(SD)+1)
		for (i in 1:nrow(SD)) {
			plotCI(x=x.pos+i-1, y=M[i,], uiw=SD[i,], add=TRUE, gap=0, pch=20, xpd=TRUE, sfra=0.001)
		}
	}
	# Add error bars if required - RELATIVE Ct VALUES 
	if	(!missing(calibrator) & conf.int) {
		# Get the values compared to the average of the mean
		data.samp.ratio	<- lapply(seq_along(data.samp), 
			function(x) 
				sapply(data.samp[[x]], "/", cali.mean[x]))
		# Calculate mean and sd of each 
		SD.ratio	<- sapply(data.samp.ratio, 
			function(xx) { 
				stdev <- sqrt(apply(xx, 2, var))
  				n     <- nrow(xx)
  				ciw   <- qt(0.975, n) * stdev / sqrt(n)
   				c(SD=ciw)
			})
		# Get the right X positions
		x.pos	<- seq(1.5, ncol(SD.ratio)*(nrow(SD.ratio)+1), nrow(SD.ratio)+1)
		# Plot error bars
		for (i in 1:nrow(SD.ratio)) {
			plotCI(x=x.pos+i-1, y=M[i,], uiw=SD.ratio[i,], add=TRUE, gap=0, pch=20, xpd=TRUE, sfra=0.001)
		}
	}
}

