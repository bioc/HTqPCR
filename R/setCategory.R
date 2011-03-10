setCategory <-
function(q,
	Ct.max	= 35,
	Ct.min	= 10,
	replicates	= TRUE,
	quantile	= 0.90,
	groups,
	flag	= TRUE,
	flag.out	= "Failed",
	verbose	= TRUE,
	plot	= FALSE,
	...)
{
	# The output
	out	<- q
	# Get the Ct values
	data	<- exprs(q)

	## LOOKING JUST AT Ct VALUES
	# Mark high/low Ct values
	featureCategory(out)[data > Ct.max]	<- "Undetermined"
	featureCategory(out)[data < Ct.min]	<- "Unreliable"
	# Print summary of how many per sample
	if (verbose) {
		# Count. Don't use apply, since not all samples contain "Unreliable"
		feats	<- featureCategory(out)
		cats	<- sort(unique(unlist(feats)))
		count	<- array(0, c(length(cats), n.samples(out)), list(cats, sampleNames(out)))
		for (i in 1:ncol(feats)) {
			tab	<- table(feats[,i])
			count[names(tab),i]	<- tab	
		}
		# Print the info
		cat("Categories after Ct.max and Ct.min filtering:\n")
		print(count)			
	}

	## INCLUDE FLAGGING BY THE SOFTWARE
	if (flag) {
		flags	<- flag(q)
		featureCategory(out)[flags %in% flag.out] <- "Unreliable"
	}

	## INCLUDE STANDARD DEVIATION ACROSS REPLICATES
	# Collapse replicated genes if required
	if (replicates) {
		split.by	<- rownames(data)
	} else {
		split.by	<- paste(rownames(data), "_well", 1:nrow(data), sep="")
	}
	data2 <- split(as.data.frame(data), split.by)
	# Also looks at standard deviation
	if (!is.null(quantile)) {
		# Divide into groups
		if(missing(groups)) {
			warning("Sample groups must be supplied to filter by standard deviation.")
			invisible(out)
		}
		l.groups	<- unique(groups)
		# Calculate mean and standard deviations across technical and biological replicates
		SD	<- AV <- N <- array(0, dim=c(length(data2), length(l.groups)), dimnames=list(names(data2), l.groups))
		for (g in l.groups) {
			SD[,g]	<- sapply(data2, function(d) sd(unlist(d[,groups==g])))
			AV[,g]	<- sapply(data2, function(d) mean(unlist(d[,groups==g])))
			N[,g]	<- sapply(data2, function(d) length(unlist(d[,groups==g])))
		}
		# Plot overall standard deviation for each sample type
		if (plot) {
			for (g in seq_along(l.groups)) {
				dev.new()
				par(mfrow=c(1,2), mar=c(4,4,2,1))
				hist(SD[,g], n=25, main=l.groups[g], xlab="Standard deviation")
				plot(AV[,g], SD[,g], pch=20, xlab="Mean", ylab="Standard deviation", main=l.groups[g], ...)
			}
		}
		# Calculate conf. intervals. NB: Assumes normal distribution!!
		for (d in seq_along(data2)) {
			gene	<- names(data2)[d]
			for (g in l.groups) {
				# Get data for gene and group
				sd	<- SD[gene, g]
				n	<- N[gene,g]
				av	<- AV[gene,g]
				# Calculate error and upper/lower bounds
				bounds	<- 0.5+quantile/2
#				ciw   <- qt(bounds, n) * sd / sqrt(n)
##				error <- qnorm(bounds)*sd/sqrt(n)
##				conf	<- av+c(-1,1)*error
				conf	<- vector("numeric", 2)
				conf[1] <- qnorm(bounds, mean=av, sd=sd)
				conf[2] <- qnorm(1-bounds, mean=av, sd=sd)
				# Mark Ct values beyond the limts
				Ct	<- data[split.by==gene,groups==g]
				index <- Ct > conf[1] | Ct < conf[2]
#				index <- Ct > av+ciw | Ct < av-ciw
				featureCategory(out)[split.by==gene,groups==g][index] <- "Unreliable"
				# "remark" the undetermined one - this has priority over unreliable
				featureCategory(out)[data > Ct.max] <- "Undetermined"
			}
		}
		# Print summary of how many per sample
		if (verbose & !is.null(groups)) {
			# Count. Don't use apply, since not all samples contain "Unreliable"
			feats	<- featureCategory(out)
			cats	<- sort(unique(unlist(feats)))
			count	<- array(0, c(length(cats), n.samples(out)), list(cats, sampleNames(out)))
			for (i in 1:ncol(feats)) {
				tab	<- table(feats[,i])
				count[names(tab),i]	<- tab	
			}
			# Print the info
			cat("Categories after standard deviation filtering:\n")
			print(count)			
		}
	}
	# Add to the history of the object
	if (nrow(getCtHistory(out))==0)
		out@history	<- data.frame(history="Manually created qPCRset object.", stringsAsFactors=FALSE)
	out@history	<- rbind(out@history, capture.output(match.call(setCategory)))
	# Return the filtered object, if desired
	invisible(out)
}

