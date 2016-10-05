normalizeCtData <-
function(q,
	norm	= "deltaCt",
	deltaCt.genes	= NULL,
	scale.rank.samples,
	rank.type	= "pseudo.median",
	Ct.max	= 35,
	geo.mean.ref,
	verbose	= TRUE)
{
	# Extract the data
	data <- exprs(q)
	data.norm <- data
	# Get the normalisation method
	method <- match.arg(norm, c("quantile", "scale.rankinvariant", "norm.rankinvariant", "deltaCt", "geometric.mean"))
	# Some general stuff that will be used by both rank.invariant methods
	if (method %in% c("scale.rankinvariant", "norm.rankinvariant")) {
		# Index to use for too high Ct values
		Ct.index <- data>Ct.max
		data.Ctmax <- data
		data.Ctmax[Ct.index] <- NA
		# Define what to rank against
		if (rank.type=="pseudo.median") {
			ref.data <- apply(data.Ctmax, 1, median, na.rm=TRUE)
		} else if (rank.type=="pseudo.mean") {
			ref.data <- apply(data.Ctmax, 1, mean, na.rm=TRUE)
		}
		# Mark + replace NA values with something temporary
		na.index <- is.na(ref.data)
		ref.data[na.index] <- 30
		# Run the rank.invariant function
		data.rankinvar <- apply(data, 2, normalize.invariantset, ref=ref.data)
	}
	# The actual normalisation
	switch(method,
		quantile = {
			# Use an internal limma function
			data.norm <- normalizeQuantiles(data)
		},
		scale.rankinvariant = {
			# Get all the rank invariant genes
			ri.genes <- sapply(data.rankinvar, "[[", "i.set")
			# Remove those with too high Ct values
			ri.genes[Ct.index] <- FALSE
			# Remove those that were all NA for potentially other reasons
			ri.genes[na.index,] <- FALSE
			# Select those to use here
			ri.count <- rowSums(ri.genes)
			if (missing(scale.rank.samples))
				scale.rank.samples <- ncol(data)-1
			ri.index <- ri.count >= scale.rank.samples
			if (sum(ri.index)==0)
				stop(paste("No rank invariant genes were found across", scale.rank.samples, "samples"))
			# Extract the corresponding Ct values; average
			ri.mean <- colMeans(data[ri.index,,drop=FALSE])
			ri.scale <- ri.mean/ri.mean[1]
			# Correct the data
			data.norm <- t(t(data)*ri.scale)
			# Print info
			if (verbose) {
				cat(c("Scaling Ct values\n\tUsing rank invariant genes:", paste(featureNames(q)[ri.index], collapse=" "), "\n"))
				cat(c("\tScaling factors:", format(ri.scale, digits=3), "\n"))
			}
		},
		norm.rankinvariant = {
			# Print info
			if (verbose)
				cat("Normalizing Ct values\n\tUsing rank invariant genes:\n")
			# Correct the data based on the calculations above
			for (i in 1:ncol(data)) {
				# Check if there are any rank invariant genes
				ri.sub <- data.rankinvar[[i]]
				ri.genes <- ri.sub[["i.set"]]
				# Remove those that don't pass the Ct.max criteria
				ri.genes[Ct.index[,i]] <- FALSE
				# Remove those that are NA for other reasons
				ri.genes[na.index] <- FALSE
				if (sum(ri.genes)==0) {
	  		warning(paste("\tNo rank invariant genes were found for sample ", sampleNames(q)[i], "; sample not normalized\n", sep=""))
	  		next
	  	}
	  	# If verbose, print some info
	  	if (verbose)
	  		cat(paste("\t", sampleNames(q)[i], ": ", sum(ri.genes), " rank invariant genes\n", sep=""))
	  # The actual correction
	  data.norm[,i] <- as.numeric(approx(ri.sub$n.curve$y, ri.sub$n.curve$x, xout=data[,i], rule=2)$y)
	 }
 		},
		deltaCt	= {
			# Which are the reference genes (endogenous controls)
			if (is.null(deltaCt.genes))
				deltaCt.genes <- unique(featureNames(q)[featureType(q)=="Endogenous Control"])
			c.index <- featureNames(q) %in% deltaCt.genes
			if (verbose) {
				cat(c("Calculating deltaCt values\n\tUsing control gene(s):", paste(deltaCt.genes, collapse=" "), "\n"))
			}
			# Run though all cards; perform internal normalisation
			for (c in 1:ncol(data)) {
				# Calculate the control genes
				refCt <- mean(data[c.index,c], na.rm=TRUE)
				refsd <- sd(data[c.index,c], na.rm=TRUE)
				# Difference for target and controls
				data.norm[,c] <- data[,c]-refCt
				# Print results
				if (verbose)
					cat(paste("\tCard ", c, ":\tMean=", format(refCt, dig=4), "\tStdev=", format(refsd, dig=3), "\n", sep=""))
			}
		},
		geometric.mean = {
			# For each column, calculate the geometric mean of Ct values<Ct.max
			geo.mean <- apply(data, 2, function(x) {
									xx <- log2(subset(x, x<Ct.max))
									2^mean(xx)})
			# Which sample to scale to
			if (missing(geo.mean.ref))
				geo.mean.ref <- 1
			# Calculate the scaling factor
			geo.scale <- geo.mean/geo.mean[geo.mean.ref]
			# Adjust the data accordingly
			data.norm <- t(t(data) * geo.scale)
			if (verbose) {
				cat(c("Scaling Ct values\n\tUsing geometric mean within each sample\n"))
				cat(c("\tScaling factors:", format(geo.scale, digits=3), "\n"))
			}
		} # switch
	)
	# Replace with the normalised Ct exprs
	exprs(q) <- data.norm
	# Add to the history of the object
	if (nrow(getCtHistory(q))==0)
		setCtHistory(q) <- data.frame(history="Manually created qPCRset object.", stringsAsFactors=FALSE)
	setCtHistory(q) <- rbind(getCtHistory(q), capture.output(match.call(normalizeCtData)))
	# Return the normalised object
	q
}

