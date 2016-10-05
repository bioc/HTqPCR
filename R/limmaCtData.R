limmaCtData <-
function(q,
	design	= NULL,
	contrasts,
	sort	= TRUE,
	stringent	= TRUE,
	ndups	= 1,
	spacing	= NULL,
	dupcor,
	...)
{
	# Get the data
	data <- exprs(q)
	featPos <- featurePos(q)
	# If there are duplicates, calculate the correlation between them, unless there are technical replicates
	if (missing(dupcor)) {
		if (ndups>1) {
			dup.cor <- duplicateCorrelation(data, ndups=ndups, spacing=spacing, design=design)
			temp <- unwrapdups(featPos, ndups=ndups, spacing=spacing)
			featPos <- apply(temp, 1, paste, collapse=";")
		} else {
			dup.cor <- NULL
		}
	}
	# Do the limma fitting
	fit <- lmFit(data, design=design, ndups=ndups, spacing=spacing, correlation=dup.cor$consensus, ...)
	if (!missing(contrasts))
		fit <- contrasts.fit(fit, contrasts=contrasts)
	fit2 <- eBayes(fit)
	# Prepare output
	out <- list()
	if (!missing(contrasts)) {
		coefs <- colnames(contrasts)
		cont <- design %*% contrasts
	} else {
		coefs <- colnames(design)
		cont <- design
	}
	# Summarize result of each test
	for (coef in coefs) {
		# The actual results
		res <- topTable(fit2, coef=coef, number=nrow(fit2), sort.by="none", ...)
		both.means <- both.cats <- array(0, c(nrow(res), 2), list(rownames(res), c("Test", "Reference")))
		# Additional info
		for (i in c(-1,1)) {
			sample <- ifelse(i==1, "Test", "Reference")
			# The test and the reference samples in the given comparison
			index <- cont[,coef]==i
			# The mean
			mean <- rowMeans(unwrapdups(data[,index], ndups=ndups, spacing=spacing))
			# Summarize data quality categories
			new.cat <- rep("OK", length(mean))
			old.cat <- unwrapdups(featureCategory(q)[,index], ndups=ndups, spacing=spacing)
			count.cat <- apply(old.cat, 1, function(x) sum(x %in% c("Undetermined", "Unreliable")))
			cutoff <- ifelse(stringent, 1, ceiling(sum(index)/2))
			new.cat[count.cat>=cutoff] <- "Undetermined"
			# Add to objects
			both.means[, sample] <- mean
			both.cats[, sample] <- new.cat
		}
		# Combine info
#		res.out <- cbind(res[,c("ID", "t", "P.Value", "adj.P.Val", "logFC")], meanTest, meanReference, categoryTest, categoryReference)
		res.out <- cbind(rownames(res), featPos, res[,c("t", "P.Value", "adj.P.Val", "logFC")], 2^(-res$logFC), both.means, both.cats)
		colnames(res.out) <- c("genes", "feature.pos", "t.test", "p.value", "adj.p.value", "ddCt", "FC", "meanTarget", "meanCalibrator", "categoryTarget", "categoryCalibrator")
		# Assign to output
		if (sort)
			res.out <- res.out[order(res.out$adj.p.value),]
		out[[coef]] <- res.out
	}
	# Brief summary across all contrasts/tests
	res <- decideTests(fit2, ...)
	rownames(res) <- rownames(topTable(fit2, sort="none", n=nrow(fit2)))
	out[["Summary"]] <- res
	# Return output
	out
}

