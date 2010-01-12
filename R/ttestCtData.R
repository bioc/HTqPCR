ttestCtData <-
function(q,
	groups	= NULL,
	calibrator, 
	alternative	= "two.sided",
	paired	= FALSE,
	replicates	= TRUE,
	sort	= TRUE,
	stringent	= TRUE,
	...)
{
	# Get the relevant data
	data	<- exprs(q)
	# Collapse replicated genes if required
	if (replicates) {
		split.by	<- rownames(data)
	} else {
		split.by	<- featurePos(q)
	}
	data2 <- split(as.data.frame(data), split.by)
	feats	<- split(featureNames(q), split.by)
	featPos	<- split(featurePos(q), split.by)
	# Various checks
	if (length(groups) != ncol(data))
		stop("Dimensions of data and groups doesn't match\n")
	groups	<- factor(groups, levels=unique(groups))
	if (length(levels(groups)) != 2)
		stop("Two factor levels required for \'groups\'\n")
	# Assign calibrator and target samples
	if (missing(calibrator)) 
		calibrator	<- groups[1]
	g1	<- groups==calibrator
	g2	<- groups!=calibrator
	# Perform the t-test. 
	t.tests	<- lapply(data2, function(x) {
		x	<- as.matrix(x)
		if (all(x==x[1,1])) {
			# Have to remove samples where all values are identical!
			list(p.value=1, statistic=NA, estimate=c(x[1,1], x[1,1]))
		} else {
			res	<- t.test(x[,g1], x[,g2], alternative=alternative, paired=paired, ...)
			# Calculate mean for each group (not available if paired=TRUE)
			res[["estimate"]]	<- c(mean(x[,g1]), mean(x[,g2]))
			res
		}})
	# Collect output
	means	<- t(sapply(t.tests, "[[", "estimate"))
	colnames(means)	<- c("meanCalibrator", "meanTarget")
	p.value	<- sapply(t.tests, "[[", "p.value")
	t.value	<- sapply(t.tests, "[[", "statistic")
	genes	<- sapply(feats, "[[", 1)
	featurePos	<- sapply(featPos, paste, collapse=";")
	# Fold change is calculated as ddCt
	cal	<- grep(calibrator, colnames(means))
	FC	<- means[, "meanTarget"] - means[, "meanCalibrator"]
	out	<- data.frame(genes, featurePos, t.value, p.value, FC, means, row.names=1:length(genes))
	# Indicate reliability of measure
	for (l in unique(groups)) {
		new.cat	<- rep("OK", length(data2))
		old.cat	<- split(featureCategory(q[,groups==l]), split.by)
		count.cat	<- sapply(old.cat, function(x) sum(x %in% c("Undetermined", "Unreliable")))
		cutoff	<- ifelse(stringent, 1, ceiling(sum(groups==l)/2))
		new.cat[count.cat>=cutoff]	<- "Undetermined"
		out[, paste("category", ifelse(l==calibrator, "Calibrator", "Target"), sep="")]	<- new.cat
	}
	# Return output, sorted by p-value if requested
	names(out)	<- c("genes", "feature.pos", "t.test", "p.value", "ddCt", colnames(means), grep("category", colnames(out), value=TRUE))
	if (sort)
		out	<- out[order(out$p.value),]
	out
}

