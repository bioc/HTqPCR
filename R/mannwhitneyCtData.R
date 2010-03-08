mannwhitneyCtData <-
function(q,
	groups	= NULL,
	calibrator, 
	alternative	= "two.sided",
	paired	= FALSE,
	replicates	= TRUE,
	sort	= TRUE,
	stringent	= TRUE,
	p.adjust	= "BH",
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
	# Perform the Mann-Whitney test. 
	mw.tests	<- lapply(data2, function(x) {
		x	<- as.matrix(x)
		res	<- wilcox.test(x[,g1], x[,g2], alternative=alternative, paired=paired, exact=FALSE, ...)
		# Calculate mean for each group 
		res[["estimate"]]	<- c(mean(x[,g1]), mean(x[,g2]))
		res
		})
	# Collect output
	means	<- t(sapply(mw.tests, "[[", "estimate"))
	colnames(means)	<- c("meanCalibrator", "meanTarget")
	p.value	<- sapply(mw.tests, "[[", "p.value")
	mw.value	<- sapply(mw.tests, function(x) paste(names(x[["statistic"]]), "=", x[["statistic"]]))
	genes	<- sapply(feats, "[[", 1)
	featurePos	<- sapply(featPos, paste, collapse=";")
	# Make adjusted p-value
	adj.p.value	<- p.adjust(p.value, method=p.adjust)
	# Fold change is calculated as ddCt, as well as 2^(-ddCT)
	cal	<- grep(calibrator, colnames(means))
	FC	<- means[, "meanTarget"] - means[, "meanCalibrator"]
	FC2	<- 2^(-FC)
	out	<- data.frame(genes, featurePos, mw.value, p.value, adj.p.value, FC, FC2, means, row.names=1:length(genes))
	# Indicate reliability of measure
	for (l in unique(groups)) {
		new.cat	<- rep("OK", length(data2))
		old.cat	<- split(featureCategory(q[,groups==l]), split.by)
		count.cat <- sapply(old.cat, function(x) sum(unlist(x) %in% c("Undetermined", "Unreliable")))
		cutoff	<- ifelse(stringent, 1, ceiling(sum(groups==l)/2))
		new.cat[count.cat>=cutoff]	<- "Undetermined"
		out[, paste("category", ifelse(l==calibrator, "Calibrator", "Target"), sep="")]	<- new.cat
	}
	# Return output, sorted by p-value if requested
	names(out)	<- c("genes", "feature.pos", "MW.value", "p.value", "adj.p.value", "ddCt", "FC", colnames(means), grep("category", colnames(out), value=TRUE))
	if (sort)
		out	<- out[order(out$p.value),]
	out
}

