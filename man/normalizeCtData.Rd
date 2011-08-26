\name{normalizeCtData}
\Rdversion{1.1}
\alias{normalizeCtData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Normalization of Ct values from qPCR data.}
\description{This function is for normalizing Ct data from high-throughput qPCR platforms like the TaqMan Low Density Arrays. Normalization can be either within or across different samples.}

\usage{
normalizeCtData(q, norm = "deltaCt", deltaCt.genes = NULL, scale.rank.samples, rank.type = "pseudo.median", Ct.max = 35, geo.mean.ref, verbose = TRUE)
}

\arguments{
  \item{q}{object of class qPCRset.}
  \item{norm}{character string with partial match allowed, the normalisation method to use. "deltaCt" (default) , "scale.rankinvariant", "norm.rankinvariant", "quantile" and "geometric.mean" are implemented. See details.}
  \item{deltaCt.genes}{character vector, the gene(s) to use for deltaCt normalization. Must correspond to some of the \code{featureNames} in \code{q} or \code{NULL}, in which case the endogenous controls from \code{featureType} are used.}
  \item{scale.rank.samples}{integer, for the "scale.rankinvariant" method, how many samples should a feature be rank invariant across to be included. Defaults to number of samples-1.}
  \item{rank.type}{string, the reference sample for the rank invariant normalisation. Either "pseudo.median" or "pseudo.mean" for using the median or mean across samples as a pseudo-reference sample.}
  \item{Ct.max}{numeric, Ct values above this will be ignored when identifying rank invariant genes.}
  \item{geo.mean.ref}{numeric, the reference sample to scale to for the "geometric.mean" method. Defaults tot he first sample..}
  \item{verbose}{boolean, should some information be printed to the prompt.}
}

\details{"quantile" will make the expression distributions across all cards more or less identical. "deltaCt" calculates the standard deltaCt values, i.e. subtracts the mean of the chosen controls from all other values on the array. "scale.rankinvariant" sorts features from each sample based on Ct values, and identifies a set of features that remain rank invariant, i.e. whose ordering is constant. The average of these rank invariant features is then used to scale the Ct values on each array individually. "norm.rankinvariant" also identifies rank invariant features between each sample and a reference, and then uses these features to generate a normalisation curve individually for each sample by smoothing. "geometric.mean" calculates the geometric mean of all Ct values below Ct.max in each sample, and scales the Ct values accordingly. 
	
For the rank invariant methods it can make a significant difference whether high Ct values, such as "40" or something else being used for undetermined Ct values is removed during the normalisation using the Ct.max parameter. "norm.rankinvariant" also depends on having enough rank invariant genes for generating a robust smoothing curve.

"quantile" is base on \code{normalizeQuantiles} from \code{limma}, and the rank invariant normalisations implement methods from \code{normalize.invariantset} in package \code{affy}.

The distribution of Ct values before/after normalisation can be assessed with the function \code{plotCtDensity}.}

\value{An object of class qPCRset like the input.}

%\references{
%Mar, Kimura, Schroder, Irvine, Hayashizaki, Suzuki, Hume and Quackenbush (2009), Data-driven normalization strategies for high-throughput quantitative RT-PCR,
%BMC Bioinformatics, 10.
%}

\author{Heidi Dvinge}

\seealso{
\code{\link[affy]{normalize.invariantset}} for the rank invariant normalisations, \code{\link[limma]{normalizequantiles}} and \code{\link{plotCtDensity}}
}

\examples{
# Load example data
data(qPCRraw)
# Perform different normalisations
dnorm <- normalizeCtData(qPCRraw, norm="deltaCt", deltaCt.genes="Gene1")
qnorm <- normalizeCtData(qPCRraw, norm="quantile")
nrnorm <- normalizeCtData(qPCRraw, norm="norm.rankinvariant")
srnorm <- normalizeCtData(qPCRraw, norm="scale.rankinvariant")
gnorm <- normalizeCtData(qPCRraw, norm="geometric.mean")
# Normalized versus raw data
cols <- rep(brewer.pal(6, "Spectral"), each=384)
plot(exprs(qPCRraw), exprs(dnorm), pch=20, col=cols, main="deltaCt normalization")
plot(exprs(qPCRraw), exprs(qnorm), pch=20, col=cols,  main="Quantile normalization")
plot(exprs(qPCRraw), exprs(nrnorm), pch=20, col=cols,  main="norm.rankinvariant")
plot(exprs(qPCRraw), exprs(srnorm), pch=20, col=cols,  main="scale.rankinvariant")
plot(exprs(qPCRraw), exprs(gnorm), pch=20, col=cols,  main="geometric.mean")
# With or without removing high Ct values
nrnorm <- normalizeCtData(qPCRraw, norm="norm.rankinvariant")
nrnorm2 <- normalizeCtData(qPCRraw, norm="norm.rankinvariant", Ct.max=40)
plot(exprs(nrnorm), exprs(nrnorm2), pch=20, col=cols, xlab="Ct.max = 35", ylab="Ct.max = 40")
# Distribution of the normalised data
par(mfrow=c(2,3), mar=c(3,3,2,1))
plotCtDensity(qPCRraw, main="Raw Ct values")
plotCtDensity(dnorm, main="deltaCt")
plotCtDensity(qnorm, main="quantile")
plotCtDensity(srnorm, main="scale.rankinvariant")
plotCtDensity(nrnorm, main="norm.rankinvariant")
plotCtDensity(gnorm, main="geometric.mean")
}

\keyword{ htest }
