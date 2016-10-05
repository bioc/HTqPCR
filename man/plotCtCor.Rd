\name{plotCtCor}
\Rdversion{1.1}
\alias{plotCtCor}
\title{Correlation between Ct values from qPCR data}
\description{Function for plotting the correlation based on Ct values between samples containing high-throughput qPCR data.}

\usage{
plotCtCor(q, col, col.range = c(0, 1), main, mar, ...)
}

\arguments{
 \item{q}{object of class qPCRset. }
%% \item{plot}{character string among "Ct" (default), "dCt" and "ddCt". See Details for a longer description.}
 \item{col}{vector of colours to use, defaults to a spectrum from red to blue/purple.}
 \item{col.range}{vector, the range of colours to use.}
 \item{main}{character string, plot title.}
 \item{mar}{vector, the size of the borrom and right hand side margins.}
 \item{\dots}{any other arguments will be passed to the \code{heatmap.2} function.}
}

\details{This function may be used to cluster the samples based on Ct values and present the result in a heatmap. Per default the colours are a rainbow scale from 0 to 1.

The correlation is calculated as 1 - the 'Pearson' method. Prior to version 1.9.1 the value plotted was the correlation directly, rather than 1-correlation.
	
A standard heatmap is drawn, but this can be modified extensively using the arguments available in the \code{heatmap.2} function.}

\value{A plot is created on the current graphics device.}

\author{Heidi Dvinge}

\seealso{\code{\link[gplots]{heatmap.2}}
}

\examples{
data(qPCRraw)
plotCtCor(qPCRraw)
plotCtCor(qPCRraw, col.range=c(0,0.6))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
