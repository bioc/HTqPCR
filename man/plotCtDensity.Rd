\name{plotCtDensity}
\Rdversion{1.1}
\alias{plotCtDensity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Distribution plot for qPCR Ct values.}
\description{Function for plotting the density distribution of Ct values from high-throughput qPCR data.}

\usage{
plotCtDensity(q, cards = TRUE, xlab = "Ct", ylab = "Density", col, main = NULL, legend = TRUE, lwd = 2, ...)
}

\arguments{
 \item{q}{object of class qPCRset. }
 \item{cards}{vector, the numbers of the cards to plot. Defaults to TRUE = all cards.}
 \item{xlab}{character string, label for the x-axis.}
 \item{ylab}{character string, label for the y-axis.}
 \item{col}{vector of colours to use, defaults to different colour for each card.}
 \item{main}{character string, plot title.}
 \item{legend}{logical, whether to include a colour legend or not.}
 \item{lwd}{numeric, the width of the lines.}
 \item{\dots}{any other arguments will be passed to the \code{matplot} function.}
}

\details{The distribution of Ct values in the qPCRset \code{q} is calculated using \code{density}.
}
	
\value{A plot is created on the current graphics device.}

\author{Heidi Dvinge}

\seealso{\code{\link{matplot}}, \code{\link{density}}
}

\examples{
# Loading the data
data(qPCRraw)
# Make plot with all samples or just a few
plotCtDensity(qPCRraw)
plotCtDensity(qPCRraw, cards=c(1,4))
}

\keyword{hplot}
