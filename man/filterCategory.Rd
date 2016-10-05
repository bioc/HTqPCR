\name{filterCategory}
\Rdversion{1.1}
\alias{filterCategory}
\title{Filter Ct values based on their feature categories.}
\description{
Ct values corresponding to selected feature categories will be replaced by NA. Generally, the feature categories indicate how reliable the values are.
}

\usage{
filterCategory(q, na.categories = c("Unreliable", "Undetermined"))
}

\arguments{
 \item{q}{a qPCRset object.}
 \item{na.categories}{character vector, with the name(s) of the feature categories where Ct values will be considered NA.}
}

%\details{
%% ~~ If necessary, more details than the description above ~~
%}

\value{A qPCRset object like the input, but with the selected Ct values replaced by NAs}

\author{Heidi Dvinge}

\seealso{
\code{\link{setCategory}} for adjusting the categories.
}

\examples{
data(qPCRraw)
qPCRraw2 <- setCategory(qPCRraw, groups=NULL)
x <- filterCategory(qPCRraw2)
summary(qPCRraw)
summary(x)
}

\keyword{ htest }
