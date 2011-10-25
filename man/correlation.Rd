\name{correlation}
\alias{cor.mlogit}
\alias{cov.mlogit}
\title{Correlation structure of the random parameters}
\description{
  Functions that extract the correlation structure of a mlogit object
}
\usage{
cor.mlogit(x)
cov.mlogit(x)
}

\arguments{
  \item{x}{an \code{mlogit} object with random parameters and
    \code{correlation=TRUE}.}
}

\value{
  A numerical matrix which returns either the correlation or the
  covariance matrix of the random parameters.
}

\author{Yves Croissant}

\keyword{regression}
