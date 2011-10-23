\name{plot.mlogit}
\alias{plot.mlogit}
\alias{plot.rpar}
\title{Plot of the distribution of estimated random parameters}
\description{
  Methods for \code{rpar} and \code{mlogit} objects which provide a
  plot of the distribution of one or all of the estimated random parameters
}
\usage{
\method{plot}{mlogit}(x, par = NULL, norm = NULL,
                      type = c("density", "probability"), ...)
\method{plot}{rpar}(x, norm = NULL, type = c("density","probability"), ...)
}
\arguments{
  \item{x}{a \code{mlogit} or a \code{rpar} object,}
  \item{type}{the function to be plotted, whether the density or the
    probability density function,}
  \item{par}{a subset of the random parameters ;  if \code{NULL}, all
    the parameters are selected,}
  \item{norm}{the coefficient's name for the \code{mlogit} method or the
    coefficient's value for the \code{rpar} method used for normalization,}
  \item{...}{further arguments, passed to \code{plot.rpar} for the
    \code{mlogit} method and to \code{plot} for the \code{rpar} method.}
}
\details{
  For the \code{rpar} method, one plot is drawn. For the \code{mlogit}
  method, one plot for each selected random parameter is drawn.
}

\seealso{
  \code{\link{mlogit}} for the estimation of random parameters logit models
  and \code{\link{rpar}} for the description of \code{rpar} objects and
  \code{\link{distribution}} for functions which return
  informations about the distribution of random parameters.
}

\author{Yves Croissant}
\keyword{regression}