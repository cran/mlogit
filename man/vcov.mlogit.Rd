\name{vcov.mlogit}
\alias{vcov.mlogit}
\title{vcov method for mlogit objects}
\description{
  The \code{vcov} method for \code{mlogit} objects extract the
  covariance matrix of the coefficients, the errors or the random
  parameters.
}
\usage{
\method{vcov}{mlogit}(object, what = c('coefficient', 'errors', 'rpar'),
                        type = c('cov', 'cor', 'sd'), reflevel = NULL, ...)
}
\arguments{
  \item{object}{a \code{mlogit} object,}
  \item{what}{indicates which covariance matrix has to be extracted :
    the default value is \code{coefficients}, in this case, \code{vcov}
    behaves as usual. If \code{what} equals \code{errors} the covariance
    matrix of the errors of the model is returned. Finally, if
    \code{what} equals \code{rpar}, the covariance matrix of the random
    parameters are extracted.}
  \item{type}{with this argument, the covariance matrix may be returned
    (the default) ; the correlation matrix and the standard deviation
    vector may also be extracted.}
  \item{reflevel}{relevent for the extraction of the errors of a
    multinomial probit model ; in this case the covariance matrix is of
    error differences is returned and, with this argument, the
    alternative used for differentiation is indicated.}
  \item{...}{further arguments.}
}
\details{
  This new interface replaces the \code{cor.mlogit} and
  \code{cov.mlogit} functions which are deprecated.
}

\seealso{
  \code{\link{mlogit}} for the estimation of multinomial logit models.
}

\author{Yves Croissant}
\keyword{regression}