\name{effects.mlogit}
\alias{effects.mlogit}
\title{Marginal effects of the covariates}
\description{
  The \code{effects} method for \code{mlogit} objects computes the
  marginal effects of the selected covariate on the probabilities of
  choosing the alternatives
}
\usage{
\method{effects}{mlogit}(object, covariate = NULL,
                        type = c("aa", "ar", "rr", "ra"), data = NULL, ...)
}
\arguments{
  \item{object}{a \code{mlogit} object,}
  \item{covariate}{the name of the covariate for which the effect should
    be computed,}
  \item{type}{the effect is a ratio of two marginal variations of the
    probability and of the covariate ; these variations can be absolute
    \code{"a"} or relative \code{"r"}. This argument is a string that
    contains two letters, the first refers to the probability, the second
    to the covariate,}
  \item{data}{a data.frame containing the values for which the effects
    should be calculated. The number of lines of this data.frame should be
    equal to the number of alternatives,}
  \item{...}{further arguments.}
}

\value{
  If the covariate is alternative specific, a $J$ times $J$ matrix is
  returned, $J$ being the number of alternatives. Each line contains the
  marginal effects of the covariate of one alternative on the
  probability to choose any alternative. If the covariate is individual
  specific, a vector of length $J$ is returned.
}

\seealso{
  \code{\link{mlogit}} for the estimation of multinomial logit models.
}

\examples{
data("Fishing", package = "mlogit")
Fish <- mlogit.data(Fishing, varying = c(2:9), shape = "wide", choice = "mode")
m <- mlogit(mode ~ price | income | catch, data = Fish)
# compute a data.frame containing the mean value of the covariates in
# the sample
z <- with(Fish, data.frame(price = tapply(price, index(m)$alt, mean),
                           catch = tapply(catch, index(m)$alt, mean),
                           income = mean(income)))
# compute the marginal effects (the second one is an elasticity
effects(m, covariate = "income", data = z)
effects(m, covariate = "price", type = "rr", data = z)
effects(m, covariate = "catch", type = "ar", data = z)
}

\author{Yves Croissant}
\keyword{regression}