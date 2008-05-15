\name{mlogit}
\alias{mlogit}
\alias{print.mlogit}
\alias{summary.mlogit}
\alias{print.summary.mlogit}
\alias{vcov.mlogit}
\alias{logLik.mlogit}
\title{Multinomial logit model}
\description{
  Estimation by maximum likelihood of the multinomial logit model, with
  choice-specific and/or individual specific variables.
}
\usage{
mlogit(formula, data, weights = NULL, ...)
\method{print}{mlogit}(x, digits = max(3, getOption("digits") - 2),
    width = getOption("width"), ...)
\method{summary}{mlogit}(object, ...)
\method{print}{summary.mlogit}(x, digits = max(3, getOption("digits") - 2),
    width = getOption("width"), ...)
\method{print}{mlogit}(x, digits = max(3, getOption("digits") - 2),
    width = getOption("width"), ...)
\method{logLik}{mlogit}(object, ...)
\method{vcov}{mlogit}(object, ...)
}

\arguments{
  \item{x, object}{an object of class \code{mlogit}}
  \item{formula}{a symbolic description of the model to be estimated,}
  \item{data}{the data,}
  \item{weights}{an optional vector of weights,}
  \item{digits}{the number of digits,}
  \item{width}{the width of the printing,}
  \item{...}{further arguments.}
}

\value{
An object of class \code{"mlogit"}, a list with elements
\code{coefficients}, \code{logLik}, \code{hessian},
\code{gradient}, \code{call}, \code{est.stat},
\code{residuals} and \code{fitted.values}.
}
\details{
  
  Let \code{J} being the number of alternatives.  The formula may
  include choice-specific and individual specific variables. For the
  latter, \code{J-1} coefficients are estimated for each
  variable. Choice and individual specific variables are separated by a
  \code{|}. For example, if \code{x1} and \code{x2} are choice specific
  and \code{z1} and \code{z2} are individual specific, the formula
  \code{y~x1+x2|z1+z2} describe a model with one coefficient for
  \code{x1} and \code{x2} and \code{J-1} coefficients for \code{z1} and
  \code{z2}. \code{J-1} intercepts are also estimated. A model without
  intercepts is defined by the formula : \code{y~x1+x2-1|z1+z2}. To
  obtain alternative specific coefficients for the choice-specific
  variable \code{x2}, use : \code{y~x1+x2+x2:alt|z1+z2} (replace
  \code{alt} by the relevant variable name if the alternative index is
  provided. Models with only choice-specific or individual-specific
  variables are respectively estimated by \code{y~x1+x2} and
  \code{y~1|z1+z2}.

  The model is estimated with the \code{maxLik} package and the
  Newton-Raphson method, using analytic gradient and hessian.
  
}

\examples{
# Heating data, from the Ecdat package
data("Heating",package="Ecdat")

# Heating is a "horizontal" data.frame with three choice-specific
# variables (ic: investment cost, oc: operating cost) and some
# individual-specific variables (income, region, rooms)
Heatingh <- mlogit.data(Heating,cvar=c(ic=3,oc=8,pb=17),
                       shape="hor.var",choice="depvar")

# a model with two choice-specific variables
summary(mlogit(depvar~ic+oc,data=Heatingh))

# same wihtout intercept
summary(mlogit(depvar~ic+oc-1,data=Heatingh))

# a model with choice-specific and individual-specific variables
summary(mlogit(depvar~ic+oc|income+rooms,data=Heatingh))

# a model with choice-specific coefficients for a choice-specific variable
summary(mlogit(depvar~ic+oc+oc:alt,data=Heatingh))

# a model with only individual-specific variables
summary(mlogit(depvar~1|income+rooms,data=Heatingh))

# the same model estimated with multinom from the nnet package
library(nnet)
summary(multinom(depvar~income+rooms,data=Heating))

}

\seealso{
  \cite{\code{mlogit.data}} to shape the data. \code{multinom}
  from package \code{nnet} performs the estimation of the multinomial
  logit model with individual specific variables
}

\author{Yves Croissant}

\references{
  
  McFadden, D. (1973) Conditional logit analysis of qualitative choice
  behavior, in P. Zarembka ed., \emph{Frontiers in Econometrics},
  New-York: Academic Press.
  
  McFadden, D. (1974) ``The Measurement of Urban Travel Demand'',
  \emph{Journal of Public Economics}, 3, pp. 303-328.
  
  Train, K. (2004) \emph{Discrete Choice Modelling, whith Simulations},
  Cambridge University Press.
}
\keyword{regression}
