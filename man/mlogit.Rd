\name{mlogit}
\alias{mlogit}
\alias{print.mlogit}

\alias{fitted.mlogit}
\alias{residuals.mlogit}
\alias{df.residual.mlogit}
\alias{terms.mlogit}
\alias{model.matrix.mlogit}
\alias{update.mlogit}

\alias{summary.mlogit}
\alias{print.summary.mlogit}
\alias{vcov.mlogit}
\alias{logLik.mlogit}
\title{Multinomial logit model}
\description{
  Estimation by maximum likelihood of the multinomial logit model, with
  alternative-specific and/or individual specific variables.
}
\usage{
mlogit(formula, data, subset, weights, na.action,
       alt.subset = NULL, reflevel = NULL, ...)
\method{print}{mlogit}(x, digits = max(3, getOption("digits") - 2),
    width = getOption("width"), ...)
\method{summary}{mlogit}(object, ...)
\method{print}{summary.mlogit}(x, digits = max(3, getOption("digits") - 2),
    width = getOption("width"), ...)
\method{print}{mlogit}(x, digits = max(3, getOption("digits") - 2),
    width = getOption("width"), ...)
\method{logLik}{mlogit}(object, ...)
\method{vcov}{mlogit}(object, ...)
\method{residuals}{mlogit}(object, outcome = TRUE, ...)
\method{fitted}{mlogit}(object, outcome = TRUE, ...)
\method{df.residual}{mlogit}(object, ...)
\method{terms}{mlogit}(x, ...)
\method{model.matrix}{mlogit}(object, ...)
\method{update}{mlogit}(object, new, ...)
}

\arguments{
  \item{x, object}{an object of class \code{mlogit}}
  \item{formula}{a symbolic description of the model to be estimated,}
  \item{new}{an updated formula for the \code{update} method,}
  \item{data}{the data,}
  \item{subset}{an optional vector specifying a subset of observations,}
  \item{weights}{an optional vector of weights,}
  \item{na.action}{a function which indicates what should happen when
  the data contains '\code{NA}'s,}
  \item{alt.subset}{a vector of character strings containing the subset of
  alternative on which the model should be estimated,}
  \item{reflevel}{the base alternative (the one for which the
    coefficients of individual-specific variables are normalized to 0),}
  \item{digits}{the number of digits,}
  \item{width}{the width of the printing,}
  \item{outcome}{a boolean which indicates, for the \code{fitted} and the
  \code{residuals} methods whether a matrix (for each choice, one value
  for each alternative) or a vector (for each choice, only a value for
  the alternative chosen) should be returned,}
  \item{...}{further arguments.}
}

\value{
  An object of class \code{"mlogit"}, a list with elements:

  \item{coefficients}{the named vector of coefficients,}
  \item{logLik}{the value of the log-likelihood,}
  \item{hessian}{the hessian of the log-likelihood at convergence,}
  \item{gradient}{the gradient of the log-likelihood at convergence,}
  \item{call}{the matched call,}
  
  \item{est.stat}{some information about the estimation (time used,
    optimisation method),}
  
  \item{freq}{the frequency of choice,}
  \item{residuals}{the residuals,}
  \item{fitted.values}{the fitted values,}
  \item{formula}{the formula (a \code{logitform} object),}
  \item{expanded.formula}{the formula (a \code{formula} object),}
  \item{model}{the model frame used,}
  \item{index}{the index of the choice and of the alternatives.}

}
\details{
  
  Let \code{J} being the number of alternatives.  The formula may
  include alternative-specific and individual specific variables. For the
  latter, \code{J-1} coefficients are estimated for each
  variable. Alternative and individual specific variables are separated by a
  \code{|}. For example, if \code{x1} and \code{x2} are alternative specific
  and \code{z1} and \code{z2} are individual specific, the formula
  \code{y~x1+x2|z1+z2} describe a model with one coefficient for
  \code{x1} and \code{x2} and \code{J-1} coefficients for \code{z1} and
  \code{z2}. \code{J-1} intercepts are also estimated. A model without
  intercepts is defined by the formula : \code{y~x1+x2-1|z1+z2}. To
  obtain alternative specific coefficients for the alternative-specific
  variable \code{x2}, use : \code{y~x1+x2+x2:alt|z1+z2} (replace
  \code{alt} by the relevant variable name if the alternative index is
  provided). Models with only alternative-specific or individual-specific
  variables are respectively estimated by \code{y~x1+x2} and
  \code{y~1|z1+z2}.

  The model is estimated with the \code{maxLik} package and the
  Newton-Raphson method, using analytic gradient and hessian.
  
}

\examples{

## Cameron and Trivedi's Microeconometrics p.493
## There are two alternative specific variables : pr (price) and ca (catch)
## and four fishing mode : beach, pier, boat, charter

data("Fishing",package="mlogit")
Fish <- mlogit.data(Fishing,varying=c(4:11),shape="wide",choice="mode")

## a pure "conditional" model without intercepts

summary(mlogit(mode~pr+ca-1,data=Fish))

## a pure "multinomial model"

summary(mlogit(mode~1|income,data=Fish))

## which can also be estimated using multinom (package nnet)

library(nnet)
summary(multinom(mode~income,data=Fishing))

## a "mixed" model

m <- mlogit(mode~pr+ca|income,data=Fish)
summary(m)

}

\seealso{
  \cite{\code{mlogit.data}} to shape the data. \code{multinom}
  from package \code{nnet} performs the estimation of the multinomial
  logit model with individual specific variables
}

\author{Yves Croissant}

\references{
  
  McFadden, D. (1973) Conditional Logit Analysis of Qualitative Choice
  Behavior, in P. Zarembka ed., \emph{Frontiers in Econometrics},
  New-York: Academic Press.
  
  McFadden, D. (1974) ``The Measurement of Urban Travel Demand'',
  \emph{Journal of Public Economics}, 3, pp. 303-328.
  
  Train, K. (2004) \emph{Discrete Choice Modelling, with Simulations},
  Cambridge University Press.
}
\keyword{regression}
