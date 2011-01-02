\name{mlogit}
\alias{mlogit}
\alias{print.mlogit}
\alias{suml}
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
mlogit(formula, data, subset, weights, na.action, start = NULL,
       alt.subset = NULL, reflevel = NULL, 
       nests = NULL, un.nest.el = FALSE, unscaled = FALSE,
       heterosc = FALSE, rpar = NULL,
       R = 40, correlation = FALSE, halton = NULL,
       random.nb = NULL, panel = FALSE, estimate = TRUE, ...)
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
  \item{data}{the data: an \code{mlogit.data} object or an ordinary
  \code{data.frame},}
  \item{subset}{an optional vector specifying a subset of observations,}
  \item{weights}{an optional vector of weights,}
  \item{na.action}{a function which indicates what should happen when
    the data contains '\code{NA}'s,}
  \item{start}{a vector of starting values,}
  \item{alt.subset}{a vector of character strings containing the subset of
  alternative on which the model should be estimated,}
  \item{reflevel}{the base alternative (the one for which the
    coefficients of individual-specific variables are normalized to 0),}
  \item{estimate}{a boolean indicating whether the model should be
    estimated or not: if not, the \code{model.frame} is returned,}
  \item{nests}{a named list of characters vectors, each names being a
  nest, the corresponding vector being the set of alternatives that
  belong to this nest,}
  \item{un.nest.el}{a boolean, if \code{TRUE}, the hypothesis of unique
    elasticity is imposed for nested logit models,}
  \item{unscaled}{a boolean, if \code{TRUE}, the unscaled version of the
  nested logit model is estimated,}
  \item{heterosc}{a boolean, if \code{TRUE}, the heteroscedastic logit
    model is estimated,}
  \item{rpar}{a named vector whose names are the random parameters and
  values the distribution : \code{'n'} for normal, \code{'l'} for
  log-normal, \code{'t'} for truncated normal, \code{'u' } for uniform,}
  \item{R}{the number of function evaluation for the gaussian quadrature
  method used if \code{heterosc=TRUE}, the number of draws of
  pseudo-random numbers if \code{rpar} is not \code{NULL},}
  \item{correlation}{only relevant if \code{rpar} is not \code{NULL}, if
  true, the correlation between random parameters is taken into
  account,}
  \item{halton}{only relevant if \code{rpar} is not \code{NULL}, if not
  \code{NULL}, halton sequence is used instead of pseudo-random
  numbers. If \code{halton=NA}, some default values are used for the
  prime of the sequence (actually, the primes are used in order) and for
  the number of elements droped. Otherwise, \code{halton} should be a
  list with elements \code{prime} (the primes used) and \code{drop} (the
  number of elements droped).}
  \item{random.nb}{only relevant if \code{rpar} is not \code{NULL}, a
    user-supplied matrix of random,}
  
  \item{panel}{only relevant if \code{rpar} is not \code{NULL} and if
  the data are repeated observations of the same unit ; if \code{TRUE},
  the mixed-logit model is estimated using panel techniques,}

  \item{digits}{the number of digits,}
  \item{width}{the width of the printing,}
  \item{outcome}{a boolean which indicates, for the \code{fitted} and the
  \code{residuals} methods whether a matrix (for each choice, one value
  for each alternative) or a vector (for each choice, only a value for
  the alternative chosen) should be returned,}
  \item{...}{further arguments passed to \code{mlogit.data} or
  \code{mlogit.optim}.}
  
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
  \item{formula}{the formula (a \code{mFormula} object),}
  \item{expanded.formula}{the formula (a \code{formula} object),}
  \item{model}{the model frame used,}
  \item{index}{the index of the choice and of the alternatives.}
}
\details{
  
  For how to use the formula argument, see \code{\link{mFormula}}.

  The \code{data} argument may be an ordinary \code{data.frame}. In this
  case, some supplementary arguments should be provided and are passed
  to \code{\link{mlogit.data}}. Note that it is not necessary to indicate the
  choice argument as it is deduced from the formula.

  The model is estimated using the \code{\link{mlogit.optim}} function. 

  The basic multinomial logit model and three important extentions of
  this model may be estimated.

  If \code{heterosc=TRUE}, the heteroscedastic logit model is
  estimated. \code{J-1} extra coefficients are estimated that represent
  the scale parameter for \code{J-1} alternatives, the scale parameter
  for the reference alternative being normalized to 1. The probabilities
  don't have a closed form, they are estimated using a gaussian
  quadrature method.

  If \code{nests} is not \code{NULL}, the nested logit model is
  estimated.
  
  If \code{rpar} is not \code{NULL}, the random parameter model is
  estimated. The probabilities are approximated using simulations with
  \code{R} draws and halton sequences are used if \code{halton} is not
  \code{NULL}. Pseudo-random numbers are drawns from a standard normal
  and the relevant transformations are performed to obtain numbers
  drawns from a normal, log-normal, censored-normal or uniform
  distribution. If \code{correlation=TRUE}, the correlation between the
  random parameters are taken into account by estimating the components
  of the cholesky decomposition of the covariance matrix. With G random
  parameters, without correlation G standard deviations are estimated,
  with correlation G * (G + 1) /2 coefficients are estimated.
  
}

\examples{

## Cameron and Trivedi's Microeconometrics p.493 There are two
## alternative specific variables : price and catch one individual
## specific variable (income) and four fishing mode : beach, pier, boat,
## charter

data("Fishing", package = "mlogit")
Fish <- mlogit.data(Fishing, varying = c(2:9), shape = "wide", choice = "mode")

## a pure "conditional" model

summary(mlogit(mode ~ price + catch, data = Fish))

## a pure "multinomial model"

summary(mlogit(mode ~ 0 | income, data = Fish))

## which can also be estimated using multinom (package nnet)

library("nnet")
summary(multinom(mode ~ income, data = Fishing))

## a "mixed" model

m <- mlogit(mode ~ price+ catch | income, data = Fish)
summary(m)

## same model with charter as the reference level

m <- mlogit(mode ~ price+ catch | income, data = Fish, reflevel = "charter")

## same model with a subset of alternatives : charter, pier, beach

m <- mlogit(mode ~ price+ catch | income, data = Fish,
            alt.subset = c("charter", "pier", "beach"))

## model on unbalanced data i.e. for some observations, some
## alternatives are missing

# a data.frame in wide format with two missing prices
Fishing2 <- Fishing
Fishing2[1, "price.pier"] <- Fishing2[3, "price.beach"] <- NA
mlogit(mode~price+catch|income, Fishing2, shape="wide", choice="mode", varying = 2:9)

# a data.frame in long format with three missing lines
data("TravelMode", package = "AER")
Tr2 <- TravelMode[-c(2, 7, 9),]
mlogit(choice~wait+gcost|income+size, Tr2, shape = "long",
       chid.var = "individual", alt.var="mode", choice = "choice")

## An heteroscedastic logit model

data("TravelMode", package = "AER")
hl <- mlogit(choice ~ wait + travel + vcost, TravelMode,
             shape = "long", chid.var = "individual", alt.var = "mode",
             method = "bfgs", heterosc = TRUE, tol = 10)

## A nested logit model

TravelMode$avincome <- with(TravelMode, income * (mode == "air"))
TravelMode$time <- with(TravelMode, travel + wait)/60
TravelMode$timeair <- with(TravelMode, time * I(mode == "air"))
TravelMode$income <- with(TravelMode, income / 10)

#Heiss p.231
nl <- mlogit(choice ~ time + timeair | income, TravelMode, shape = "long",
             alt.var = "mode", print.level = 0, method = "bfgs",
             nests = list(public = c("train", "bus"), other = c("air","car")),
             tol = 1)

# Hensher and Greene (2002), table 1 p.8-9 model 5
TravelMode$incomeother <- with(TravelMode, ifelse(mode \%in\% c('air', 'car'), income, 0))
nl <- mlogit(choice~gcost+wait+incomeother, TravelMode,
             shape='long', alt.var='mode',
             nests=list(plublic=c('train', 'bus'), other=c('car','air')))

# same with a comon nest elasticity (model 1)
nl2 <- update(nl, un.nest.el = TRUE)

## a mixed logit model
rpl <- mlogit(mode ~ price+ catch | income, Fishing, varying = 2:9,
              shape = 'wide', rpar = c(price= 'n', catch = 'n'),
              correlation = TRUE, halton = NA,
              R = 10, tol = 10, print.level = 0)
summary(rpl)
rpar(rpl)
cor.mlogit(rpl)
cov.mlogit(rpl)
rpar(rpl, "catch")
summary(rpar(rpl, "catch"))

}

\seealso{
  
  \code{\link{mlogit.data}} to shape the data. \code{\link{multinom}} from
  package \code{nnet} performs the estimation of the multinomial logit
  model with individual specific variables. \code{\link{mlogit.optim}}
  for details about the optimization function.
  
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
