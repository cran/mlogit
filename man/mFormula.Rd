\name{mFormula}
\alias{mFormula}
\alias{is.mFormula}
\alias{mFormula.formula}
\alias{model.matrix.mFormula}
\alias{model.frame.mFormula}
\title{Model formula for logit models }
\description{
  Two kinds of variables are used in logit models: alternative specific and
  individual specific variables. \code{mFormula} provides a relevant
  class to deal with this specificity and suitable methods to extract
  the elements of the model.
}
\usage{
mFormula(object)
\method{mFormula}{formula}(object)
is.mFormula(object)
\method{model.matrix}{mFormula}(object, data, ...)
\method{model.frame}{mFormula}(formula, data, ..., lhs = NULL, rhs = NULL)
}
\arguments{
  
  \item{object}{for the \code{mFormula} function, a formula, for the
    \code{update} and \code{model.matrix} methods, a \code{mFormula}
    object,}
  \item{formula}{a \code{mFormula} object,}
  \item{data}{a \code{data.frame},}
%  \item{id}{the indvidual index in case of panel data,}
  \item{lhs}{see \code{Formula}}
  \item{rhs}{see \code{Formula}}
  \item{...}{further arguments.}
    
}

\value{

  an object of class \code{mFormula}.
}


\details{
  
  Let \code{J} being the number of alternatives.  The formula may
  include alternative-specific and individual specific variables. For
  the latter, \code{J-1} coefficients are estimated for each
  variable. For the former, only one (generic) coefficient or \code{J}
  different coefficient may be estimated.

  A \code{mFormula} is a formula for which the right hand side may
  contain three parts: the first one contains the alternative specific
  variables with generic coefficient, \emph{i.e.} a unique coefficient
  for all the alternatives ; the second one contains the individual
  specific variables for which one coefficient is estimated for all the
  alternatives except one of them ; the third one contains the
  alternative specific variables with alternative specific coefficients.
  The different parts are separeted by a ``\code{|}'' sign. If a
  standard formula is writen, it is assumed that there are only
  alternative specific variables with generic coefficients.
  
  The intercept is necessarely alternative specific (a generic intercept
  is not identified because only utility differences are
  relevant). Therefore, it deals with the second part of the formula. As
  it is usual in \code{R}, the default behaviour is to include an
  intercept. A model without an intercept (which is hardly meaningfull)
  may be specified by includint \code{+0} or \code{-1} in the second rhs
  part of the formula. \code{+0} or \code{-1} in the first and in the
  third part of the formula are simply ignored. 
 
Specific methods are provided to build correctly the model matrix and to
update the formula. The \code{mFormula} function is not intended to be
use directly. While using the \cite{\code{mlogit}} function, the first
argument is automaticaly coerced to a \code{mFormula} object.

}

\examples{
data("Fishing", package = "mlogit")
Fish <- mlogit.data(Fishing, varying = c(2:9), shape = "wide", choice =
"mode")

# a formula with to alternative specific variables (price and catch) and
# an intercept
f1 <- mFormula(mode ~ price + catch)
head(model.matrix(f1, Fish), 2)

# same, with an individual specific variable (income)
f2 <- mFormula(mode ~ price + catch | income)
head(model.matrix(f2, Fish), 2)

# same, without an intercept
f3 <- mFormula(mode ~ price + catch | income + 0)
head(model.matrix(f3, Fish), 2)

# same as f2, but now, coefficients of catch are alternative specific
f4 <- mFormula(mode ~ price | income | catch)
head(model.matrix(f4, Fish), 2)
}

\author{Yves Croissant}
\keyword{models}