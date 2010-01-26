\name{hmftest}
\alias{hmftest}
\alias{hmftest.formula}
\alias{hmftest.mlogit}
\title{Hausman-McFadden Test}
\description{
  Test the IIA hypothesis (independence of irrelevant alternatives) for
  a multinomial logit model.
}
\usage{
hmftest(x, ...)
\method{hmftest}{mlogit}(x, z, ...)
\method{hmftest}{formula}(x, alt.subset, ...)
}
\arguments{
  
  \item{x}{an object of class \code{mlogit} or a formula,}

  \item{z}{an object of class \code{mlogit} or a subset of alternatives
    for the \code{mlogit} method. This should be the same model as
    \code{x} estimated on a subset of alternatives,}

  \item{alt.subset}{a subset of alternatives,}

  \item{...}{further arguments passed to \code{mlogit} for the
    \code{formula} method.}
  
}
\value{
  an object of class \code{"htest"}.
}
\details{

  This is an implementation of the Hausman's consistency test for
  multinomial logit models. If the independance of irrelevant
  alternatives applies, the probability ratio of every two alternatives
  depends only on the characteristics of these
  alternatives. Consequentely, the results obtained on the estimation
  with all the alternatives or only on a subset of them are consistent,
  but more efficient in the first case. On the contrary, only the
  results obtained from the estimation on a relevant subset are
  consistent. To compute this test, one needs a model estimated with all
  the alternatives and one model estimated on a subset of
  alternatives. This can be done by providing two objects of class
  \code{mlogit}, one object of class \code{mlogit} and a character
  vector indicating the subset of alternatives, or a formula and a
  subset of alternatives.

  
}
\references{

  Hausman, J.A. and D. McFadden (1984), A Specification Test for the
  Multinomial Logit Model,
  \emph{Econometrica}, \bold{52}, pp.1219--1240.

}

\examples{

## from Greene's Econometric Analysis p. 731

data("TravelMode",package="AER")
TravelMode <- mlogit.data(TravelMode,choice="choice",shape="long",
                          alt.var="mode",id.var="individual")

## Create a variable of income only for the air mode

TravelMode$avinc <- with(TravelMode,(mode=='air')*income)

## Estimate the model on all alternatives, with car as the base level
## like in Greene's book.

#x <- mlogit(choice~wait+gcost+avinc,TravelMode,reflevel="car")
x <- mlogit(choice~wait+gcost+avinc,TravelMode)

## Estimate the same model for ground modes only (the variable avinc
## must be dropped because it is 0 for every observation

g <- mlogit(choice~wait+gcost,TravelMode,reflevel="car",
            alt.subset=c("car","bus","train"))

## Compute the test

hmftest(x,g)

}
\author{Yves Croissant}
\keyword{htest}