\name{Beef}
\docType{data}
\alias{Beef}
\title{Choice of a beef product }
\description{
  panel data

 \emph{number of observations} :  5 scenario for 47 individuals

 \emph{observation} :  individuals

   \emph{country} :  United States

}
\usage{data(Beef)}
\format{A dataframe containing :
  \describe{
    \item{id}{individuals identifier}
    \item{product}{one of \code{fresh}, \code{lean}, \code{dietLean},
      \code{organic} and \code{noPurchase}}
    \item{scenario}{scenario identifier}
    \item{choice}{one if the product is chosen, 0 otherwise}
    \item{price}{the price of the product}
  }
}
\source{
  
  Jae Bong Chang and Jayson L. Lusk (2010) \dQuote{Mixed logit models:
    accuracy and software choice}, \emph{Journal of Applied
    Econometrics}.

  Chang, J.B., J.L. Lusk, and F.B. Norwood (2009) \dQuote{How Closely Do
Hypothetical Surveys and Laboratory Experiments Predict Field Behavior?}
\emph{American Journal of Agricultural Economics}, \bold{91},
pp. 518-34.

}
\references{
        Journal of Applied Econometrics data archive : \url{http://jae.wiley.com/jae/}.
}

\keyword{datasets}
