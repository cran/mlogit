\name{MobilePhones}
\docType{data}
\alias{MobilePhones}
\title{Stated Preferences survey for mobile phones }
\description{
  a cross-section from 2003

 \emph{number of observations} :  11184

 \emph{observation} :  individuals

   \emph{country} :  Netherland

}
\usage{data(MobilePhones)}
\format{A dataframe containing :
  \describe{
    \item{alt}{the alternative, denoted by 1 or 2}
    \item{choice}{1 if the alternative is chosen, 0 otherwise}
    \item{price}{purchase price in euros (100, 135 or 170)}
    \item{mincost}{cost of a minute in euros (0.25, 0.30 or 0.35)}
    \item{extras}{extra features of the telephone : a factor with levels
      \code{games}, \code{internet} (which means games and internet) and
      \code{camera} (which means games and internet and camera)}
    \item{network}{a factor with levels \code{KPNVodaphone} and
      \code{other}}
    \item{sms}{the cost of an sms in euros (0.17 or 0.23)}
    \item{design}{a factor with levels \code{basic} or \code{trendy}}
  }
}
\source{
  Sandor Zsolt and Philip Hans Franses (2009) \dQuote{Consumer price
    evaluations through choice experiments}, \emph{Journal of Applied
    Econometrics}, \bold{24}, 517--535.
}
\references{
      Journal of Applied Econometrics data archive : \url{http://jae.wiley.com/jae/}.
}

\keyword{datasets}
