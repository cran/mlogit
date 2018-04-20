\name{Train}
\docType{data}
\alias{Train}
\title{Stated Preferences for Train Traveling }
\description{
  a cross-section from 1987

 \emph{number of observations} :  2929

 \emph{observation} :  individuals

   \emph{country} :  Netherland

}
\usage{data(Train)}
\format{A dataframe containing :
  \describe{
    \item{id}{individual identifient}
    \item{choiceid}{choice identifient}
    \item{choice}{one of 'A' or 'B'}
    \item{price\_z}{price of proposition z (z = 'A', 'B') in cents of guilders}
    \item{time\_z}{travel time  of proposition z (z = 'A', 'B') in minutes}
    \item{comfort\_z}{comfort  of proposition z (z = 'A', 'B'), 0, 1 or 2 in decreasing comfort order}
    \item{change\_z}{number of changes for proposition z (z = 'A', 'B')}
  }
}
\source{
    Ben--Akiva, M., D.  Bolduc and M.  Bradley (1993) \dQuote{Estimation
      of travel choice models with randomly distributed values of time},
    \emph{Transportation Research Record}, \bold{1413}, 88--97.

    Meijer, Erik  and  Jan  Rouwendal (2006) \dQuote{Measuring welfare
      effects in models with random coefficients}, \emph{Journal of
      Applied Econometrics}, \bold{21}(2), 227-244.

}
\references{
      Journal of Applied Econometrics data archive : \url{http://jae.wiley.com/jae/}.
}

\keyword{datasets}
