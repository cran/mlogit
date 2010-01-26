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
    \item{choice}{one of choice1, choice2}
    \item{pricez}{price of proposition z (z=1,2) in cents of guilders}
    \item{timez}{travel time  of proposition z (z=1,2) in minutes}
    \item{comfortz}{comfort  of proposition z (z=1,2), 0, 1 or 2 in decreasing comfort order}
    \item{changez}{number of changes for proposition z (z=1,2)}
  }
}
\source{
    Meijer, Erik  and  Jan  Rouwendal (2005) \dQuote{Measuring welfare effects in models with random coefficients}, \emph{Journal of Applied Econometrics}, \bold{forthcoming}.

    Ben--Akiva, M.,  D.   Bolduc and  M.  Bradley (1993) \dQuote{Estimation of travel choice models with randomly distributed values of time}, \emph{Transportation Research Record}, \bold{1413}, 88--97.

    Carson, R.T.,  L.   Wilks and  D.  Imber (1994) \dQuote{Valuing the preservation of Australia's Kakadu conservation zone}, \emph{Oxford Economic Papers}, \bold{46}, 727--749.

}
\references{
      Journal of Applied Econometrics data archive : \url{http://jae.wiley.com/jae/}.
}

\keyword{datasets}
