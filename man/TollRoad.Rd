\name{TollRoad}
\docType{data}
\alias{TollRoad}
\title{Stated Preferences survey for a toll road }
\description{
  a panel from 1999-2000

 \emph{number of observations} :  1448 observations for 548 individuals

 \emph{observation} :  individuals

   \emph{country} :  United-States

}
\usage{data(TollRoad)}
\format{A dataframe containing :
  \describe{
    \item{id}{the individual id}
    
    \item{src}{the source of the data, one of \code{br} : revealed
      preference survey conducted by the Brookings Institution,
      \code{bs} : stated preferences survey conducted by the Brookings
      Institution and \code{cal} : revealed preferences survey conducted
      by the Californian Polytechnic State University }
    \item{route}{the route chosen, one of \code{express} (the toll-road)
      or \code{freeway} (the free road)}
    \item{toll.alt}{the monetary cost of the road (a for the free road)}
    \item{time.alt}{the median time of the trip on both highways for the
      given schedule}
    \item{reliability.alt}{the reliability of the trip length on both
      highways for the given schedule, measured by the difference
      between the 80th and the 50th percentile of the trip length}
    \item{occupance}{the number of people in the car}
    \item{size}{the household size}
    \item{sex}{one of \code{male} or \code{female}}
    \item{flexibility}{does the respondent declare having a flexible
      arrival time, a factor with levels \code{yes} or \code{no}}
    \item{distance}{trip distance in miles}
    \item{commute}{a long-commute for trips longer than 45 minutes (a
      factor with levels \code{yes} or \code{no})}
    \item{age3050}{a factor with levels \code{yes} if the respondent is
      between 30 and 50 years old, \code{no} otherwise}
    \item{income}{a factor with levels \code{low}, \code{medium} and \code{high}}
  }
}
\source{
  Kenneth A. Small, Clifford Winston, Jia Yan (2005) \dQuote{Uncovering
      the distribution of motorists' preferences for travel time and
      reliability}, \emph{Econometrica}, \bold{73(4)}, 1367-1382.
}
\references{
      Econometrica data archive
}

\keyword{datasets}
