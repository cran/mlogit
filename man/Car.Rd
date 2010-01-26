\name{Car}
\docType{data}
\alias{Car}
\title{Stated Preferences for Car Choice }
\description{
  a cross-section

 \emph{number of observations} :  4654

 \emph{observation} :  individuals

   \emph{country} :  United States

}
\usage{data(Car)}
\format{A dataframe containing :
  \describe{
    \item{choice}{choice of a vehicule amoung 6 propositions}
    \item{college}{college education ?}
    \item{hsg2}{size of household greater than 2 ?}
    \item{coml5}{commulte lower than 5 miles a day ?}
    \item{typez}{body type, one of regcar (regular car), sportuv (sport utility vehicule), sportcar, stwagon (station wagon), truck, van, for each proposition z from 1 to 6}
    \item{fuelz}{fuel for proposition z, one of gasoline, methanol, cng (compressed natural gas), electric.}
    \item{pricez}{price of vehicule divided by the logarithme of income}
    \item{rangez}{hundreds of miles vehicule  can travel between refuelings/rechargings}
    \item{accz}{acceleration, tens of seconds required to reach 30 mph from stop}
    \item{speedz}{highest attainable speed in hundreds of mph}
    \item{pollutionz}{tailpipe emissions as fraction of those for new gas vehicule}
    \item{sizez}{0 for a mini, 1 for a subcompact, 2 for a compact and 3 for a mid--size or large vehicule}
    \item{spacez}{fraction of luggage space in comparable new gas vehicule}
    \item{costz}{cost per mile of travel (tens of cents) : home recharging for electric vehicule, station refueling otherwise}
    \item{stationz}{fraction of stations that can refuel/recharge vehicule}
  }
}
\source{
    McFadden, Daniel  and  Kenneth  Train (2000) \dQuote{Mixed MNL models for discrete response}, \emph{Journal of Applied Econometrics}, \bold{15(5)}, 447--470.

}
\references{
      Journal of Applied Econometrics data archive : \url{http://jae.wiley.com/jae/}.
}

\keyword{datasets}
