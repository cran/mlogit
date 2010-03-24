\name{Electricity}
\docType{data}
\alias{Electricity}
\title{Stated preference data for the choice of electricity suppliers}
\description{
  panel data

 \emph{number of observations} :  4308

 \emph{observation} :  households

   \emph{country} :  United States

}
\usage{data(Electricity)}
\format{A dataframe containing :
  \describe{
    \item{choice}{the choice of the individual, one of 1, 2, 3, 4,}
    \item{id}{the individual index,}
    \item{pfi}{fixed price at a stated cents per kWh, with the price
      varying over suppliers and experiments, for scenario i=(1, 2, 3,
      4),}
    \item{cli}{the length of contract that the supplier offered, in
years (such as 1 year or 5 years.) During this contract period, the
supplier guaranteed the prices and the buyer would have to pay a penalty
if he/she switched to another supplier. The supplier could offer no
contract in which case either side could stop the agreement at any
time. This is recorded as a contract length of 0}
\item{loci}{is the supplier a local company,}
\item{wki}{is the supplier a well-known company}
    \item{todi}{a time-of-day rate under which the price is 11 cents per
  kWh from 8am to 8pm and 5 cents per kWh from 8pm to 8am. These TOD
  prices did not vary over suppliers or experiments: whenever the
  supplier was said to offer TOD, the prices were stated as above.}
    \item{seasi}{a seasonal rate under which the price is 10 cents per kWh in the
  summer, 8 cents per kWh in the winter, and 6 cents per kWh in the
  spring and fall. Like TOD rates, these prices did not vary.  Note
  that the price is for the electricity only, not transmission and
  distribution, which is supplied by the local regulated utility.}
  }
}
\references{
  Kenneth Train's home page :  \url{http://elsa.berkeley.edu/~train/}.
}

\keyword{datasets}
