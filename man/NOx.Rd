\name{NOx}
\docType{data}
\alias{NOx}
\title{Technologies to reduce NOx emissions}
\description{
  a cross-section

 \emph{number of observations} :  9480

   \emph{country} :  United States

\emph{economic topic} : microeconomics

\emph{econometrics topic} : discrete choice

}
\usage{data("NOx")}
\format{A dataframe containing :
  \describe{
    \item{chid}{the plant id,}
    \item{alt}{the alternative,}
    \item{id}{the owner id,}
    \item{choice}{the chosen alternative,}
    \item{available}{a dummy indicating that the alternative is available,}
    \item{env}{the regulatory environment, one of \code{'regulated'}, \code{'deregulated'} and \code{'public'},}
    \item{post}{dummy for post-combustion polution control technology,}
    \item{cm}{dummy for combustion modification technology,}
    \item{lnb}{dummy for low NOx burners technology,}
    \item{age}{age of the plant (in deviation from the mean age).}
    \item{vcost}{variable cost,}
    \item{kcost}{capital cost,}
  }
}
\source{
    American Economic Association data archive : \url{http://aeaweb.org/aer/}.

}
\references{
    Fowlie, Meredith (2010) \dQuote{Emissions Trading, Electricity Restructuring, and Investment in Pollution Abatement}, \emph{American Economic Review}, \bold{100(3)}, 837-869.

}

\keyword{datasets}
