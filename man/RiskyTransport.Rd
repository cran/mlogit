\name{RiskyTransport}
\docType{data}
\alias{RiskyTransport}
\title{Risky Transportation Choices }
\description{
  a cross-section

 \emph{number of observations} :  5405

   \emph{country} :  Sierra Leone

}
\usage{data(RiskyTransport)}
\format{A dataframe containing :
  \describe{
    \item{id}{individual id}
    \item{choice}{1 for the chosen mode}
    \item{mode}{one of \code{Helicopter}, \code{WaterTaxi}, \code{Ferry} and \code{Hovercraft}}
    \item{cost}{the generalised cost of the transport mode}
    \item{risk}{the fatality rate, numbers of death per 100,000 trips}
    \item{weight}{weights}
    \item{seats}{}
    \item{noise}{}
    \item{crowdness}{}
    \item{convloc}{}
    \item{clientele}{}
    \item{chid}{choice situation id}
    \item{african}{\code{yes} if born in Africa, \code{no} otherwise}
    \item{lifeExp}{declared life expectancy}
    \item{dwage}{declared hourly wage}
    \item{iwage}{imputed hourly wage}
    \item{educ}{level of education, one of \code{low} and \code{high}}
    \item{fatalism}{self-ranking of the degree of fatalism}
    \item{gender}{gender, one of \code{female} and \code{male}}
    \item{age}{age}
    \item{haveChildren}{\code{yes} if the traveler has children, \code{no} otherwise}
    \item{swim}{\code{yes} if the traveler knows how to swim, \code{no} itherwise}
  }
}
\source{
    American Economic Association data archive : \url{http://aeaweb.org/aer/}.

}
\references{
    Gianmarco, Leon  and  Edward  Miguel (2017) \dQuote{Risky Transportation Choices and the Value of a Statistical Life}, \emph{American Economic Journal: Applied Economics}, \bold{9(1)}, 202-28.

}

\keyword{datasets}
