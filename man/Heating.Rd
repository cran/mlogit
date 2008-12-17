\name{Heating}
\docType{data}
\alias{Heating}
\title{Heating System Choice in California Houses }
\description{
  a cross-section

 \emph{number of observations} :  900

 \emph{observation} :  households

   \emph{country} :  California

}
\usage{data(Heating)}
\format{A dataframe containing :
  \describe{
    \item{idcase}{id}
    \item{depvar}{heating system, one of gc (gas central), gr (gas room), ec (electric central), er (electric room), hp (heat pump)}
    \item{ic.z}{installation cost for heating system z (defined for the 5 heating systems)}
    \item{oc.z}{annual operating cost for heating system z (defined for the 5 heating systems)}
    \item{pb.z}{ratio oc.z/ic.z }
    \item{income}{annual income of the household}
    \item{agehed}{age of the household head}
    \item{rooms}{numbers of rooms in the house}
  }
}
\references{
      Kenneth Train's home page :  \url{http://elsa.berkeley.edu/~train/}.
}

\keyword{datasets}
