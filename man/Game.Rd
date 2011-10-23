\name{Game}
\docType{data}
\alias{Game}
\alias{Game2}
\title{Ranked data for gaming platforms}
\description{
  a cross-section

  \emph{number of observations} :  91
  
  \emph{observation} :  individuals
  
  \emph{country} :  Netherlands

}
\usage{data(Game)}
\format{A dataframe containing :
  \describe{
    \item{ch.Platform}{where \code{platform} is one of \code{Xbox},
      \code{PlayStation}, \code{PSPortable}, \code{GameCube},
      \code{GameBoy} and \code{PC}. This variables contain the ranking
      of the platforms from 1 to 6,}
    \item{own.Platform}{these 6 variables are dummies which indicate
      whether the given plaform is already owned by the respondent,}
    \item{age}{the age of the respondent,}
    \item{hours}{hours per week spent on gaming.}
  }
}

\details{

  The data are also provided in long format (use in this case
  \code{data(Game2)}). In this case, the alternative and the choice
  situation are respectively indicated in the \code{platform} and
  \code{chid} variables.

}

\source{
  
  Denis Fok, Richard Paap, and Bram van Dijk (2010) \dQuote{A
    Rank-Ordered Logit Model with Unobserved Heterogeneity in Ranking
    Capabilities}, Journal of Applied Econometrics, forthcoming

}
\references{
      Journal of Applied Econometrics data archive : \url{http://jae.wiley.com/jae/}.
}

\keyword{datasets}
