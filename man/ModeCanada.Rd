\name{ModeCanada}
\docType{data}
\alias{ModeCanada}
\title{Mode Choice for the Montreal-Toronto Corridor}
\description{
  a cross-section

 \emph{number of observations} :  3880

 \emph{observation} :  individuals

}
\usage{data(ModeCanada)}
\format{A dataframe containing :
  \describe{
    \item{case}{the individual index}
    \item{alt}{the alternative, one of train, car, bus and air,}
    \item{choice}{one if the mode is chosen, zero otherwise,}
    \item{cost}{monetary cost,}
    \item{ivt}{in vehicule time,}
    \item{ovt}{out vehicule time,}
    \item{frequency}{frequency,}
    \item{income}{income,}
    \item{urban}{urban,}
    \item{noalt}{the number of alternatives available.}
  }
}



\references{

  Bhat, Chandra R. (1995) \dQuote{A heteroscedastic extreme value model
  of intercity travel mode choice}, \emph{Transportation Research Part
  B}, \bold{29(6)}, 471-483.
  
  Koppelman Franck S.  and Chieh-Hua Wen (2001) \dQuote{The paired
    combinatorial logit model:properties, estimation and application},
    \emph{Transportation Research Part B}, 75-89.
  
  Wen, Chieh-Hua and Franck S. Koppelman (2001) \dQuote{The generalized
    nested logit model}, \emph{Transportation Research Part B}, 627-641.
  
}

\examples{
data("ModeCanada", package = "mlogit")
bususers <- with(ModeCanada, case[choice == 1 & alt == "bus"])
ModeCanada <- subset(ModeCanada, !case \%in\% bususers)
ModeCanada <- subset(ModeCanada, nchoice == 4)
ModeCanada <- subset(ModeCanada, alt != "bus")
ModeCanada$alt <- ModeCanada$alt[drop = TRUE]
KoppWen00 <- mlogit.data(ModeCanada, shape='long', chid.var = 'case',
                         alt.var = 'alt', choice='choice',
                         drop.index=TRUE)
pcl <- mlogit(choice~freq+cost+ivt+ovt, KoppWen00, reflevel='car',
              nests='pcl', constPar=c('iv_train_air'))

% use the subset of users who have the four alternatives available and who
% don't take the bus

% bususers <- with(ModeCanada, case[choice == 1 & alt == "bus"])
% ModeCanada <- subset(ModeCanada, !case %in% bususers)
% ModeCanada <- subset(ModeCanada, nchoice == 4)
% ModeCanada <- subset(ModeCanada, alt != "bus")
% ModeCanada$alt <- ModeCanada$alt[drop = TRUE]
% KoppWen00 <- mlogit.data(ModeCanada, shape='long', chid.var = 'case',
%                          alt.var = 'alt', choice='choice',
%                          drop.index=TRUE)
% mlogit(choice~freq+cost+ivt+ovt, KoppWen00, reflevel='car')
% z <- mlogit(choice~freq+cost+ivt+ovt, KoppWen00, reflevel='car', nests="pcl")

}


\keyword{datasets}
