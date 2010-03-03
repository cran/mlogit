\name{mlogit.data}
\alias{mlogit.data}
\title{data.frame for logit model}
\description{
  shape a \code{data.frame} in a suitable form for the use of the \code{mlogit} function.
}
\usage{
mlogit.data(data, choice, shape = c("wide","long"), varying = NULL,
            sep=".",alt.var = NULL, chid.var = NULL, alt.levels = NULL,
            id = NULL, opposite = NULL, drop.index = FALSE, ...)
}

\arguments{

  \item{data}{a \code{data.frame},}

  \item{choice}{the variable indicating the choice made: it can be
  either a logical vector, a numerical vector with 0 where the
  alternative is not chosen, a factor with level 'yes' when the
  alternative is chosen}
  
  \item{shape}{the shape of the \code{data.frame}: whether \code{long} if each
    row is an alternative or \code{wide} if each row is an observation,}

  \item{varying}{the indexes of the variables that are alternative
    specific,}
  
  \item{sep}{the seperator of the variable name and the alternative name
  (only relevant for a \code{wide} \code{data.frame}),}
  
  \item{alt.var}{the name of the variable that contains the alternative
    index (for a \code{long} \code{data.frame} only) or the name under which the
    alternative index will be stored (the default name is \code{alt}),}
  
  \item{chid.var}{the name of the variable that contains the choice index
    or the name under which the choice index will be stored,}

  \item{alt.levels}{the name of the alternatives: if null, for
  a \code{wide} data.frame, they are guessed from the variable names and
  the choice variable (both should be the same), for a \code{long}
  \code{data.frame}, they are guessed from the \code{alt.var} argument,}

  \item{id}{the name of the variable that contains the individual index
  if any,}

  \item{opposite}{returns the opposite of the specified variables,}
  
  \item{drop.index}{should the index variables be dropped from the
  \code{data.frame},}
  
  \item{...}{further arguments passed to \code{reshape}.}
}

\value{
  
  A \code{mlogit.data} object, which is a \code{data.frame} in
  \code{long} format, \emph{i.e.} one line for each alternative.  It has
  a \code{index} attribute, which is a \code{data.frame} that contains
  the index of the choice made (\code{'chid'}), the index of the
  alternative (\code{'alt'}) and, if any, the index of the individual
  (\code{'id'}).  The choice variable is a boolean which indicates the
  choice made. This function use \code{reshape} if the \code{data.frame}
  is in \code{wide} format.
  
}

\examples{

# ModeChoice is a long data.frame 

data("TravelMode", package = "AER")
TM <- mlogit.data(TravelMode, choice = "choice", shape = "long",
                 alt.levels = c("air", "train", "bus", "car"))

# Same but the alt variable called mode is provided

TM <- mlogit.data(TravelMode ,choice = "choice", shape = "long",
                  alt.var = "mode")

# Same but the chid variable called individual is provided

TM <- mlogit.data(TravelMode, choice = "choice",
                  shape = "long", id.var = "individual",
                  alt.levels = c("air", "train", "bus", "car"))

# Same but with two own provided variables

TM <- mlogit.data(TravelMode, choice = "choice", shape = "long",
                 id.var = "individual", alt.var = "mode")

# Same but with two own provided variables which are deleted from the
  data.frame

TM <- mlogit.data(TravelMode, choice = "choice", shape = "long",
                 id.var = "individual", alt.var = "mode", drop.index = TRUE)


# Train is a wide data.frame with columns 'choiceid' is the choice
#  index, the alternatives are named "ch1" and "ch2", the opposite of
#  the variables is returned

data("Train", package = "mlogit")
Train <- mlogit.data(Train, choice = "choice", shape = "wide",
                     varying = 4:11, alt.levels = c("ch1", "ch2"), sep = "",
                     opposite = c("price", "time", "change", "comfort"))


# Car is a wide data.frame

data("Car", package = "mlogit")
Car <- mlogit.data(Car, varying = 5:70, shape = "wide", sep = "",
                   choice = "choice", alt.levels = 1:6)

data("HC", package = "mlogit")
HC <- mlogit.data(HC, choice = "depvar", varying=c(2:8, 10:16), shape="wide")

}
  
\author{Yves Croissant}

\seealso{\cite{\code{reshape}}}

\keyword{attribute}
