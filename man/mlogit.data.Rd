\name{mlogit.data}
\alias{mlogit.data}
\title{data.frame for logit model}
\description{
  shape a \code{data.frame} in a suitable form for the use of the \code{mlogit} function.
}
\usage{
mlogit.data(data, choice, shape = c("wide","long"), varying = NULL,
            sep=".",alt.var = NULL, id.var = NULL,
            alt.levels = NULL, opposite = NULL, ...)
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
  
  \item{id.var}{the name of the variable that contains the choice index
    or the name under which the choice index will be stored,}

  \item{alt.levels}{the name of the alternatives: if null, for
  a \code{wide} data.frame, they are guessed from the variable names and
  the choice variable (both should be the same), for a \code{long}
  \code{data.frame}, they are guessed from the \code{alt.var} argument,}

  \item{opposite}{returns the opposite of the specified variables,}

  \item{...}{further arguments passed to \code{reshape}.}
}

\value{
  
  A \code{data.frame} in \code{long} format, \emph{i.e.} one line for
  each alternative.  The first two columns contain the choice index and
  the alternative index. The choice variable is a boolean which
  indicates the choice made. This function use \code{reshape} if the
  \code{data.frame} is in \code{wide} format.
  
}

\examples{

# ModeChoice is a long data.frame 

data("TravelMode",package="AER")
TM <- mlogit.data(TravelMode,choice="choice",shape="long",
                 alt.levels=c("air","train","bus","car"))

# Same but the alt variable called mode is provided

TM <- mlogit.data(TravelMode,choice="choice",shape="long",
                  alt.var="mode")

# Same but the chid variable called individual is provided

TM <- mlogit.data(TravelMode,choice="choice",shape="long",
                 id.var="individual",alt.levels=c("air","train","bus","car"))

# Same but with two own provided variables

TM <- mlogit.data(TravelMode,choice="choice",shape="long",
                 id.var="individual",alt.var="mode")


# Train is a wide data.frame with columns 'choiceid' is the choice
#  index, the alternatives are named "ch1" and "ch2", the opposite
# of the variables is returned

data("Train",package="Ecdat")
Train <- mlogit.data(Train,choice="choice",shape="wide",
                     varying=4:11,alt.levels=c("ch1","ch2"),
                     opposite=c("price","time","change","comfort"), sep="")

# Car is a wide data.frame, indexes names are provided

data("Car",package="Ecdat")
Car <-
mlogit.data(Car,varying=5:70,shape="wide",sep="",choice="choice",
            id.var="Mr",alt.var="Car",alt.levels=1:6)

  }

\author{Yves Croissant}

\seealso{\cite{\code{reshape}}}

\keyword{attribute}
