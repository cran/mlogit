\name{mlogit.data}
\alias{mlogit.data}
\title{data.frame for logit model}
\description{
  shape a data.frame in a suitable form for the use of the \code{mlogit} function.
}
\usage{
mlogit.data(x, choice, cvar = NULL, shape = "vert",
            alt = NULL, chid = NULL, opposite = NULL, ...)
}

\arguments{

  \item{x}{a \code{data.frame},}

  \item{choice}{the variable indicating the choice made,}
  
  \item{cvar}{a named vector indicating in which columns start a
    choice-specific variable if data are stored ``vertically'',}

  \item{shape}{if equal to \code{hor.var} or \code{hor.alt}, there is
    one row for each choice and each column contains a variable for an
    alternative (the columns are ordered by variable in the first case
    and by alternative in the second case). If \code{vert}, there is one
    row for each alternative,}

  \item{alt}{if shape is equal to \code{hor.var} or \code{hor.alt},
    it may be \code{NULL}. In this case, the levels of the choice
    variables are assumed to be in the same order as the subsequent
    columns representing the choice-specific variables. Otherwise, they
    can be indicated in the proper order. If shape is \code{vert}, it
    can be whether the name of the variable containing the alternative
    index, or a vector containing the names of the alternative in the
    proper order,}
  \item{chid}{the choice index, only relevant if \code{shape='vert'},}
  \item{opposite}{returns the opposite of the specified variables,}
  \item{...}{further arguments.}
}

\value{
  
  A \code{data.frame} with a line for every alternative.
  
}
\details{

  The first column is a variable that identifies the choice, whether
  auto-generated and called \code{chid} or provided with the \code{chid}
  argument. The second column identifies the alternatives. It is
  provided with the \code{alt} argument or is auto-generated. In this
  later case, its name is \code{alt}.
}

\examples{

# ModeChoice is a vertical shaped data.frame 

data("ModeChoice",package="Ecdat")
ModeChoice <- mlogit.data(ModeChoice,choice="mode",shape="vert",
                 alt=c("air","train","bus","car"))

# Same but with an alt variable called altm

data("ModeChoice",package="Ecdat")
ModeChoice$altm <- rep(c("air","train","bus","car"),210)
ModeChoice <- mlogit.data(ModeChoice,choice="mode",shape="vert",
                 alt="altm")

# Same but with an chid variable called ind 

data("ModeChoice",package="Ecdat")
ModeChoice$ind <- rep(1:210,each=4)
ModeChoice <- mlogit.data(ModeChoice,choice="mode",shape="vert",
                 chid="ind",alt=c("air","train","bus","car"))

# Same but with two own provided variables

data("ModeChoice",package="Ecdat")
ModeChoice$ind <- rep(1:210,each=4)
ModeChoice$altm <- rep(c("air","train","bus","car"),210)
ModeChoice <- mlogit.data(ModeChoice,choice="mode",shape="vert",
                 chid="ind",alt="altm")


# Train is an horizontal shaped data.frame with columns ordered by
# alternatives, 'choiceid' is the choice index

data("Train",package="Ecdat")
Train <- mlogit.data(Train,choice="choice",shape="hor.alt",
                 cvar=c(price=4,time=5,change=6,comfort=7),
                 opposite=c("price","time","change","comfort"),chid="choiceid")

# Car is an horizontal data.frame with columns ordered by variables

data("Car",package="Ecdat")
Car <- mlogit.data(Car,cvar=c(type=5,fuel=11,price=17,range=23,acc=29,speed=35,
                 pollution=41,size=47,space=53,cost=59,station=65),
                 shape="hor.var",choice="choice")

  }

\author{Yves Croissant}

\keyword{attribute}
