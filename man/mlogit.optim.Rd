\name{mlogit.optim}
\alias{mlogit.optim}
\title{Non-linear minimization routine}
\description{
  This function performs efficiently the optimization of the likelihood
  functions for multinomial logit models
}
\usage{
mlogit.optim(f, start, method = c("bfgs", "nr", "bhhh"), iterlim = 2000,
             tol = 1E-06, print.level = 1, constPar = NULL, ...)
}
\arguments{
  \item{f}{the likelihood function to be maximized,}
  \item{start}{the initial value of the vector of coefficients,}
  \item{method}{the method used, one of \code{'nr'} for Newton-Ralphson,
  \code{'bhhh'} for Berndt-Hausman-Hall-Hall and \code{'bfgs'},}
  \item{iterlim}{the maximum number of iterations,}
  \item{tol}{the value of the criteria for the stoping rule,}
  \item{print.level}{one of (0, 1, 2), the details of the printing
    messages. If \code{'print.level=0'}, no information about the
    optimization process is provided, if \code{'print.level=1'} the
    value of the likelihood, the step and the stoping criteria is
    printing, if \code{'print.level=2'} the vectors of the parameters
    and the gradient are also printed.}
  \item{constPar}{a numeric or a character vector which indicates that
    some parameters should be treated as constant,}
  \item{...}{further arguments passed to \code{f}.}
}

\value{
  a list that contains the followings elements :

  \item{optimum}{the value of the function at the optimum, with
  attributes: \code{gradi} a matrix that contains the contribution of
  each individual to the gradient, \code{gradient} the gradient and, if
  \code{method='nr'} \code{hessian} the hessian,}

  \item{coefficients}{the vector of the parameters at the optimum,}

  \item{est.stat}{a list that contains some information about the
  optimization : \code{'nb.iter'} the number of iterations, \code{'eps'}
  the value of the stoping criteria, \code{'method'} the method of
  optimization method used, \code{'message'}}

}

\details{
  
  The optimization is performed by updating, at each iteration, the
  vector of parameters by the amount step * direction, where step is a
  positive scalar and direction = H^-1 * g, where g is the gradient and
  H^-1 is an estimation of the inverse of the hessian. The choice of
  H^-1 depends on the method chosen :
  
  if \code{method='nr'}, H is the hessian (i.e. is the second derivates
  matrix of the likelihood function),
  
  if \code{method = 'bhhh'}, H is the outer-product of the individual
  contributions of each individual to the gradient,

  if \code{method = 'bfgs'}, H^-1 is updated at each iteration using a
    formula that uses the variations of the vector of parameters and the
    gradient. The initial value of the matrix is the inverse of the
    outer-product of the gradient (i.e. the bhh estimator of the
    hessian).
  
  The initial step is 1 and, if the new value of the function is less
  than the previous value, it is divided by two, until a higher value is
  obtained.
  
  The routine stops when the gradient is sufficiently close to 0. The
  criteria is g * H^-1 * g which is compared to the \code{tol}
  argument. It also may stops if the number of iterations equals
  \code{iterlim}.
  
  The function \code{f} has a \code{initial.value} argument which is the
  initial value of the likelihood. The function is then evaluated a
  first time with a step equals to one. If the value is lower than the
  initial value, the step is divided by two until the likelihood
  increases. The gradient is then computed and the function returns as
  attributes the gradient is the step. This method is more efficient
  than other functions available for R :

  For the \code{optim} and the \code{maxLik} functions, the function and
  the gradient should be provided as separate functions. But, for
  multinomial logit models, both depends on the probabilities which are
  the most time-consuming elements of the model to compute.

  For the \code{nlm} function, the fonction returns the gradient as an
  attribute. The gradient is therefore computed at each iteration, even
  when the function is computed with a step that is unable to increase
  the value of the likelihood.

  Previous versions of \code{mlogit} depended on the \code{'maxLik'}
  package. We kept the same interface, namely the \code{start},
  \code{method}, \code{iterlim}, \code{tol}, \code{print.level} and
  \code{constPar} arguments.

  The default method is \code{'bfgs'}, which is known to perform well,
  even if the likelihood function is not well behaved and the default
  value for \code{print.level=1}, which means moderate printing.
  
  A special default behavior is performed if a simple multinomial logit
  model is estimated. Indeed, for this model, the likelihood function is
  concave, the analytical hessian is simple to write and the
  optimization is straightforward. Therefore, in this case, the default
  method is \code{'nr'} and \code{print.level=0}.
  
}

\author{Yves Croissant}
\keyword{regression}
