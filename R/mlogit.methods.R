##----------------------------
## methods for mlogit objects |
##----------------------------
##    * fitted                |
##    * residuals             |
##    * df.residual           |
##    * terms                 |
##    * model.matrix          |
##    * model.response        |
##    * update                |
##    * print                 |
##    * vcov                  |
##    * logLik                |
##    * summary               |
##    * print.summary         |
##----------------------------

fitted.mlogit <- function(object, outcome = TRUE, ...){
  if (!outcome){
    result <- object$fitted.values
  }
  else{
    index <- attr(object$model, "index")
    J <- length(levels(index[[2]]))
    y <- matrix(model.response.mlogit(object),ncol=J,byrow=T)
    result <- apply(y*object$fitted.values,1,sum)
  }
  result
}

residuals.mlogit <- function(object, outcome = TRUE, ...){
  if (!outcome){
    result <- object$residuals
  }
  else{
    J <- ncol(object$residuals)
    y <- matrix(model.response(object$model),ncol=J,byrow=T)
    result <- apply(y*object$residuals,1,sum)
  }
  result
}

df.residual.mlogit <- function(object, ...){
  n <- length(residuals(object))
  K <- length(coef(object))
  n-K
}

terms.mlogit <- function(x, ...){
  terms(x$formula)
}

model.matrix.mlogit <- function(object, ...){
  model.matrix(object$formula,object$model)
}

model.response.mlogit <- function(object, ...){
  y.name <- paste(deparse(object$formula[[2]]))
  object$model[[y.name]]
}

update.mlogit <- function (object, new, ...){
  call <- object$call
  if (is.null(call))
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(new))
    call$formula <- update(formula(object), new)
  if(length(extras) > 0) {
    existing <- !is.na(match(names(extras), names(call)))
    ## do these individually to allow NULL to remove entries.
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if(any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  eval(call, parent.frame())
}

print.mlogit <- function (x, digits = max(3, getOption("digits") - 2), width = getOption("width"), ...){
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}

vcov.mlogit <- function(object,...){
  solve(-object$hessian)
}

logLik.mlogit <- function(object,...){
  object$logLik
}

summary.mlogit <- function (object,...){
  b <- coef(object)
  std.err <- sqrt(diag(vcov(object)))
  z <- b/std.err
  p <- 2*(1-pnorm(abs(z)))
  CoefTable <- cbind(b,std.err,z,p)
  colnames(CoefTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
  object$CoefTable <- CoefTable
  if (has.intercept(object$formula)){
    object$lratio <- lratio(object)
    object$mfR2 <- mfR2(object)
  }
  if (!is.null(object$rpar)){
    rpar <- object$rpar
    object$summary.rpar <- t(sapply(rpar,summary))
  }

  class(object) <- c("summary.mlogit","mlogit")
  return(object)
}

print.summary.mlogit <- function(x,digits= max(3, getOption("digits") - 2),width=getOption("width"),...){

  cat("\nCall:\n")
  print(x$call)

  cat("\n")
  
  cat("Frequencies of alternatives:")
  print(prop.table(x$freq),digits=digits)
  
  cat("\n")
  print(x$est.stat)
  
  cat("\nCoefficients :\n")
  printCoefmat(x$CoefTable,digits=digits)
  cat("\n")
  cat(paste("Log-Likelihood: ",signif(x$logLik,digits),"\n",sep=""))

  if (has.intercept(x$formula)){
  
    cat("McFadden R^2: ",signif(x$mfR2,digits),"\n")
  
    cat("Likelihood ratio test : ",names(x$lratio$statistic),
        " = ",signif(x$lratio$statistic,digits),
        " (p.value=",format.pval(x$lratio$p.value,digits=digits),")\n",sep="")
  }

  if (!is.null(x$summary.rpar)){
    cat("\nrandom coefficients\n")
    print(x$summary.rpar)
  }
  invisible(x)
}
