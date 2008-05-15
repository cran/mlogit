mlogit <- function(formula, data, weights = NULL,...){
  alt.name <- names(data)[[2]]
  formula <- make.formula(formula,alt.name)
  start.time <- proc.time()
  cl <- match.call()
  mf <- model.frame(formula,data)
  y <- model.response(mf)
#  choiceid <- data[["choiceid"]]
  X <- make.X(formula,data,alt.name)
  namesX <- colnames(X)

  if (!is.null(weights)) weights <- data[[weights]]/mean(data[[weights]])
  alt <- data[[alt.name]]
  X <- split(as.data.frame(X),alt)
  X <- lapply(X,as.matrix)
  y <- split(y,alt)
  n <- length(y[[1]])
  K <- ncol(X[[1]])
  f <- function(param) mlogit.lnl(param, X, y, weights = NULL)
  g <- function(param) mlogit.grad(param, X, y, weights = NULL)
  h <- function(param) mlogit.hess(param, X, y, weights = NULL)
#  f.nlm <- function(param) mlogit.nlm(param, X, y, weights =NULL)

  result <- maxLik(f,g,h,start=rep(0,K), method="nr", ...)
#  result2 <- nlm(f.nlm,rep(0,K))
  coef <- result$estimate
  P <- mlogit.P(coef,X)
  fitted.values <- as.matrix(as.data.frame(P))
  residuals <- as.matrix(as.data.frame(y))-fitted.values
  logLik <- result$maximum
  hessian <- result$hessian
  convergence.OK <- result$code<=2
  grad.conv <- result$gradient
  names(coef) <- rownames(hessian) <- colnames(hessian) <- namesX
  elaps.time <- proc.time() - start.time
  nb.iter <- result$iterations
  eps <- grad.conv%*%solve(-result$hessian)%*%grad.conv

  est.stat <- list(elaps.time = elaps.time,
                   nb.iter = nb.iter,
                   eps = eps,
                   method = result$type,
                   message = result$message
                   )
  class(est.stat) <- "est.stat"
  result <- list(coefficients = coef, logLik = logLik,
                 hessian = hessian, gradient = grad.conv,
                 call = cl, est.stat = est.stat,
                 residuals = residuals, fitted.values = fitted.values)
  class(result) <- "mlogit"
  result
}

mlogit.P <- function(param,X){
  eXb <- lapply(X,function(x) exp(crossprod(t(x),param)))
  seXb <- suml(eXb)
  P <- lapply(eXb,function(x){ v <- x/seXb; as.vector(v)})
  P
}
  
mlogit.lnl <- function(param, X, y, weights = NULL){
  if (is.null(weights)) weights <- 1
  P <- mlogit.P(param,X)
  Pch <- suml(mapply("*",P,y,SIMPLIFY=FALSE))
  sum(weights*log(Pch))
}

mlogit.grad <- function(param, X, y, weights = NULL){
  if (is.null(weights)) weights <- 1
  P <- mlogit.P(param,X)
  PX <- suml(mapply("*",X,P,SIMPLIFY=FALSE))
  Xch <- suml(mapply("*",X,y,SIMPLIFY=FALSE))
  weights*(Xch-PX)
}

mlogit.hess <- function(param, X, y, weights = NULL){
  if (is.null(weights)) weights <- 1
  P <- mlogit.P(param,X)
  PX <- suml(mapply("*",X,P,SIMPLIFY=FALSE))
  Pch <- suml(mapply("*",P,y,SIMPLIFY=FALSE))
  Xch <- suml(mapply("*",X,y,SIMPLIFY=FALSE))
  XmPX <- lapply(X,function(x) x-PX)
  -suml( mapply(function(x,y) crossprod(x*y,y),P,XmPX,SIMPLIFY=FALSE))
}


## mlogit.nlm <- function(param, X, y, weights =NULL){
##   if (is.null(weights)) weights <- 1
##   P <- mlogit.P(param,X)
##   Pch <- suml(mapply("*",P,y,SIMPLIFY=FALSE))
##   lnl <- -sum(weights*log(Pch))
##   PX <- suml(mapply("*",X,P,SIMPLIFY=FALSE))
##   Xch <- suml(mapply("*",X,y,SIMPLIFY=FALSE))
##   gradient <- -apply(weights*(Xch-PX),2,sum)
##   XmPX <- lapply(X,function(x) x-PX)
##   hessian <- suml( mapply(function(x,y) crossprod(x*y,y),P,XmPX,SIMPLIFY=FALSE))
##   attr(lnl,"gradient") <- gradient
##   attr(lnl,"hessian") <- hessian
##   lnl
## }


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
  class(object) <- c("summary.mlogit","mlogit")
  return(object)
}

print.summary.mlogit <- function(x,digits= max(3, getOption("digits") - 2),width=getOption("width"),...){

  cat("\nCall:\n")
  print(x$call)

  cat("\n")
  print(x$est.stat)
  
  cat("\nCoefficients :\n")
  printCoefmat(x$CoefTable,digits=digits)
  cat("\n")
  cat(paste("Log-Likelihood: ",signif(x$logLik,digits),"\n",sep=""))
  invisible(x)
}

print.est.stat <- function(x, ...){
  et <- x$elaps.time[3]
  i <- x$nb.iter[1]
  eps <- as.numeric(x$eps)
  s <- round(et,0)
  h <- s%/%3600
  s <- s-3600*h
  m <- s%/%60
  s <- s-60*m
  tstr <- paste(h,"h:",m,"m:",s,"s",sep="")
  cat(paste(x$method,"\n"))
  cat(paste(x$message,"\n"))
  cat(paste(i,"iterations,",tstr,"\n"))
  cat(paste("g'(-H)^-1g =",sprintf("%5.3G",eps),"\n"))
}

make.formula <- function(formula, alt.name){
  ind.spec <- NULL
  rhs <- formula[[3]]
  saveformula <- formula
  if(length(rhs)>1 && rhs[[1]]=="|"){
#    ind.spec <- paste("alt:(",deparse(rhs[[3]]),")",sep="")
    ind.spec <- paste(alt.name,":(",deparse(rhs[[3]]),")",sep="")
    choice.spec <- paste(deparse(rhs[[2]]))
    y <- as.character(formula[[2]])
    choice.spec <- paste(choice.spec,"+",ind.spec)
    formula <- as.formula(paste(y," ~ ",choice.spec,sep=""))
  }
  formula
}


make.X <- function(formula,data,alt.name){
  if (attr(terms(formula),"intercept")==1){
    upform <- as.formula(paste(".~",alt.name,"+.+1"))
    formula <- update(formula,upform)
  }
  else{
    formula <- update(formula,.~.+1)
  }
  mf <- model.frame(formula,data)
  X <- model.matrix(formula,data)[,-1,drop=F]
  y <- model.response(mf)
#  choiceid <- data[["choiceid"]]
  lev1 <- levels(data[[alt.name]])[1]
  motif <- paste(alt.name,lev1,":",sep="")
  varkeep <- regexpr(motif,colnames(X))<0
  X <- X[,varkeep]
  X
}


mlogit.data <- function(x, choice, cvar = NULL, shape = "vert", alt = NULL, chid = NULL, opposite = NULL, ...){
  # alt is whether the alternative variable or a vector of alternatives
  if (is.null(alt)){
    if (shape=="vert") stop("alt shouldn't be NULL when shape equal vert")
    alt.levels <- levels(x[[choice]])
    J <- length(alt.levels)
    alt.is <- "void"
    alt.name <- "alt"
  }
  else{
    if (length(alt)==1){
      alt.name <- alt
      if (!is.factor(x[[alt.name]])) x[[alt.name]] <- factor(x[[alt.name]])
      alt.levels <- levels(x[[alt.name]])
      J <- length(alt.levels)
      alt.is <- "variable"
    }
    else{
      alt.name <- "alt"
      alt.levels <- alt
      J <- length(alt.levels)
      alt.is <- "levels"
    }
  }
  if (shape=="vert") n <- nrow(x)/J else n <- nrow(x)
  
  # chid is the choice variable

  if (is.null(chid)) chid.name <- "chid" else chid.name <- chid
  
  if (!(shape %in% c("hor.var","hor.alt","vert")))
    stop("shape must be one of hor.var, hor.alt or vert")

  horizontal <- substr(shape,1,3) == "hor"
  choice.name <- choice
  choice <- x[[choice]]
  if (horizontal){
    classt <- substr(shape,5,7)
    data <- list()
    cvarindex <- c()
    if (is.null(chid)) choiceid <- rep(as.factor(rownames(x)),each=J) else choiceid <- rep(x[[chid.name]],each=J)
    K <- length(cvar)
    if (is.null(cvar)) stop("there should be at least one choice specific variable")
    if (is.list(cvar)) cvar <- unlist(cvar)
    for (i in 1:length(cvar)){
      value <- cvar[i]
      name <- names(cvar)[i]
      if (classt=="var"){
        varirange <- value:(value+J-1)
      }
      else{
        varirange <- value+K*(0:(J-1))
      }
      v <- x[,varirange]
      v <- as.vector(t(as.matrix(v)))
      data[[name]] <- v
      cvarindex <- c(cvarindex,varirange)
    }
    alt <- factor(rep(alt.levels,n),levels=alt.levels)
    data <- data.frame(data)
    x <- cbind(index=1:n,x[,-cvarindex,drop=FALSE])
    xr <- x
    for (i in 1:(J-1)){
      xr <- rbind(xr,x)
    }
    xr <- xr[order(xr$index),-1,drop=FALSE]
    x <- cbind(xr,data)
    x[[choice.name]] <- x[[choice.name]] == alt
  }
  else{
    if (alt.is=="variable"){
      alt <- x[[alt.name]]
    }
    if (alt.is=="levels"){
      alt <- factor(rep(alt.levels,n),levels=alt.levels)
    }
    if (alt.is=="void"){
      alt <- x[[alt.name]]
    }
    if (is.null(chid)){
      choiceid <- rep(1:n,each=J)
    }
    else{
      choiceid <- x[[chid.name]]
    }
    if (!is.logical(x[[choice.name]])) x[[choice.name]] <- as.logical(x[[choice.name]])
  }
  x <- cbind(ch=choiceid,alt=alt,x)
  if (alt.is=="variable"){
    i <- which(names(x)==alt.name)
    x <- x[-i]
  }

  if (!is.null(chid)){
    i <- which(names(x)==chid.name)
    x <- x[-i]
  }
  names(x)[1:2] <- c(chid.name,alt.name)
    
  rownames(x) <- 1:(n*J)
  if (!is.null(opposite)){
    for (i in opposite){
      x[[i]] <- -x[[i]]
    }
  }
  x
}

suml <- function(x){
  n <- length(x)
  if (!is.null(dim(x[[1]]))){
    d <- dim(x[[1]])
    s <- matrix(0,d[1],d[2])
    for (i in 1:n){
      s <- s+x[[i]]
    }
  }
  else{
    s <- rep(0,length(x[[n]]))
    for (i in 1:n){
      s <- s+x[[i]]
    }
  }
  s
}
