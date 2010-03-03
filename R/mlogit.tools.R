num.gradient <- function(f, param, ...){
  m <- match.call(expand.dots = TRUE)
  m[[1]] <- as.name(m[[2]])
  m[[2]] <- NULL
  eps <- 1E-10
  nc <- c()
  for (i in 1:length(param)){
    params <- param
    parami <- param
    params[i] <- params[i] + eps
    parami[i] <- parami[i] - eps
    m$param <- params
    lnls <- eval(m, parent.frame())
    m$param <- parami
    lnli <- eval(m, parent.frame())
    nc <- c(nc, (lnls-lnli)/(2*eps))
  }  
  nc 
}

num.hessian <- function(f, param, ...){
  m <- match.call(expand.dots = TRUE)
  m[[1]] <- as.name(m[[2]])
  m[[2]] <- NULL
  m$gradient <- TRUE
  eps <- 1E-8
  K <- length(param)
  nc <- c()
  for (i in 1:K){
    params <- param
    parami <- param
    params[i] <- params[i] + eps
    parami[i] <- parami[i] - eps
    m$param <- params
    gs <- attr(eval(m, parent.frame()), "gradient")
    m$param <- parami
    gi <- attr(eval(m, parent.frame()), "gradient")
    nc <- c(nc, (gs-gi)/(2*eps))
  }  
  matrix(nc, K, K)
}

print.est.stat <- function(x, ...){
  et <- x$elaps.time[3]
  i <- x$nb.iter[1]
  halton <- x$halton
  method <- x$method
  eps <- as.numeric(x$eps)
  if (!is.null(x$type) && x$type != "simple"){
    R <- x$nb.draws
    cat(paste("Simulated maximum likelihood with", R, "draws\n"))
  }
  s <- round(et,0)
  h <- s%/%3600
  s <- s-3600*h
  m <- s%/%60
  s <- s-60*m
  cat(paste(method, "method\n"))
  tstr <- paste(h, "h:", m, "m:", s, "s", sep="")
  cat(paste(i,"iterations,",tstr,"\n"))
  if (!is.null(halton)) cat("Halton's sequences used\n")
  cat(paste("g'(-H)^-1g =", sprintf("%5.3G", eps),"\n"))
  cat(paste(x$message, "\n"))
}

suml <- function(x){
  n <- length(x)
  if (!is.null(dim(x[[1]]))){
    d <- dim(x[[1]])
    s <- matrix(0,d[1],d[2])
    for (i in 1:n){
      x[[i]][is.na(x[[i]])] <- 0
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

numderiv <- function(f, param, ...){
  m <- match.call()
  m[[1]] <- as.name(m[[2]])
  m[[2]] <- NULL
  eps <- 1E-4
  nc <- c()
  for (i in 1:length(param)){
    params <- param
    parami <- param
    params[i] <- params[i]+eps
    parami[i] <- parami[i]-eps
    m$param <- params
    lnls <- eval(m, parent.frame())
    m$param <- parami
    lnli <- eval(m, parent.frame())
    nc <- c(nc, (lnls-lnli)/(2*eps))
  }
  nc
}

mlogit.optim <- function(f, start,
                         method = c('bfgs', 'nr', 'bhhh'),
                         iterlim = 2000,
                         tol = 1E-06,
                         print.level = 1,
                         constPar = NULL,
                         ...){
  # construct a call for the function
  param <- start
  method <- match.arg(method)
  callT <- match.call(expand.dots = TRUE)
  f <- callT
  optimoptions <- c('iterlim', 'tol', 'method', 'print.level', 'constPar')
  m <- match(optimoptions, names(callT), 0L)
  if (sum(m)) f <- f[-m]
  f[[1]] <- as.name(f[[2]])
  K <- length(param)
  fixed <- rep(FALSE, K)
  if (!is.null(constPar)) fixed[constPar] <- TRUE
  f$gradient <- TRUE
  if (method == 'nr') f$hessian <- TRUE else f$hessian <- FALSE
  f[[2]] <- NULL
  names(f)[2] <- 'param'
  chi2 <- 1E+10
  i <- 0
  # eval a first time the function, the gradient and the hessian

###################################
##  Test the gradient
##   nd <- f
##   nd[["f"]] <- nd[[1]]
##   nd[[1]] <- as.name("numderiv")
##   x <- eval(nd, parent.frame())
##   print(x)
###################################

  x <- eval(f, parent.frame())
  if (print.level > 0)
    cat(paste("Initial value of the function :", as.numeric(x), "\n"))
  g <- attr(x, "gradient")
  if (method == 'nr')   H <- attr(x, "hessian")
  if (method == 'bhhh') H <- crossprod(attr(x, "gradi"))
  if (method == 'bfgs') Hm1 <- solve(crossprod(attr(x, "gradi")))
  d <- rep(0, K)
  while(abs(chi2) > tol && i < iterlim){
    if (method == "bfgs") d[!fixed] <- - as.vector(Hm1[!fixed, !fixed] %*% g[!fixed])
    else d[!fixed] <- - as.vector(solve(H[!fixed, !fixed], g[!fixed]))
    i <- i + 1
    oldx <- x
    oldg <- g
    f$param <- param
    f$direction <- d
    f$initial.value <- x
    x <- eval(f, parent.frame())
    g <- attr(x, "gradient")
    step <- attr(x, "step")
    param <- param + step * d
    if (method == 'nr')   H <- attr(x, "hessian")
    if (method == 'bhhh') H <-  crossprod(attr(x, "gradi"))
    if (method == 'bfgs'){
      incr <- step * d
      y <- g - oldg
      Hm1 <- Hm1 +
        outer( incr, incr) *
          (sum(y * incr) + as.vector( t(y) %*% Hm1 %*% y)) / sum(incr * y)^2 -
            (Hm1 %*% outer(y, incr) + outer(incr, y) %*% Hm1)/ sum(incr * y)
    }
    chi2 <- -  crossprod(d[!fixed], oldg[!fixed])
    if (print.level > 0){
      chaine <- paste("iteration ",i,", step = ",step,
                      ", lnL = ",round(x,8),", chi2 = ",
                      round(chi2,8),"\n",sep="")
    cat(chaine)
    }
    if (print.level > 1){
      resdet <- rbind(param = param, gradient = g)
      print(round(resdet,3))
      cat("--------------------------------------------\n")
    }
  }
  if (i >= iterlim) message = "maximum number of iterations reached"
  else message = "optimum reached"
  if (method == 'bfgs') H <- solve(Hm1)
  rownames(H) <- colnames(H) <-
    names(attr(x, 'gradient')) <- colnames(attr(x, 'gradi')) <- names(param)
  est.stat = structure(list(elaps.time = NULL, nb.iter = i, eps = chi2,
    method = method, message = message), class = 'est.stat')
  result <- list(optimum = x,
                 coefficients = param,
                 est.stat = est.stat
                 )
  result
}
