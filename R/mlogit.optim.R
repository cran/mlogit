#' Non-linear minimization routine
#' 
#' This function performs efficiently the optimization of the
#' likelihood functions for multinomial logit models
#' 
#' @name mlogit.optim
#' @param logLik the likelihood function to be maximized,
#' @param start the initial value of the vector of coefficients,
#' @param method the method used, one of `'nr'` for Newton-Ralphson,
#'     `'bhhh'` for Berndt-Hausman-Hall-Hall and `'bfgs'`,
#' @param iterlim the maximum number of iterations,
#' @param tol the value of the criteria for the gradient,
#' @param ftol the value of the criteria for the function,
#' @param steptol the value of the criteria for the step,
#' @param print.level one of (0, 1, 2), the details of the printing
#'     messages.  If `'print.level = 0'`, no information about the
#'     optimization process is provided, if `'print.level = 1'` the
#'     value of the likelihood, the step and the stoping criteria is
#'     printing, if `'print.level = 2'` the vectors of the parameters
#'     and the gradient are also printed.
#' @param constPar a numeric or a character vector which indicates
#'     that some parameters should be treated as constant,
#' @param ... further arguments passed to `f`.
#'
#' @details
#'
#' The optimization is performed by updating, at each iteration, the
#' vector of parameters by the amount step * direction, where step is
#' a positive scalar and direction = H^-1 * g, where g is the gradient
#' and H^-1 is an estimation of the inverse of the hessian. The choice
#' of H^-1 depends on the method chosen :
#' 
#' if `method = 'nr'`, H is the hessian (*i.e.* is the second
#' derivates matrix of the likelihood function),
#' 
#' if `method = 'bhhh'`, H is the outer-product of the individual
#' contributions of each individual to the gradient,
#' 
#' if `method = 'bfgs'`, H^-1 is updated at each iteration using a
#' formula that uses the variations of the vector of parameters and
#' the gradient. The initial value of the matrix is the inverse of the
#' outer-product of the gradient (i.e. the bhhh estimator of the
#' hessian).
#' 
#' The initial step is 1 and, if the new value of the function is less
#' than the previous value, it is divided by two, until a higher value
#' is obtained.
#' 
#' The routine stops when the gradient is sufficiently close to 0. The
#' criteria is g * H^-1 * g which is compared to the `tol`
#' argument. It also may stops if the number of iterations equals
#' `iterlim`.
#' 
#' The function `f` has a `initial.value` argument which is the
#' initial value of the likelihood. The function is then evaluated a
#' first time with a step equals to one. If the value is lower than
#' the initial value, the step is divided by two until the likelihood
#' increases. The gradient is then computed and the function returns
#' as attributes the gradient is the step.  This method is more
#' efficient than other functions available for `R`:
#' 
#' For the `optim` and the `maxLik` functions, the function and the
#' gradient should be provided as separate functions. But, for
#' multinomial logit models, both depends on the probabilities which
#' are the most time-consuming elements of the model to compute.
#' 
#' For the `nlm` function, the fonction returns the gradient as an
#' attribute. The gradient is therefore computed at each iteration,
#' even when the function is computed with a step that is unable to
#' increase the value of the likelihood.
#' 
#' Previous versions of `mlogit` depended on the `'maxLik'` package.
#' We kept the same interface, namely the `start`, `method`,
#' `iterlim`, `tol`, `print.level` and `constPar` arguments.
#' 
#' The default method is `'bfgs'`, which is known to perform well,
#' even if the likelihood function is not well behaved and the default
#' value for `print.level = 1`, which means moderate printing.
#' 
#' A special default behavior is performed if a simple multinomial
#' logit model is estimated. Indeed, for this model, the likelihood
#' function is concave, the analytical hessian is simple to write and
#' the optimization is straightforward. Therefore, in this case, the
#' default method is `'nr'` and `print.level = 0`.
#' 
#' @return
#' 
#' a list that contains the followings elements :
#'
#' - optimum: the value of the function at the optimum, with
#' attributes: `gradi` a matrix that contains the contribution of each
#' individual to the gradient, `gradient` the gradient and, if `method
#' = 'nr', `hessian` the hessian,
#' - coefficients: the vector of the parameters at the optimum,
#' - est.stat: a list that contains some information about the
#' optimization : `'nb.iter'` the number of iterations, `'eps'` the
#' value of the stoping criteria, `'method'` the method of
#' optimization method used, `'message'
#' 
#' @author Yves Croissant
#' @keywords regression
#' @export
mlogit.optim <- function(logLik, start,
                         method = c('bfgs', 'nr', 'bhhh'),
                         iterlim = 2000,
                         tol = 1E-06,
                         ftol = 1E-08,
                         steptol = 1E-10,
                         print.level = 0,
                         constPar = NULL,
                         ...){

    method <- match.arg(method)
    param <- start 
    callT <- match.call(expand.dots = TRUE)
    optimoptions <- c('iterlim', 'tol', 'method', 'print.level', 'constPar', 'ftol', 'steptol')
    chi2 <- 1E+10
    i <- 0
    K <- length(param)
    d <- rep(0, K)
    
    # construct a vector of fixed parameters
    fixed <- rep(FALSE, K)
    names(fixed) <- names(start)
    if (! is.null(constPar)) fixed[constPar] <- TRUE
    # construct a call for the function
    f <- callT
    # if the model is updated based on a multinomial logit model, change the method to bfgs
    # if (method == 'nr' && as.character(f[[1]]) != "lnl.mlogits") method <- 'bfgs'
    m <- match(optimoptions, names(callT), 0L)
    if (sum(m)) f <- f[-m]
    f[[1]] <- as.name(f[[2]])
    f$gradient <- TRUE
    #  f$steptol <- steptol
    f$stptol <- steptol
    if (method == 'nr') f$hessian <- TRUE else f$hessian <- FALSE
    f[[2]] <- NULL
    names(f)[2] <- 'param'
    # eval a first time the function, the gradient and the hessian
    x <- eval(f, parent.frame())
    # set to TRUE to check the analytical gradient
    if (FALSE){
        nd <- f
        nd[["f"]] <- nd[[1]]
        nd[[1]] <- as.name("numderiv")
        nd$gradient <- FALSE
        bla <- eval(nd, parent.frame())
        print(cbind(coef = f$param, numeric = bla, analytic = attr(x, "gradient")))
        cat("____________\n");
    }

    if (print.level > 0)
        cat(paste("Initial value of the function :", as.numeric(x), "\n"))
    g <- attr(x, "gradient")

    if (method == 'nr')   H <- attr(x, "hessian")[! fixed, ! fixed]
    if (method == 'bhhh') H <- crossprod(attr(x, "gradi")[, ! fixed])
    if (method == 'bfgs') Hm1 <- solve(crossprod(attr(x, "gradi")[, ! fixed]))

    repeat{
        # save the previous values of the function and of the gradient
        oldx <- x
        oldg <- g

        # Compute the direction, ie d = H^-1 g
    
        # For the predict method, I don't want the solve
        if (iterlim > 0){
            if (method == "bfgs") d[! fixed] <- - as.vector(Hm1 %*% g[! fixed])
            else d[! fixed] <- - as.vector(solve(H, g[! fixed]))

        }
        i <- i + 1
        
        if (i > iterlim){
            # exit if the iteration limit is reached
            code <- 4
            break
        }

        # indicate in the call the previous parameters vector, the
        # direction and the value of the function
        f$param <- param
        f$direction <- d
        f$initial.value <- oldx
        # eval the function and compute the gradient and the hessian
        x <- eval(f, parent.frame())
        if (is.null(x)){
             # x is null if steptol is reached
            code = 3
            break
        }
        if (abs(x - oldx) < ftol){
            code = 2
            break
        }
        g <- attr(x, "gradient")
        step <- attr(x, "step")
        param[!fixed] <- param[!fixed] + step * d[!fixed]

        if (method == 'nr')   H <- attr(x, "hessian")[!fixed, !fixed]
        if (method == 'bhhh') H <-  crossprod(attr(x, "gradi")[, !fixed])
        if (method == 'bfgs'){
            incr <- step * d
            y <- g - oldg
            Hm1 <- Hm1 +
                outer( incr[!fixed], incr[!fixed]) *
                (sum(y[!fixed] * incr[!fixed]) +
                 as.vector( t(y[!fixed]) %*% Hm1 %*% y[!fixed])) /
                sum(incr[!fixed] * y[!fixed])^2 -
                (Hm1 %*% outer(y[!fixed], incr[!fixed])
                    + outer(incr[!fixed], y[!fixed]) %*% Hm1)/
                sum(incr[!fixed] * y[!fixed])
        }
        # compute the quadratic form of the gradient
        chi2 <- -  crossprod(d[!fixed], oldg[!fixed])

        # print some informations about the iteration
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
        
        if (abs(chi2) < tol){
            # exit if the quadratic form of the gradient is small enough
            code = 1
            break
        }
        
    }
    if (code == 3) x <- oldx
    names(attr(x, 'gradient')) <- colnames(attr(x, 'gradi')) <- names(param)
    attr(x, "fixed") <- fixed
    est.stat = structure(list(elaps.time = NULL, nb.iter = i, eps = chi2,
                              method = method, code = code), class = 'est.stat')
    result <- list(optimum = x,
                   coefficients = param,
                   est.stat = est.stat
                   )
    result
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

print.est.stat <- function(x, ...){
    et <- x$elaps.time[3]
    i <- x$nb.iter[1]
    halton <- x$halton
    method <- x$method
    if (!is.null(x$type) && x$type != "simple"){
        R <- x$nb.draws
        cat(paste("Simulated maximum likelihood with", R, "draws\n"))
    }
    s <- round(et,0)
    h <- s %/% 3600
    s <- s - 3600 * h
    m <- s %/% 60
    s <- s - 60 * m
    cat(paste(method, "method\n"))
    tstr <- paste(h, "h:", m, "m:", s, "s", sep="")
    cat(paste(i,"iterations,",tstr,"\n"))
    if (!is.null(halton)) cat("Halton's sequences used\n")
    if (!is.null(x$eps)) cat(paste("g'(-H)^-1g =", sprintf("%5.3G", as.numeric(x$eps)),"\n"))
    if (is.numeric(x$code)){
        msg <- switch(x$code,
                      "1" = "gradient close to zero",
                      "2" = "successive function values within tolerance limits",
                      "3" = "last step couldn't find higher value",
                      "4" = "iteration limit exceeded"
                      )
        cat(paste(msg, "\n"))
    }
    else cat(paste(x$code, "\n"))
}
