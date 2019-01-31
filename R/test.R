#' Hausman-McFadden Test
#' 
#' Test the IIA hypothesis (independence of irrelevant alternatives)
#' for a multinomial logit model.
#' 
#' @name hmftest
#' @aliases hmftest hmftest.formula hmftest.mlogit
#' @param x an object of class `mlogit` or a formula,
#' @param z an object of class `mlogit` or a subset of alternatives
#'     for the `mlogit` method. This should be the same model as `x`
#'     estimated on a subset of alternatives,
#' @param alt.subset a subset of alternatives,
#' @param ... further arguments passed to `mlogit` for the `formula`
#'     method.
#' @export
#' @return an object of class `"htest"`.
#'
#' @details
#' This is an implementation of the Hausman's consistency test for
#' multinomial logit models. If the independance of irrelevant
#' alternatives applies, the probability ratio of every two
#' alternatives depends only on the characteristics of these
#' alternatives. Consequentely, the results obtained on the estimation
#' with all the alternatives or only on a subset of them are
#' consistent, but more efficient in the first case. On the contrary,
#' only the results obtained from the estimation on a relevant subset
#' are consistent. To compute this test, one needs a model estimated
#' with all the alternatives and one model estimated on a subset of
#' alternatives. This can be done by providing two objects of class
#' `mlogit`, one object of class `mlogit` and a character vector
#' indicating the subset of alternatives, or a formula and a subset of
#' alternatives.
#'
#' @author Yves Croissant
#' @references
#' 
#' Hausman, J.A. and D. McFadden (1984), A Specification Test for the
#' Multinomial Logit Model, *Econometrica*, **52**, pp.1219--1240.
#' 
#' @keywords htest
#' @examples
#' 
#' ## from Greene's Econometric Analysis p. 731
#' 
#' data("TravelMode",package="AER")
#' TravelMode <- mlogit.data(TravelMode,choice="choice",shape="long",
#'                           alt.var="mode",chid.var="individual")
#' 
#' ## Create a variable of income only for the air mode
#' 
#' TravelMode$avinc <- with(TravelMode,(mode=='air')*income)
#' 
#' ## Estimate the model on all alternatives, with car as the base level
#' ## like in Greene's book.
#' 
#' #x <- mlogit(choice~wait+gcost+avinc,TravelMode,reflevel="car")
#' x <- mlogit(choice~wait+gcost+avinc,TravelMode)
#' 
#' ## Estimate the same model for ground modes only (the variable avinc
#' ## must be dropped because it is 0 for every observation
#' 
#' g <- mlogit(choice~wait+gcost,TravelMode,reflevel="car",
#'             alt.subset=c("car","bus","train"))
#' 
#' ## Compute the test
#' 
#' hmftest(x,g)
hmftest <- function(x,...){
  UseMethod("hmftest")
}

#' @rdname hmftest
#' @method hmftest formula
#' @export
hmftest.formula <- function(x, alt.subset, ...){
  formula <- x
  x <- mlogit(formula,...)
  x$call$data <- match.call()$data
  xs <- mlogit(formula, alt.subset=alt.subset, ...)
  hmftest(x,xs)
}

#' @rdname hmftest
#' @method hmftest mlogit
#' @export
hmftest.mlogit <- function(x, z, ...){
  if (is.character(z)) xs <- update(x,alt.subset=z)
  if (class(z)=="mlogit") xs <- z
  coef.x <- coef(x)
  coef.s <- coef(xs)
  un <- names(coef.x) %in% names(coef.s)
  diff.coef <- coef.s-coef.x[un]
  diff.var <- vcov(xs)-vcov(x)[un,un]
  hmf <- as.numeric(diff.coef%*%solve(diff.var)%*%diff.coef)
  names(hmf) <- "chisq"
  df <- sum(un)
  names(df) <- "df"
  pv <- pchisq(hmf,df=df,lower.tail=FALSE)
  res <- list(data.name = x$call$data,
              statistic = hmf,
              p.value =pv,
              parameter = df,
              method = "Hausman-McFadden test",
              alternative = "IIA is rejected")
  class(res) <- "htest"
  res
}  


mfR2 <- function(x){
##   ll <- logLik(x)
##   data.name <- x$call$data
##   choice.name <- as.character(x$call$formula[[2]])
##   data <- eval(data.name,envir=parent.frame())
##   alt <- data[[2]]
##   choice <- data[[choice.name]]
##   eff <- table(alt[choice])
##   n <- sum(eff)
##   llo <- sum(eff*log(eff/n))
##   1-ll/llo
  logLik0 <- attr(x$logLik, 'null')
  1-x$logLik/logLik0
}

lratio <- function(object){
  freq <- object$freq
  llo <- sum(freq*log(prop.table(freq)))
  data.name <- object$call$data
  stat <- -2*(llo-logLik(object))
  names(stat) <- "chisq"
  parameter <- length(coef(object))-length(freq)+1
  names(parameter) <- "df"
  pval <- pchisq(stat,df=parameter,lower.tail=FALSE)
  lrtest <- list(statistic = stat,
                 data.name = data.name,
                 p.value = pval,
                 parameter = parameter,
                 method = "likelihood ratio test")
  class(lrtest) <- "htest"
  lrtest
}

irrelevant.args.warning <- function(object, args){
    if (any(!(names(object) %in% args))){
        irrelevant.args <- !(names(object) %in% args)
        be <- ifelse(sum(irrelevant.args) > 1, "are", "is")
        irrelevant.args <- paste(names(object)[irrelevant.args], collapse = ", ")
        warning(paste("arguments", irrelevant.args, be, "irrelevant and", be, "ignored", sep = " "))
    }
}

#' The three tests for mlogit models
#' 
#' Three tests for mlogit models: specific methods for the Wald test
#' and the likelihood ration test and a new function for the score
#' test
#' 
#' @name scoretest
#' @importFrom lmtest lrtest lrtest.default waldtest waldtest.default
#' @aliases scoretest scoretest.mlogit scoretest.default
#'     waldtest.mlogit waldtest lrtest.mlogit lrtest
#' @param object an object of class `mlogit` or a formula,
#' @param ... two kinds of arguments can be used. If `mlogit`
#'     arguments are introduced, initial model is updated using these
#'     arguments. If `formula` or other `mlogit` models are
#'     introduced, the standard behavior of [lmtest::waldtest()] and
#'     [lmtest::lrtest()] is followed.
#' @details The `scoretest` function and `mlogit` method for
#'     `waldtest` and `lrtest` from the `lmtest` package provides the
#'     infrastructure to compute the three tests of hypothesis for
#'     `mlogit` objects.
#' 
#' The first argument must be a `mlogit` object. If the second one is a
#' fitted model or a formula, the behaviour of the three functions is the one
#' of the default methods of `waldtest` and `lrtest`: the two
#' models provided should be nested and the hypothesis tested is that the
#' constrained model is the `right' model.
#' 
#' If no second model is provided and if the model provided is the
#' constrained model, some specific arguments of `mlogit` should be
#' provided to descibe how the initial model should be updated. If the
#' first model is the unconstrained model, it is tested versus the
#' `natural' constrained model; for example, if the model is a
#' heteroscedastic logit model, the constrained one is the multinomial
#' logit model.
#' @return an object of class `htest`.
#' @export
#' @author Yves Croissant
#' @keywords htest
#' @examples
#' library("mlogit")
#' library("lmtest")
#' data("TravelMode", package = "AER")
#' ml <- mlogit(choice ~ wait + travel + vcost, TravelMode,
#'              shape = "long", chid.var = "individual", alt.var = "mode")
#' hl <- mlogit(choice ~ wait + travel + vcost, TravelMode,
#'              shape = "long", chid.var = "individual", alt.var = "mode",
#'              method = "bfgs", heterosc = TRUE)
#' lrtest(ml, hl)
#' waldtest(hl)
#' scoretest(ml, heterosc = TRUE)
scoretest <- function(object, ...){
    UseMethod("scoretest")
}

#' @rdname scoretest
#' @export
scoretest.mlogit <- function(object, ...){
    objects <- list(object, ...)
    margs <- c('nests', 'un.nest.el', 'unscaled', 'heterosc', 'rpar',
               'R', 'correlation', 'halton', 'random.nb', 'panel', 'constPar')
    mlogit.args <- objects[names(objects) %in% margs]
    if (!is.null(names(objects))) objects <- objects[!(names(objects) %in% margs)]
    nmodels <- length(objects)
    start.values <- c(coef(object))
    m <- list(nests = NULL, un.nest.el = FALSE, unscaled = FALSE, heterosc = FALSE,
              rpar = NULL, R = 40, correlation = FALSE, halton = NULL,
              random.nb = NULL, panel = FALSE)
    m[names(mlogit.args)] <- mlogit.args
    
    # if several models are provided, just use the default method
    if (nmodels > 1){
        return(scoretest.default(object, ...))
    }
    heterosc.logit <- (m$heterosc)
    nested.logit <- (! is.null(m$nests) || ! is.null(object$nests))
    mixed.logit <- (! is.null(m$rpar) || m$correlation)
    if (heterosc.logit + nested.logit + mixed.logit == 0)
        stop("an unconstrained model should be described")
    if (heterosc.logit + nested.logit + mixed.logit > 1)
        stop("only one unconstrained model should be described")
    if (heterosc.logit){
        alt.hyp <- "heteroscedastic model"
        data.name <- "heterosc = TRUE"
    }
    if (nested.logit){
        init.nested.model <- ! is.null(object$call$nests)
        if (init.nested.model){
            if (is.null(object$call$un.nest.el) || !object$call$un.nest.el){
                stop("irrelevant model for a score test")
            }
            J <- length(object$nests)
            start.values <- c(coef(object), rep(coef(object)[length(coef(object))], J - 1))
            data.name <- "un.nest.el = FALSE"
            alt.hyp <- "unique nest elasticity"
        }
        else{
            alt.hyp <- ifelse(m$un.nest.el, "nested model with a unique nest elasticity",
                              "nested model")
            nest.list <- c()
            for (i in 1:length(m$nests)){
                anest <- paste("c(\'",paste(m$nests[[i]],collapse="\',\'"),"\')", sep="")
                anest <- paste(names(m$nests)[i], " = ", anest, sep = "")
                nest.list <- c(nest.list, anest)
            }
            data.name = paste("nests = list(", paste(nest.list, collapse = ", "), ")", sep = "")
            data.name <- paste(names(m$nests), collapse = ", ")
        }
    }
    if (mixed.logit){
        init.mixed.model <- ! is.null(object$call$rpar)
        if (init.mixed.model){
            if (! is.null(object$call$correlation) && object$call$correlation)
                stop("not a relevant model for a score test")
            alt.hyp <- "uncorrelated random effects"
            data.name <- "correlation = TRUE"
        }
        else{
            if (m$correlation) alt.hyp <- "no correlated random effects"
            else alt.hyp <- "no uncorrelated random effects"
            data.name <- paste(names(m$rpar), paste("\'",as.character(m$rpar),"\'", sep = ""),
                               collapse = ",", sep = "=")
            data.name <- paste("rpar", "(", data.name, ")", sep = "")
        }
        if (init.mixed.model){
            ncoef <- names(coef(object))
            J <- length(object$rpar)
            K <- ncol(model.matrix(object))
            sd <- coef(object)[grep("sd.", ncoef)]
            rd.el <- K + (1:(J * (J + 1) / 2))
            diag.el <- K + cumsum(1:J)
            start.values <- c(start.values[1:K], rep(0, length(rd.el)))
            start.values[diag.el] <- sd
        }
    }
    
    mc <- match.call()
    mc[[1]] <- as.name('update')
    mc[c('iterlim', 'method', 'start', 'print.level')] <- list(0, 'bfgs', start.values, 0)
    newmodel <- eval(mc, parent.frame())
    # gradient used to be a vector, now a matrix (the following ifelse should may be removed
    if (is.matrix(newmodel$gradient))
        gradvect <- apply(newmodel$gradient, 2, sum) else gradvect <- newmodel$gradient
    # fixed coefficients should be removed to compute the statistic
    fixed <- attr(newmodel$coefficients, "fixed")
    stat <- - sum(gradvect[! fixed] * solve(newmodel$hessian[! fixed, ! fixed], gradvect[! fixed]))
    names(stat) <- "chisq"
    df <- c(df = length(coef(newmodel)) - length(coef(object)))
    pval <- pchisq(stat, df = df, lower.tail = FALSE)
    result <- list(statistic = stat,
                   parameter = df,
                   p.value = pval,
                   data.name = data.name,
                   method = "score test",
                   alternative = alt.hyp
                   )
    class(result) <- 'htest'
    result
}

#' @rdname scoretest
#' @export
scoretest.default <- function(object, ...){
    new <- list(...)[[1]]
    cls <- class(object)[1]
    nmodels <- length(new)
    if (! inherits(new, 'formula') & ! inherits(new, cls))
        stop("the updating argument doesn't have a correct class")
    if (inherits(new, cls)){
        ncoefs <- names(coef(new))
        new <- formula(formula(new))
    }
    else ncoefs <- names(coef(update(object, new, iterlim = 0)))
    start <- numeric(length = length(ncoefs))
    names(start) <- ncoefs
    supcoef <- ! ncoefs %in% names(coef(object))
    start[names(coef(object))] <- coef(object)
    newmodel <- update(object, new, start= start, iterlim = 0)
    data.name <- paste(deparse(formula(newmodel)))
    alt.hyp <- "unconstrained model"
    if (is.matrix(newmodel$gradient))
        gradvect <- apply(newmodel$gradient, 2, sum) else gradvect <- newmodel$gradient
    stat <- - sum(gradvect * solve(newmodel$hessian, gradvect))
    names(stat) <- "chisq"
    df <- c(df = length(coef(newmodel)) - length(coef(object)))
    pval <- pchisq(stat, df = df, lower.tail = FALSE)
    result <- list(statistic = stat,
                   parameter = df,
                   p.value = pval,
                   data.name = data.name,
                   method = "score test",
                   alternative = alt.hyp
                   )
    class(result) <- 'htest'
    result
}

#' @rdname scoretest
#' @export
waldtest.mlogit <- function(object, ...){
    objects <- list(object, ...)
    margs <- c('nests', 'un.nest.el', 'unscaled', 'heterosc', 'rpar',
               'R', 'correlation', 'halton', 'random.nb', 'panel')
    mlogit.args <- objects[names(objects) %in% margs]
    if (!is.null(names(objects))) objects <- objects[!(names(objects) %in% margs)]
    nmodels <- length(objects)
    specific.computation <- FALSE
    # if several models are provided, just use the default method
    if (nmodels > 1){
        return(waldtest.default(object, ...))
    }
    
    K <- length(colnames(model.matrix(object)))
    L <- length(object$freq)
    
    # guess the nature of the fitted model
    mixed.logit <- !is.null(object$call$rpar)
    heterosc.logit <- !is.null(object$call$heterosc) && object$call$heterosc
    nested.logit <- !is.null(object$call$nests)

    ## Heteroscedastic logit model
    # the hypothesis is that J-1 parameters = 1
    if (heterosc.logit){
        su <- (K+1):(K+L-1)
        q <- rep(1, length(su))
        hyp <- "homoscedasticity"
    }
    
    ## Nested logit Models
    if (nested.logit){
        J <- length(coef(object)) - K
        # First check whether the fitted model has a unique nest
        # elasticity or not
        if (is.null(object$call$un.nest.el)) un.nest.el <- FALSE
        else un.nest.el <- object$call$un.nest.el
        
        # If the fitted model has a unique nest elasticity, the only
        # relevant test is no nests : mlogit.args should be nests=NULL
        # or nothing. A warning is returned in case of supplementary
        # arguments

        if (un.nest.el){
            if (!is.null(mlogit.args$nests)) stop("the nest argument should be NULL")
            irrelevant.args.warning(mlogit.args, "nests")
            su <- K + 1
            q <- 1
            hyp <- "no nests"
        }
        
        # If the nests elasticities are different, two possible tests :
        # 1. no nests (mlogit.args = (nests = NULL)) or nothing. stop if
        # !is.null(nests) and warning if other arguments than nests.
        # 2. unique nest elasticity (mlogit.args = (un.nest.el =
        # TRUE)). stop if un.nest.el = FALSE and warning if other arguments
        # are provided.
        
        if (! un.nest.el){
            if (! is.null(mlogit.args$nests)) stop("the nest argument should be NULL")
            if (! is.null(mlogit.args$un.nest.el) && mlogit.args$un.nest.el){
                irrelevant.args.warning(mlogit.args, "un.nest.el")
                su <- (K + 1):length(coef(object))
                R <- matrix(0, nrow = length(coef(object)), ncol = length(su) - 1)
                for (i in 1:ncol(R)){
                    R[K + 1, i] <- 1
                    R[K + 1 + i, i] <- -1
                }
                Rb <- crossprod(R, coef(object))
                VRV <- t(R) %*% vcov(object) %*% R
                stat <- as.numeric(crossprod(Rb,solve(VRV, Rb)))
                df <- c(df = length(su) - 1)
                specific.computation <- TRUE
                hyp <- "unique nest elasticity"
            }
            else{
                if (length(mlogit.args) == 0 |
                    ("nests" %in% names(mlogit.args) & is.null(mlogit.args$nests))){
                    irrelevant.args.warning(mlogit.args, "nests")
                    su <- (K+1):(K+J)
                    q <- rep(1, length(su))
                    hyp <- "no nests"
                }
                else{
                    stop("irrelevant constrained model")
                }
            }
        }
    }
    
    ## Mixed logit model
    if (mixed.logit){
        ncoefs <- names(coef(object))
        J <- length(object$rpar)
        # First check whether the random effects are correlated or not
        if (is.null(object$call$correlation)) correlation <- FALSE
        else correlation <- object$call$correlation
        # If the fitted model is uncorrelated, the only relevant test is
        # no random effects, mlogit.args = (rpar = NULL) ; stop if rpar is
        # not NULL and warning if supplementary arguments are provided
        if (! correlation){
            if (!is.null(mlogit.args$rpar)) stop("rpar should be NULL")
            irrelevant.args.warning(mlogit.args, "rpar")
            su <- grep("sd.", ncoefs)
            hyp <- "no random effects"
        }
        else{
         # if the fitted model is correlated, two possible tests :
         # 1. uncorrelated random effects : mlogit.args = (correlation =
         # FALSE), stop if (correlation = TRUE) and warning if
         # supplementary aguments are provided
         # 2. no random effects : mlogit.args = (rpar = NULL), stop if
         # rpar not NULL and a warning if supplementary arguments are
         # provided
            rd.el <- grep("chol.", ncoefs)
            diag.el <- rd.el[cumsum(1:J)]
            if (! is.null(mlogit.args$correlation) && mlogit.args$correlation)
                stop("irrelevant constrained model")
            if (! is.null(mlogit.args$correlation) && !mlogit.args$correlation){
                irrelevant.args.warning(mlogit.args, "correlation")
                su <- rd.el[! (rd.el %in% diag.el)]
                hyp <- "uncorrelated random effects"
            }
            else{
                if (!is.null(mlogit.args$rpar)) stop("rpar should be NULL")
                su <- rd.el
                hyp <- "no random effects"
            }
        }
        q <- rep(0, length(su))
    }
    
    if (! specific.computation){
        if (is.null(q)) wq <- coef(object)[su] else wq <- coef(object)[su] - q
        stat <- as.numeric(crossprod(wq,
                                     crossprod(solve(vcov(object)[su, su]),
                                               wq)))
        df <- c(df = length(su))
    }
    names(stat) <- 'chisq'
    pval <- pchisq(stat, df = df, lower.tail = FALSE)
    result <- list(statistic = stat,
                   parameter = df,
                   p.value = pval,
                   data.name = hyp,
                   method = "Wald test"
 #                 alternative = "unconstrainted model"
                   )
    class(result) <- 'htest'
    result
}

#' @rdname scoretest
#' @export
lrtest.mlogit <- function(object, ...){
    dots <- list(...)
    if (length(dots) == 0){
        model2 <- update(object, heterosc = FALSE, rpar = NULL,
                         start = NULL, nests = NULL,
                         gleontief = FALSE, method = 'nr',
                         constPar = NULL)
        lrtest.default(object, model2)
    }
    else lrtest.default(object, ...)
}

