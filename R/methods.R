## #' @rdname miscmethods.mlogit
## #' @export
## coef.mlogit <- function(object,
##                         subset = c("all", "iv", "sig", "sd", "sp", "chol"),
##                         fixed = FALSE, ...){
##     whichcoef <- match.arg(subset)
##     result <- object$coefficients
##     ncoefs <- names(result)
##     # first remove the fixed coefficients if required
##     if (! fixed) result <- result[! attr(result, "fixed")]
##     attr(result, "fixed") <- NULL
##     if (whichcoef == "all") selcoef <- 1:length(result)
##     else selcoef <- grep(whichcoef, ncoefs)
##     result[selcoef]
## }

## #' @rdname miscmethods.mlogit
## #' @method tidy mlogit
## #' @export
## tidy.mlogit <- function (x, conf.int = FALSE, conf.level = 0.95, ...) 
## {
##     check_ellipses("exponentiate", "tidy", "mlogit", ...)
##     s <- summary(x)
##     ret <- as_tidy_tibble(s$CoefTable, new_names = c("estimate", 
##         "std.error", "statistic", "p.value"))
##     if (conf.int) {
##         ci <- broom_confint_terms(x, level = conf.level)
##         ret <- dplyr::left_join(ret, ci, by = "term")
##     }
##     ret
## }

## #' @importFrom generics tidy
## #' @export
## generics::tidy


#coef.mlogit <- micsr:::coef.micsr
#select_coef <- micsr:::select_coef

## #' @rdname miscmethods.mlogit
## #' @export
## df.residual.mlogit <- function(object, ...){
##     n <- length(residuals(object))
##     K <- length(coef(object))
##     n - K
## }
## #' @param confint,conflevel see tidy
## #' @param subset an optional vector of coefficients to extract for the
## #'     `coef` method,
## #' @param fixed if `FALSE` (the default), constant coefficients are
## #'     not returned,

# !!! on vire model.response.mlogit de @aliases ci-dessous


## #' Methods for mlogit objects
## #'
## #' Miscellaneous methods for `mlogit` objects.
## #' 
## #' @name miscmethods.mlogit
## #' @aliases residuals.mlogit df.residual.mlogit terms.mlogit
## #'     model.matrix.mlogit update.mlogit
## #'     print.mlogit summary.mlogit print.summary.mlogit predict.mlogit
## #'     fitted.mlogit coef.mlogit coef.summary.mlogit
## #' @param x,object an object of class `mlogit`
## #' @param digits the number of digits,
## #' @param width the width of the printing,
## #' @param new an updated formula for the `update` method,
## #' @param outcome a boolean which indicates, for the `fitted` and the
## #'     `residuals` methods whether a matrix (for each choice, one
## #'     value for each alternative) or a vector (for each choice, only
## #'     a value for the alternative chosen) should be returned,
## #' @param type one of `outcome` (probability of the chosen
## #'     alternative), `probabilities` (probabilities for all the
## #'     alternatives), `parameters` for individual-level random
## #'     parameters for the fitted method, how the correlated random
## #'     parameters should be displayed : `"chol"` for the estimated
## #'     parameters (the elements of the Cholesky decomposition matrix),
## #'     `"cov"` for the covariance matrix and `"cor"` for the
## #'     correlation matrix and the standard deviations,
## #' @param n,m see [dfidx::idx()]
## #' @param ... further arguments.
## NULL

 
#' @rdname miscmethods.mlogit
#' @export
residuals.mlogit <- function(object, outcome = TRUE, ...){
    if (! outcome){
        result <- object$residuals
    }
    else{
        J <- ncol(object$residuals)
        y <- matrix(model.response(object$model), ncol = J, byrow = T)
        result <- apply(y * object$residuals, 1, sum)
    }
    result
}

#' @rdname miscmethods.mlogit
#' @export
terms.mlogit <- function(x, ...){
    terms(x$formula)
}

#' @rdname miscmethods.mlogit
#' @export
model.matrix.mlogit <- function(object, ...){
    model.matrix(object$model)
}

## #' @rdname miscmethods.mlogit
## #' @method model.response mlogit
## #' @export
## model.response.mlogit <- function(object, ...){
##     y.name <- paste(deparse(object$formula[[2]]))
##     object$model[[y.name]]
## }

#' @rdname miscmethods.mlogit
#' @export
update.mlogit <- function (object, new, ...){
    call <- object$call
    if (is.null(call))
        stop("need an object with call component")
    extras <- match.call(expand.dots = FALSE)$...
    if (! missing(new))
        call$formula <- update(formula(object), new)
    if(length(extras) > 0) {
        existing <- ! is.na(match(names(extras), names(call)))
        ## do these individually to allow NULL to remove entries.
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if(any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }
    eval(call, parent.frame())
}

#' @rdname miscmethods.mlogit
#' @export
print.mlogit <- function (x, digits = max(3, getOption("digits") - 2),
                          width = getOption("width"), ...){
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

## #' @rdname miscmethods.mlogit
## #' @export
## logLik.mlogit <- function(object,...){
##     object$logLik
## }

#' @rdname miscmethods.mlogit
#' @export
summary.mlogit <- function (object, ..., type = c("chol", "cov", "cor")){
    type <- match.arg(type)
    fixed <- attr(object$coefficients, "fixed")
    #    b <- coef(object)[! fixed]
    b <- coef(object)
    std.err <- sqrt(diag(vcov(object)))
    z <- b / std.err
    p <- 2 * (1 - pnorm(abs(z)))
    CoefTable <- cbind(b, std.err, z, p)
    colnames(CoefTable) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
    if (type != "chol"){
        sumvcov <- summary(vcov(object, what = "rpar", type = type))
        CoefTable[grep("chol.", rownames(CoefTable)), ] <- sumvcov
        rownames(CoefTable)[grep("chol.", rownames(CoefTable))] <- rownames(sumvcov)
    }
    object$CoefTable <- CoefTable
    if (has.intercept(object)){
        object$lratio <- lratio(object)
        object$mfR2 <- mfR2(object)
    }
    if (! is.null(object$rpar)){
        rpar <- object$rpar
        object$summary.rpar <- t(sapply(rpar, summary))
    }
    class(object) <- c("summary.mlogit", "mlogit")
    return(object)
}

#' @rdname miscmethods.mlogit
#' @method print summary.mlogit
#' @export
print.summary.mlogit <- function(x, digits = max(3, getOption("digits") - 2),
                                 width = getOption("width"), ...){
    cat("\nCall:\n")
    print(x$call)
    cat("\n")
    cat("Frequencies of alternatives:")
    print(prop.table(x$freq), digits = digits)
    cat("\n")
    print(x$est.stat)
    cat("\nCoefficients :\n")
    printCoefmat(x$CoefTable, digits = digits)
    cat("\n")
    cat(paste("Log-Likelihood: ", signif(x$logLik[1], digits), "\n", sep = ""))
    if (has.intercept(x)){
        cat("McFadden R^2: ", signif(x$mfR2, digits), "\n")
        cat("Likelihood ratio test : ", names(x$lratio$statistic),
            " = ", signif(x$lratio$statistic, digits),
            " (p.value = ", format.pval(x$lratio$p.value, digits = digits), ")\n", sep = "")
    }
    if ( !is.null(x$summary.rpar)){
        cat("\nrandom coefficients\n")
        print(x$summary.rpar)
    }
    invisible(x)
}

#' @rdname miscmethods.mlogit
#' @method idx mlogit
#' @export
idx.mlogit <- function(x, n = NULL, m = NULL){
    idx(model.frame(x), n = n, m = m)
}

#' @rdname miscmethods.mlogit
#' @method idx_name mlogit
#' @export
idx_name.mlogit <- function(x, n = NULL, m = NULL){
    idx_name(model.frame(x), n = n, m = m)
}


#' @rdname miscmethods.mlogit
#' @export
fitted.mlogit <- function(object, type = c("outcome", "probabilities",
                                           "linpred", "parameters"),
                          outcome = NULL, ...){
    if (! is.null(outcome)){
        if (outcome) result <- object$fitted
        else result <- object$probabilities
    }
    else{
        type <- match.arg(type)
        result <- switch(type,
                        outcome = object$fitted,
                        probabilities = object$probabilities,
                        linpred = object$linpred,
                        parameters = object$indpar)
    }
    result
}

    
#' @rdname miscmethods.mlogit
#' @method coef summary.mlogit
#' @export
coef.summary.mlogit <- function(object, ...){
    result <- object$CoefTable
    result
}

## #' Fish <- mlogit.data(Fishing, varying = c(2:9), shape = "wide", choice = "mode")



#' vcov method for mlogit objects
#' 
#' The `vcov` method for `mlogit` objects extract the covariance
#' matrix of the coefficients, the errors or the random parameters.
#' 
#' This new interface replaces the `cor.mlogit` and `cov.mlogit`
#' functions which are deprecated.
#'
#' @name vcov.mlogit
#' 
#' @aliases vcov.mlogit print.vcov.mlogit summary.vcov.mlogit
#' print.summary.vcov.mlogit formula.mlogit
#' @param object a `mlogit` object (and a `vcov.mlogit` for the
#' summary method),
#' @param x a `vcov.mlogit` or a `summary.vcov.mlogit` object,
#' @param what indicates which covariance matrix has to be extracted : the
#' default value is `coefficients`, in this case, `vcov` behaves as
#' usual. If `what` equals `errors` the covariance matrix of the
#' errors of the model is returned. Finally, if `what` equals `rpar`,
#' the covariance matrix of the random parameters are extracted,
#' @param type with this argument, the covariance matrix may be returned (the
#' default) ; the correlation matrix with the standard deviation on the
#' diagonal may also be extracted,
#' @param reflevel relevent for the extraction of the errors of a multinomial
#' probit model ; in this case the covariance matrix is of error differences is
#' returned and, with this argument, the alternative used for differentiation
#' is indicated,
#' @param vcov,subset,fixed,grep,invert see the `micsr` method
#' @param digits the number of digits,
#' @param width the width of the printing,
#' @param ... further arguments.
#' @export
#' @author Yves Croissant
#' @seealso  [mlogit()] for the estimation of multinomial logit
#' models.
#' @keywords regression
## vcov.mlogit <- function(object,
##                         what = c('coefficient', 'errors', 'rpar'),
##                         subset = c("all", "iv", "sig", "sd", "sp", "chol", "covariates", "vcov"),
##                         type = c('cov', 'cor', 'sd'),
##                         reflevel = NULL, ...){
##     whichcoef <- match.arg(subset)
##     what <- match.arg(what)
##     type <- match.arg(type)
##     fixed <- attr(object$coefficients, "fixed")
##     ncoefs <- names(object$coefficients)
vcov.mlogit <- function(object,
                        ...,
                        vcov = NULL,
                        what = c('coefficient', 'errors', 'rpar'),
                        subset = NA,
                        fixed = FALSE,
                        grep = NULL,
                        invert = TRUE,
                        type = c('cov', 'cor', 'sd'),
                        reflevel = NULL){
    what <- match.arg(what)
    type <- match.arg(type)
#    fixed <- attr(object$coefficients, "fixed")
    ncoefs <- names(object$coefficients)

    # for the coefficients, we have to check the problem for fixed
    # coefficients
    if (what == 'coefficient'){
        oclass <- class(object)
        class(object) <- setdiff(oclass, "mlogit")
        result <- vcov(object, ..., vcov = vcov, subset = subset, fixed = fixed,
                       grep = grep, invert = invert)
        class(object) <- oclass
    }
    
    if (what == 'errors'){
        if (! is.null(object$omega)){
            if (is.null(reflevel)){
                if (is.list(object$omega)) result <- object$omega[[1]]
                else result <- object$omega
            }
            else result <- object$omega[[reflevel]]
        }
        result <- switch(type,
                         cov = result,
                         cor = result / tcrossprod(sqrt(diag(result))),
                         sd = sqrt(diag(result))
                         )
    }
    if (what == 'rpar'){
        if (is.null(object$rpar)) stop('no random parameters')
        nrpar <- names(object$rpar)
        if (is.null(attr(object$rpar, "covariance"))){
            # No correlated parameters
            result <- stdev(object)
            if (type != 'sd'){
                V  <- matrix(0, length(result), length(result),
                             dimnames = list(names(result), names(result)))
                if (type == 'cor') diag(V) <- 1
#                if (type == 'vcov') diag(V) <- result ^ 2
#                if (type == 'cov') V <- result ^ 2
                if (type == 'cov') diag(V) <- result ^ 2
                result <- V
            }
        }
        else{
            # correlated parameters
            # NEW_COEF
#            coefs <- coef(object, subset = "chol")
            chol_coefs <- grep("chol", names(coef(object)))
            coefs <- coef(object)[chol_coefs]
            ncoefs <- names(coefs)
            # compute the vcov matrix of random parameters
            result <- tcrossprod(ltm(coefs, to = "ltm"))
            # compute the variance of the vcov matrix of random
            # parameters
            vcovstruct <- chol2vcov(object, type = type)
            if (type == "cov") attr(result, "cov") <- vcovstruct
            if (type == 'cor'){
                sd <- sqrt(diag(result))
                result <- result / tcrossprod(sqrt(diag(result)))
                diag(result) <- sd
                attr(result, "cov") <- vcovstruct
            }
            attr(result, "type") <- type
            ## if (type == 'cov'){
            ##     result <- diag(result)
            ##     attr(result, "cov") <- diag(ltm(vcovstruct, to = "ltm"))
            ## }
            if (type == 'sd') result <- sqrt(diag(result))
            nrparcor <- rownames(attr(object$rpar, "covariance"))
            if (is.null(dim(result))) names(result) <- nrparcor
            else dimnames(result) <- list(nrparcor, nrparcor)
        }
    }
    structure(result, class = c('vcov.mlogit', class(result)), type = type)
}

## #' @rdname miscmethods.mlogit
## #' @method formula mlogit
## #' @export
## formula.mlogit <- function(x, ...){
##     x$formula
## }


#' @rdname vcov.mlogit
#' @method print vcov.mlogit
#' @export
print.vcov.mlogit <- function(x, ...){
    attr(x, "cov") <- attr(x, "type") <- NULL
    print(unclass(x))
}

#' @rdname vcov.mlogit
#' @method summary vcov.mlogit
#' @export
summary.vcov.mlogit <- function(object, ...){
    if (is.null(attr(object, "cov")))
        stop("summary.vcov.mlogit only implemented for random parameters")
    if (is.matrix(object)){
        coefs <- ltm(object, to = "vec")
        nrpar <- rownames(object)
        K <- length(nrpar)
        type <- attr(object, "type")
        diag <- ifelse(type == "cov", "var", "sd")
        nstruct <- names_rpar(nrpar, prefix = type, diag = diag, unique = TRUE)
    }
    else{
        coefs <- object
        nstruct <- names(coefs)
    }
    std.err <- sqrt(attr(object, "cov"))
    b <- coefs
    z <- b / std.err
    p <- 2 * pnorm(abs(z), lower.tail = FALSE)
    # construct the object of class summary.plm
    coefficients <- cbind("Estimate"   = b,
                          "Std. Error" = std.err,
                          "z-value"    = z,
                          "Pr(>|z|)"   = p)
    rownames(coefficients) <- nstruct
    if (is.matrix(object)){
        diagpos <- (1:K) * ( (1:K) + 1) / 2
        coefficients <- coefficients[c(diagpos, (1:nrow(coefficients))[- diagpos]), ]
    }
    structure(coefficients, class = "summary.vcov.mlogit")
    }   

#' @rdname vcov.mlogit
#' @method print summary.vcov.mlogit
#' @export
print.summary.vcov.mlogit <- function(x, digits = max(3, getOption("digits") - 2),
                                      width = getOption("width"), ...){
    printCoefmat(x, digits = digits)
}

chol2vcov <- function(x, type = c("cov", "cor")){
    type <- match.arg(type)
    # Take a mlogit object as argument and returns a vector of
    # variance for the structural parameters
    # First get the position of the coefficients of the Cholesky
    # decomposition
    cholspos <- grep("chol.", names(coef(x)))
    # Then get these coefficients
    coefs <- coef(x)[cholspos]
    # and the covariance matrix of these coefficients
    vcovs <- vcov(x)[cholspos, cholspos]
    # compute the matrix of derivatives
    Dcholcov <- function(x){
        # x is a Cholesky matrix (lower triangular + diagonal),
        # entered as a matrix or as a vector ; Dchol returns the
        # matrix of derivatives of the structural parameters (variance
        # and covariance) respective to the estimated parameters (the
        # element of the Cholesky decomposition).
        if (! is.matrix(x)) x <- ltm(x, to = "ltm")
        K <- nrow(x)    
        Delta <- matrix(0, nrow = K * (K + 1) / 2, ncol = K * (K + 1) / 2)
        dims <- c(0, (1:K) * (1:K + 1) / 2)
        Delta[1, 1] <- x[1, 1]
        if (K > 1){
            for (i in 2:K){
                pos <- (dims[i] + 1):dims[i + 1]
                betas <- ltm(ltm(x, to = "vec")[1:dims[i + 1]], to = "ltm")
                Delta[pos, pos] <- betas
                for (j in 1:(i - 1)){
                    Delta[dims[i] + j, (dims[j] + 1):dims[j + 1]] <- x[i, 1:j]
                }
            }
        }
        dblrows <- (1:K) * ( (1:K) + 1) / 2
        Delta[dblrows, ] <- Delta[dblrows, ] * 2
        Delta
    }
    Dcovcor <- function(x){
        y <- ltm(x, to = "ltm")
        sd <- sqrt(diag(y))
        y <- y / outer(sd, sd)
        diag(y) <- sd
        x <- ltm(x, to = "vec")
        y <- ltm(y, to = "vec")
        dims <- length(x)
        K <- - 0.5 + sqrt(1 + 8 * dims) / 2
        diags <- (1:K) * ((1:K) + 1) / 2
        rows <- Reduce("c", lapply(1:K, function(x) 1:x))
        cols <- rep(1:K, 1:K)
        M <- matrix(0, dims, dims)
        for (i in 1:dims){
            if (cols[i] == rows[i]) M[i, i] <- 1 / 2 * y[i] / x[i]
            else{
                M[i, i] <- 1  * y[i] / x[i]
                first <- rows[i]
                second <- cols[i]      
                M[i, rows == first  & cols == first ] <- - 1 / 2  * y[i] / x[i]
                M[i, rows == second & cols == second] <- - 1 / 2  * y[i] / x[i]
            }   
        }
        M
    }
    Derchols <- Dcholcov(ltm(coefs, to = "ltm"))
    # estimate the covariance matrix of the structural parameters
    result <- Derchols %*% vcovs %*% t(Derchols)
    if (type == "cor"){
        coefs <- tcrossprod(ltm(coefs, to = "ltm"))
        Dercov <- Dcovcor(ltm(coefs, to = "ltm"))
        result <- Dercov %*% result %*% t(Dercov)
    }
    # coerce it to a vector and set the relevant names
    result <- diag(result)
    names(result) <- names(coef(x))[cholspos]
    result
}

#' Compute the log-sum or inclusive value/utility
#' 
#' The `logsum` function computes the inclusive value, or inclusive
#' utility, which is used to compute the surplus and to estimate the two steps
#' nested logit model.
#' 
#' @name logsum
#' @param coef a numerical vector or a `mlogit` object, from which the
#' `coef` vector is extracted,
#' @param X a matrix or a `mlogit` object from which the
#' `model.matrix` is extracted,
#' @param formula a formula or a `mlogit` object from which the
#' `formula` is extracted,
#' @param data a `data.frame` or a `mlogit` object from which the
#' `model.frame` is extracted,
#' @param type either `"group"` or `"global"` : if a `group`
#' argument has been provided in the `mlogit.data`, the inclusive values
#' are by default computed for every group, otherwise, a unique global
#' inclusive value is computed for each choice situation,
#' @param output the shape of the results: if `"chid"`, the results is a
#' vector (if `type = "global"`) or a matrix (if `type = "region"`)
#' with row number equal to the number of choice situation, if `"obs"` a
#' vector of length equal to the number of lines of the data in long format is
#' returned.
#' @details
#' The inclusive value, or inclusive utility, or log-sum is the log of the
#' denominator of the probabilities of the multinomial logit model. If a
#' `"group"` variable is provided in the `"mlogit.data"` function,
#' the denominator can either be the one of the multinomial model or those of
#' the lower model of the nested logit model.
#' 
#' If only one argument (`coef`) is provided, it should a `mlogit`
#' object and in this case, the `coefficients` and the `model.matrix`
#' are extracted from this model.
#' 
#' In order to provide a different `model.matrix`, further arguments could
#' be used. `X` is a `matrix` or a `mlogit` from which the
#' `model.matrix` is extracted. The `formula`-`data` interface
#' can also be used to construct the relevant `model.matrix`.
#' @return either a vector or a matrix.
#' @export
#' @author Yves Croissant
#' @seealso  [mlogit()] for the estimation of a multinomial logit
#' model.
#' @keywords regression
logsum <- function(coef, X = NULL, formula = NULL, data = NULL,
                   type = NULL, output = c("chid", "obs")){
    # the model.matrix is from model x
    # the coef is from model y

    # extract the coefs
    if (is.numeric(coef)) beta <- coef
    else{
        if (inherits(coef, "mlogit")) beta <- coef(coef)
        else stop("coef should be either a numeric or a mlogit object")
    }
    
    # extract the model.matrix

    # X, formula and data is NULL, in this case, extract the
    # model.matrix from the coef object
    if (is.null(X) & (is.null(data))){
        if (inherits(coef, "mlogit")){
            .idx <- idx(coef)
            X <- model.matrix(coef)
        }
        else stop("only one argument is provided, it should be a mlogit object")           
    }
    else{
        # X is provided, in this case, the index is (by priority) the
        # index attribute of X if it exists, otherwise, it is coef's
        # and conformity should be checked
        if (! is.null(X)){
            if (! inherits(X, "matrix") & ! inherits(X, "mlogit"))
                stop("X should be either a matrix or a mlogit object")
            if (is.matrix(X)){
                if (! is.null(attr(X, "index"))) idx <- attr(X, "index")
                else{
                    if (inherits(coef, "mlogit")){
                        .idx <- idx(coef)
                        if (nrow(.idx) != nrow(X)) stop("X has no index and its dimension is uncorrect")
                    }
                    else stop("no index in for the coef and the X argument")
                }
            }
            if (inherits(X, "mlogit")){
                .idx <- idx(X)
                X <- model.matrix(X)
            }
        }
        else{
            if (is.null(data)) stop("the X or data argument should be provided")
            else{
                # data is provided, if it is a mlogit object, extract
                # the model.frame
                if (inherits(data, "mlogit")) data <- model.frame(data)
                # if it is an ordinary data.frame, coerce it using the
                # index extracted from the coef argument
                if (! is.data.frame(data)) stop("data should be a data.frame")
                if ((! inherits(data, "mlogit.data")) & (! inherits(data, "dfidx"))){
                    if (is.null(attr(coef, "index")))
                        stop("no index available to compute the model.matrix")
                    else{
                        .idx <- idx(coef)
                        if (nrow(.idx) != nrow(data)) stop("uncompatible dimensions")
                        else{
                            data <- structure(data, index = .idx,
                                              class = c("mlogit.data", "data.frame"))
                        }
                    }
                }
                else .idx <- idx(data)

                if (! is.null(formula)) X <- model.matrix(formula, data)
                else{
                    if (inherits(coef, "mlogit")){
                        mf <- update(coef, data = data, estimate = FALSE)
                        .idx <- idx(mf)
                        X <- model.matrix(mf)
                    }
                    else stop("no formula provided to compute the model.matrix")
                }
            }
        }
    }
    output <- match.arg(output)
    if (! is.null(type)){
        if (! type %in% c("group", "global"))
            stop("type should be one of 'group' or 'local'")
    }
    .idx$nb <- 1:nrow(.idx)
    coefsubset <- intersect(names(beta), colnames(X))
    X <- X[, coefsubset, drop = FALSE]
    .idx$linpred <- as.numeric(crossprod(t(X[, coefsubset, drop = FALSE]), beta[coefsubset]))

    group_name <- idx_name(.idx, 2, 2)
    alt_name <- idx_name(.idx, 2, 1)
    chid_name <- idx_name(.idx, 1, 1)
    if (! is.null(group_name) & (is.null(type) || type == "group")){
        iv <- log(tapply(exp(.idx$linpred), list(idx(.idx, 1), .idx[[group_name]]), sum))
        if (output == "obs"){
            iv <- data.frame(chid = rep(rownames(iv), each = ncol(iv)),
                             group = rep(colnames(iv), nrow(iv)),
                             iv = as.numeric(t(iv)))
            names(iv)[1:2] <- c(chid_name, group_name)
            iv <- merge(.idx, iv)
            iv <- iv[order(iv$nb), "iv"]
        }
    }
    else{
        #iv <- log(with(idx, tapply(exp(linpred), chid, sum))) bug
        # (for unbalanced set of choices) fixed thanks to Malick
        # Hossain
        iv <- log(tapply(exp(.idx$linpred), idx(.idx, 1), sum, na.rm = TRUE))
        if (output == "obs"){
            iv <- data.frame(chid = rownames(iv), iv = as.numeric(iv))
            names(iv)[1] <- c(chid_name)
            iv <- merge(.idx, iv)
            iv <- iv[order(iv$nb), "iv"]
        }
    }
    iv        
}

ltm <- function(x, to = c("vec", "mat", "ltm")){
    to <- match.arg(to)
    result <- x
    if (is.null(dim(x))){
        if (to != "vec"){
            z <- length(x)
            K <- - 0.5 + sqrt(1 + 8 * z) / 2
            if (abs(K - floor(K)) > 1E-07) stop("wrong length")
            result <- matrix(0, nrow = K, ncol = K)
            result[! lower.tri(result)] <- x
            result <- t(result)
            if (to == "mat") result[upper.tri(result)] <- t(result)[upper.tri(result)]
        }
    }
    else{
        if (! identical(x, t(x))){
            # the matrix is not symetric, its upper triangular
            # elements should be 0
            if (any(x[upper.tri(x)] != 0))
                stop("the matrix is not symetric, it should have only zero ont the upper triangular part")
            if (to == "mat") result[upper.tri(x)] <- t(x)[upper.tri(x)]
            if (to == "vec") result <- t(x)[! lower.tri(x)]
        }
        else{
            result <- x
            if (to == "ltm") result[upper.tri(result)] <- 0
            if (to == "vec") result <- t(x)[! lower.tri(x)]
        }
    }
    result
}
                  
