#' @rdname miscmethods.mlogit
#' @export
predict.mlogit <- function(object, newdata = NULL, returnData = FALSE, ...){
    # if no newdata is provided, use the mean of the model.frame
    if (is.null(newdata)){
        newdata <- mean(object$model)
        #YC2020/03/05 coerce it to a data.frame so that it gets transformed 
        newdata <- as.data.frame(newdata)
    }
    # if newdata is not a mlogit.data, it is coerced below
    if ((! inherits(newdata, "mlogit.data")) & (! inherits(newdata, "dfidx"))){
        if (FALSE){
            rownames(newdata) <- NULL
            lev <- colnames(object$probabilities)
            J <- length(lev)
            choice.name <- attr(model.frame(object), "choice")
            if (nrow(newdata) %% J)
                stop("the number of rows of the data.frame should be a multiple of the number of alternatives")
            attr(newdata, "index") <- data.frame(chid = rep(1:(nrow(newdata) %/% J ), each = J), alt = rep(lev, J))
            attr(newdata, "class") <- c("mlogit.data", "data.frame")
            if (is.null(newdata[['choice.name']])){
                newdata[[choice.name]] <- FALSE
                newdata[[choice.name]][1] <- TRUE # probit and hev requires that one (arbitrary) choice is TRUE
            }
        }
        else{
            # New stuff, coerce manually to an dfidx object
            lev <- colnames(object$probabilities)
            J <- length(lev)
            choice.name <- paste(deparse(formula(object)[[2]]))
            if (nrow(newdata) %% J)
                stop("the number of rows of the data.frame should be a multiple of the number of alternatives")
            nbid <- nrow(newdata) %/% J
            newdata$idx <- data.frame(chid = factor(rep(1:nbid, each = J)), alt = factor(rep(lev, nbid), levels = lev))
            class(newdata$idx) <- c("idx", "data.frame")
            attr(newdata$idx, "ids") <- c(1, 2)
            attr(newdata, "clseries") <- c("xseries_mlogit", "xseries")
            attr(newdata, "class") <- c("dfidx_mlogit", "dfidx", "data.frame")
            #YC2020/03/05 this seems incorect, replace 'choice.name' by choice.name
   #        if (is.null(newdata[['choice.name']])){
            if (is.null(newdata[[choice.name]])){
                newdata[[choice.name]] <- FALSE
                newdata[[choice.name]][1] <- TRUE # probit and hev requires that one (arbitrary) choice is TRUE
            }
        }
    }
    # if the updated model requires the use of mlogit.data, suppress all
    # the relevant arguments
    m <- match(c("choice", "shape", "varying", "sep",
                 "alt.var", "chid.var", "alt.levels",
                 "opposite", "drop.index", "id", "ranked"),
               names(object$call), 0L)
    if (sum(m) > 0) object$call <- object$call[ - m]
    # update the model and get the probabilities
    newobject <- update(object, start = coef(object, fixed = TRUE), data = newdata, iterlim = 0, print.level = 0)
#    newobject <- update(object, start = coef(object), data = newdata, iterlim = 0, print.level = 0)
    result <- newobject$probabilities
    if (nrow(result) == 1){
        result <- as.numeric(result)
        names(result) <- colnames(object$probabilities)
    }
    if (returnData) attr(result, "data") <- newdata
    result
}

# !!! on vire model.response.mlogit de @aliases

#' Methods for mlogit objects
#'
#' Miscellaneous methods for `mlogit` objects.
#' 
#' @name miscmethods.mlogit
#' @aliases residuals.mlogit df.residual.mlogit terms.mlogit
#'     model.matrix.mlogit update.mlogit
#'     print.mlogit  summary.mlogit print.summary.mlogit
#'     predict.mlogit fitted.mlogit coef.mlogit
#'     coef.summary.mlogit
#' @param x,object an object of class `mlogit`
#' @param digits the number of digits,
#' @param width the width of the printing,
#' @param new an updated formula for the `update` method,
#' @param newdata a `data.frame` for the `predict` method,
#' @param outcome a boolean which indicates, for the `fitted` and the
#'     `residuals` methods whether a matrix (for each choice, one
#'     value for each alternative) or a vector (for each choice, only
#'     a value for the alternative chosen) should be returned,
#' @param type one of `outcome` (probability of the chosen
#'     alternative), `probabilities` (probabilities for all the
#'     alternatives), `parameters` for individual-level random
#'     parameters for the fitted method, how the correlated random
#'     parameters should be displayed : `"chol"` for the estimated
#'     parameters (the elements of the Cholesky decomposition matrix),
#'     `"cov"` for the covariance matrix and `"cor"` for the
#'     correlation matrix and the standard deviations,
#' @param returnData for the `predict` method, if `TRUE`, the data is
#'     returned as an attribute,
#' @param n,m see [dfidx::idx()]
#' @param ... further arguments.
#' @rdname miscmethods.mlogit
NULL

#' Marginal effects of the covariates
#' 
#' The `effects` method for `mlogit` objects computes the marginal
#' effects of the selected covariate on the probabilities of choosing the
#' alternatives
#' 
#' @name effects.mlogit
#' @param object a `mlogit` object,
#' @param covariate the name of the covariate for which the effect should be
#' computed,
#' @param type the effect is a ratio of two marginal variations of the
#' probability and of the covariate ; these variations can be absolute
#' `"a"` or relative `"r"`. This argument is a string that contains
#' two letters, the first refers to the probability, the second to the
#' covariate,
#' @param data a data.frame containing the values for which the effects should
#' be calculated. The number of lines of this data.frame should be equal to the
#' number of alternatives,
#' @param ... further arguments.
#' @return If the covariate is alternative specific, a \eqn{J \times J} matrix is
#' returned, \eqn{J} being the number of alternatives. Each line contains the
#' marginal effects of the covariate of one alternative on the probability to
#' choose any alternative. If the covariate is individual specific, a vector of
#' length \eqn{J} is returned.
#' @export
#' @author Yves Croissant
#' @seealso  [mlogit()] for the estimation of multinomial logit
#' models.
#' @keywords regression
#' @examples
#' 
#' data("Fishing", package = "mlogit")
#' Fish <- dfidx(Fishing, varying = 2:9, choice = "mode")
#' m <- mlogit(mode ~ price | income | catch, data = Fish)
#' # compute a data.frame containing the mean value of the covariates in
#' # the sample
#' z <- with(Fish, data.frame(price = tapply(price, idx(m, 2), mean),
#'                            catch = tapply(catch, idx(m, 2), mean),
#'                            income = mean(income)))
#' # compute the marginal effects (the second one is an elasticity
#' ## IGNORE_RDIFF_BEGIN
#' effects(m, covariate = "income", data = z)
#' ## IGNORE_RDIFF_END
#' effects(m, covariate = "price", type = "rr", data = z)
#' effects(m, covariate = "catch", type = "ar", data = z)
effects.mlogit <- function(object, covariate = NULL,
                           type = c("aa", "ar", "rr", "ra"),
                           data = NULL, ...){
    type <- match.arg(type)
    if (is.null(data)){
        P <- predict(object, returnData = TRUE)
        data <- attr(P, "data")
        attr(P, "data") <- NULL
    }
    else P <- predict(object, data)
    newdata <- data
    J <- length(P)
    alt.levels <- names(P)
    pVar <- substr(type, 1, 1)
    xVar <- substr(type, 2, 2)
#    cov.list <- lapply(attr(formula(object), "rhs"), as.character)
    nrhs <- length(formula(object))[2]
    cov.list <- vector(length = 3, mode = "list")
    for (i in 1:nrhs) cov.list[[i]] <-
                          attr(terms(formula(object), rhs = i), "term.labels")
    rhs <- sapply(cov.list, function(x) length(na.omit(match(x, covariate))) > 0)
    rhs <- (1:length(cov.list))[rhs]
    eps <- 1E-5
    if (rhs %in% c(1, 3)){
        if (rhs == 3){
            theCoef <- paste(alt.levels, covariate, sep = ":")
            theCoef <- coef(object)[theCoef]
        }
        else theCoef <- coef(object)[covariate]
        me <- c()
        for (l in 1:J){
            newdata[l, covariate] <- data[l, covariate] + eps
            newP <- predict(object, newdata)
            me <- rbind(me, (newP - P) / eps)
            newdata <- data
        }
        if (pVar == "r") me <- t(t(me) / P)
        if (xVar == "r") me <- me * matrix(rep(data[[covariate]], J), J)
        dimnames(me) <- list(alt.levels, alt.levels)
    }
    if (rhs == 2){
        newdata[, covariate] <- data[, covariate] + eps
        newP <- predict(object, newdata)
        me <- (newP - P) / eps
        if (pVar == "r") me <- me / P
        if (xVar == "r") me <- me * data[[covariate]]
        names(me) <- alt.levels
    }
    if (inherits(me, "xseries")) attr(me, "idx") <- NULL
    unclass(me)
}

get_covariates <- function(x, rhs = 1){
    .form <- formula(x)
    if (rhs == 2 & length(.form)[2] == 1) result <- NULL
    else {
        z <- formula(.form, lhs = 0, rhs = rhs)
        result <- strsplit(as.character(z)[[2]], "\\s*\\+\\s*")[[1]]
        if (length(result) == 1){
            if (result == "0") result <- NULL
        }
    }
    result
}

pred_coefs <- function(coefs, model, data = NULL){
    if (is.null(data)){
        preds <- update(model, param = coefs)
    } else {
        if (inherits(data, "dfidx_mlogit")){
            class(data) <- setdiff(class(data), "dfidx_mlogit")
        }
        preds <- update(model, param = coefs, data = data)

    }
    as.numeric(t(preds))
    ## preds %>% 
    ##     as_tibble %>%
    ##     tidyr::pivot_longer(1:4) %>%
    ##     pull(value)
}

preds_se <- function(model){
    jacob <- numDeriv::jacobian(pred_coefs, coef(model), model = model)
    sqrt(apply(jacob, 1, function(x) x %*% vcov(model) %*% x))
}

slope_coefs <- function(coefs, model, covar, alt = NULL, eps = NULL){
    dta <- model.frame(model)
    probs_init <- pred_coefs(coefs, model, dta)
    id <- idx(model, 2, 1)
    if (is.null(eps))  eps <- sd(dta[[covar]], na.rm = TRUE) * sqrt(.Machine$double.eps)
    dta2 <- dta
    if (! is.null(alt)){
        dta2[[covar]] <- dta2[[covar]] + eps * (id == alt)
    } else {
        dta2[[covar]] <- dta2[[covar]] + eps
    }
    probs_fin <- pred_coefs(coefs, model, dta2)
    (probs_fin - probs_init) / eps
}


sd_slope <- function(model, covar, alt = NULL){
    jacob <- numDeriv::jacobian(slope_coefs, coef(model),
                                model = model, covar = covar, alt = alt)
    sqrt(apply(jacob, 1, function(x) x %*% vcov(model) %*% x))
}

#' Compute the predictions and their standard deviations
#' 
#' `preds` is an alternative to the `predict` method
#' 
#' @name preds
#' @param object a fitted model
#' @param ... further arguments for the summary method
#' @details
#' the standard errors are computed using the delta method.
#' @return a `preds` object, which inherits the `dfidx` class and
#'     therefore contains an index column.
#' @importFrom stats sd
#' @export
#' @author Yves Croissant
#' @seealso [mlogit()] for the estimation of a multinomial logit
#'     model.
#' @keywords regression
#' @export
preds <- function(object){
    preds <- pred_coefs(coef(object), object)
    se <- preds_se(object)
    mf <- model.frame(object)
    id1 <- idx(mf, 1)
    id2 <- idx(mf, 2)
    result <- data.frame(id1, id2, preds, se)
    names(result) <- c(idx_name(mf, 1), idx_name(mf, 2), "estimate", "std.error")
    result <- dfidx(result)
    if (inherits(mf, "tbl_df")){
        class(result) <- c("preds", "tbl_dfidx", "dfidx", "tbl_df", "tbl", "data.frame")
    } else {
        class(result) <- c("preds", "dfidx", "data.frame")
    }
    result
}

#' Compute the slopes and their standard deviations
#' 
#' `preds` is an alternative to the `predict` method
#' 
#' @name slps
#' @param object a fitted model
#' @param ... further arguments for the summary method
#' @details the standard errors are computed using the delta method.
#' @return a `slps` object, which inherits the `tbl_df` `tbl` and
#'     `data.frame` classes
#' @export
#' @author Yves Croissant
#' @seealso [mlogit()] for the estimation of a multinomial logit
#'     model.
#' @keywords regression
#' @export
slps <- function(object){
    is.tibble <- inherits(model.frame(object), "tbl_df")
    slps <- data.frame()
    alt <- idx(object, 2, 1)
    levs <- levels(alt)
    id <- unique(idx(object, 1, 1))
    ids <- rep(id, each = length(levs))
    alts <- rep(levs, length(id))
    covariates_1 <- get_covariates(object)
    covariates_2 <- get_covariates(object, 2)
    if (! is.null(covariates_2)){
        for (i in covariates_2){
            slp_est <- slope_coefs(coef(object), object, covar = i)
            slp_sd <- sd_slope(object, covar = i)
            .statistic <- slp_est / slp_sd
            .pvalue <- 2 * pnorm(abs(.statistic), lower.tail = FALSE)
            aslp <- data.frame(term = i,
                               id = ids,
                               alt = alts,
                               estimate = slp_est,
                               std.error = slp_sd,
                               statistic = .statistic,
                               p.value = .pvalue)
            slps <- rbind(slps, aslp)
        }
    }
    if (! is.null(covariates_1)){
        for (i in covariates_1){
            for (l in levs){
                slp_est <- slope_coefs(coef(object), object, covar = i, alt = l)
                slp_sd <- sd_slope(object, covar = i, alt = l)
                .statistic <- slp_est / slp_sd
                .pvalue <- 2 * pnorm(abs(.statistic), lower.tail = FALSE)
                aslp <- data.frame(term = paste(i, l, sep = "_"),
                                   id = ids,
                                   alt = alts,
                                   estimate = slp_est,
                                   std.error = slp_sd,
                                   statistic = .statistic,
                                   p.value = 2 * pnorm(abs(.pvalue), lower.tail = FALSE))
                slps <- rbind(slps, aslp)
            }
        }
    }
    if (is.tibble){
        class(slps) <- c("slps", "tbl_df", "tbl", "data.frame")
    } else {
        class(slps) <- c("slps", "data.frame")
    }        
    slps
}

#' @rdname slps
#' @method summary slps
#' @importFrom stats aggregate
#' @export
summary.slps <- function(object, ...){
    ## dplyr::summarise(object,
    ##           estimate = mean(estimate),
    ##           std.error = sqrt(sum(std.error ^ 2) / nrow(object)),
    ##           statistic = estimate / std.error,
    ##           p.value = 2 * pnorm(abs(statistic), lower.tail = FALSE),
    ##           .by = c(term, alt)) %>% print()
    .est <- aggregate(estimate ~ term + alt, object, mean)
    .sd <- aggregate(std.error ~ term + alt, object, function(x) sqrt(sum(x ^ 2) / nrow(object)))
    result <- merge(.est, .sd, by = c("term", "alt"))
    result <- result[order(result$term, result$alt), ]
    is.tibble <- inherits(object, "tbl_df")
    result$statistic <- result$estimate / result$std.error
    result$p.value <- 2 * pnorm(abs(result$statistic), lower.tail = FALSE)
    if (is.tibble){
        class(result) <- c("tbl_df", "tbl", "data.frame")
    } else {
        class(result) <- "data.frame"
    }
    result
}

#' @rdname preds
#' @method summary preds
#' @export
summary.preds <- function(object, ...){
    is.tibble <- inherits(object, "tbl_df")
    alt <- idx(object, 2, 1)
    cls <- class(object)
    class(object) <- setdiff(cls, c("dfidx", "preds"))
    object$alt <- alt
    e <- tapply(object$estimate, object$alt, mean)
    s <- tapply(object$std.error, object$alt, function(x) sqrt(sum(x ^ 2) / nrow(object)))
    ## result <- data.frame(
    ## object %>% add_column(alt = alt) %>% 
    ## summarise(estimate = mean(estimate),
    ##               std.error = sqrt(sum(std.error ^ 2) / nrow(object)),
    ##               .by = alt)
    result <- data.frame(alt = names(e), estimate = unname(e), std.error = unname(s))
    if (is.tibble) class(result) <- c("tbl_df", "tbl", "data.frame")
    result
}
