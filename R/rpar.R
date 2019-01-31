#' random parameter objects
#' 
#' `rpar` objects contain the relevant information about estimated
#' random parameters. The homonymous function extract on `rpar` object
#' from a `mlogit` object.
#' 
#' `mlogit` objects contain an element called `rpar` which contain a
#' list of `rpar` objects, one for each estimated random
#' parameter. The `print` method prints the name of the distribution
#' and the parameter, the `summary` behave like the one for numeric
#' vectors.
#' 
#' @name rpar
#' @aliases rpar print.rpar summary.rpar
#' @param x,object a `mlogit` object,
#' @param par the name or the index of the parameters to be extracted
#'     ; if `NULL`, all the parameters are selected,
#' @param norm the coefficient used for normalization if any,
#' @param ... further arguments.
#' @param digits the number of digits
#' @param width the width of the printed output
#' @return a `rpar` object, which contains:
#'
#' - dist: the name of the distribution,
#' - mean: the first parameter of the distribution,
#' - sigma: the second parameter of the distribution,
#' - name: the name of the parameter.
#' 
#' @export
#' @author Yves Croissant
#' @seealso [mlogit()] for the estimation of a random parameters logit
#'     model.
#' @keywords regression
rpar <- function(x, par = NULL, norm = NULL, ...){
    if (is.null(par)) par <- names(x$rpar)
    if (length(par) == 1){
        result <- x$rpar[[par]]
        if (!is.null(norm)) result$norm <- abs(coef(x)[norm])
    }    
    else{
        result <- x$rpar[par]
        if (!is.null(norm))
            lapply(result, function(x){x$norm <- abs(coef(x)[norm]); x})
    }
    result
}

#' @rdname rpar
#' @export
print.rpar <- function(x, digits = max(3, getOption("digits") - 2),
                       width = getOption("width"), ...){
    dist <- switch(x$dist,
                   "n"  = "normal",
                   "ln" = "log-normal",
                   "cn" = "censored normal",
                   "t"  = "triangular",
                   "u"  = "uniform",
                   "zbu" = "uniform",
                   "zbt" = "triangular"
                   )
    npar1 <- switch(x$dist,
                    "n"  = "mean",
                    "ln" = "meanlog",
                    "cn" = "mean",
                    "t"  = "center",
                    "u"  = "center",
                    "zbu" = "center",
                    "zbt" = "center"
                    )
    
    npar2 <- switch(x$dist,
                    "n"  = "sd",
                    "ln" = "sdlog",
                    "cn" = "sd",
                    "t"  = "span",
                    "u"  = "span",
                    "zbu" = NA,
                    "zbt" = NA
                    )
    par1 <- x$mean
    par2 <- x$sigma
    cat(paste(dist, " distribution with parameters ",round(par1, 3),
              " (", npar1, ")", " and ", round(par2, 3),
              " (",npar2,")", "\n",sep = ""))
}

#' @rdname rpar
#' @export
summary.rpar <- function(object, ...){
    norm <- object$norm
    rg <- rg.rpar(object, norm)
    Q1 <- qrpar(object, norm)(0.25)
    M <-  qrpar(object, norm)(0.5)
    Q3 <- qrpar(object, norm)(0.75)
    m <- mean(object, norm)
    c('Min.' = as.numeric(rg[1]), '1st Qu.' = as.numeric(Q1), 'Median' = as.numeric(M),
      'Mean' = as.numeric(m), '3rd Qu.' = as.numeric(Q3), 'Max.' = as.numeric(rg[2]))
}

#' Plot of the distribution of estimated random parameters
#' 
#' Methods for `rpar` and `mlogit` objects which provide a plot of the
#' distribution of one or all of the estimated random parameters
#' 
#' For the `rpar` method, one plot is drawn. For the `mlogit` method,
#' one plot for each selected random parameter is drawn.
#'
#' @name plot.mlogit
#' @aliases plot.mlogit plot.rpar
#' @importFrom graphics lines plot polygon segments title
#' @param x a `mlogit` or a `rpar` object,
#' @param type the function to be plotted, whether the density or the
#'     probability density function,
#' @param par a subset of the random parameters ; if `NULL`, all the
#'     parameters are selected,
#' @param norm the coefficient's name for the `mlogit` method or the
#'     coefficient's value for the `rpar` method used for
#'     normalization,
#' @param ... further arguments, passed to `plot.rpar` for the
#'     `mlogit` method and to `plot` for the `rpar` method.
#' @export
#' @author Yves Croissant
#' @seealso [mlogit()] the estimation of random parameters logit
#'     models and [rpar()] for the description of `rpar` objects and
#'     [distribution] for functions which return informations about
#'     the distribution of random parameters.
#' @keywords regression
plot.mlogit <- function(x, par = NULL, norm = NULL, type = c("density", "probability"), ...){
    if (is.null(x$rpar)) stop("the plot method is only relevant for random parameters")
    if (is.null(par)) par <- names(x$rpar)
    if (!is.null(norm)) norm = abs(coef(x)[norm])
    rpar <- x$rpar[par]
    K <- length(rpar)
    if (K > 1){
        nrow <- 1 + (K > 2) + (K > 6)
        ncol <- 1 + (K > 1) + (K > 4)
        opar <- par(mfrow = c(nrow, ncol))
        for (i in names(rpar)){
            plot(rpar(x, i), norm = norm, type = type, ...)
        }
        par(opar)
    }
    else{
        plot(rpar(x, 1), norm = norm, type = type, ...)
    }
}

#' @rdname plot.mlogit
#' @export
plot.rpar <- function(x, norm = NULL, type = c("density", "probability"), ...){
    type <- match.arg(type)
    if (type == "density") f <- drpar
    if (type == "probability") f <- prpar
    marg <- .05
    law <- x$dist
    rg <- rg(x)
    low <- rg[1]
    np <- x$name
    neg.values <- ifelse(low < 0,TRUE,FALSE)
    up <- rg[2]
    if (!is.finite(low)) low <- qrpar(x, norm = norm)(0.005)
    if (!is.finite(up)) up <- qrpar(x, norm = norm)(0.995)
    ptstot <- gnrpoints(low, up, 1000)
    ytot <- do.call(f, list(x = x, norm = norm))(ptstot)
    ymax <- max(ytot)*(1 + marg)
    plot(ptstot, ytot, type = "n", ann = F, xaxs = "i", yaxs = "i",
         las = 1, ylim = c(0, ymax),
         xlim = c(low - marg *(up - low), up + marg * (up - low)))
    ma <- paste("Distribution of", np)
    if (neg.values){
        pourc0 <- prpar(x)(0)
        ma <- paste(ma, ":", round(pourc0 * 100, 0),"% of 0")
        if (type == "density"){
            if (low < 0){
                ptsneg <- gnrpoints(low, 0, 10)
                yneg <- do.call(f, list(x = x, norm = norm))(ptsneg)
                polygon(c(low, ptsneg, 0),c(0, yneg, 0), col = "lightblue", border = NA)
            }
        }
        else{
            segments(low - marg * (up - low), pourc0, 0, pourc0, lty = "dotted")
            segments(0, 0, 0, pourc0, lty = "dotted")
        }
    }
    lines(ptstot, ytot)
    if (law == "u" && type == "density"){
        segments(up, 0, up, drpar(x)(up))
        segments(low, 0, low, drpar(x)(low))
    }
    title(main = ma)
}

#' Correlation structure of the random parameters
#' 
#' Functions that extract the correlation structure of a mlogit object
#' 
#' @name cor.mlogit
#' @aliases cor.mlogit cov.mlogit
#' @param x an `mlogit` object with random parameters and `correlation
#'     = TRUE`.
#' @details These functions are deprecated, use [vcov][vcov.mlogit].
#'     instead.
#' @return A numerical matrix which returns either the correlation or
#'     the covariance matrix of the random parameters.
#' @export
#' @author Yves Croissant
#' @keywords regression
cor.mlogit <- function(x){
    if (is.null(x$rpar) || is.null(attr(x$rpar, 'covariance')))
        stop('cor.mlogit only relevant for random models with correlation')
    cor.mlogit <- cov.mlogit(x)
    K <- nrow(cor.mlogit)
    sd.mlogit <- sqrt(diag(cor.mlogit))
    for (i in 1:K){
        for (j in 1:K){
            cor.mlogit[i,j] <- cor.mlogit[i,j] / sd.mlogit[i] / sd.mlogit[j]
        }
    }
    cor.mlogit
}

#' @rdname cor.mlogit
#' @export
cov.mlogit <- function(x){
    if (is.null(x$rpar) || is.null(attr(x$rpar, 'covariance')))
        stop('cov.mlogit only relevant for random models with correlation')
    attr(x$rpar, 'covariance')
}

# if a normalization coefficient is used, m2norm and s2norm transform
# the two parameters of the distribution
m2norm <- function(m, dist, norm){
    switch(dist,
           "n"  = m / norm,
           "ln" = m - log(norm),
           "t"  = m / norm,
           "cn" = m / norm,
           "u"  = m / norm,
           "zbu" = m / norm,
           "zbt" = m / norm
           )
}

s2norm <- function(s, dist, norm){
    switch(dist,
           "n"  = s / norm,
           "ln" = s,
           "t"  = s / norm,
           "cn" = s / norm,
           "u"  = s / norm,
           "zbu" = s / norm,
           "zbt" = s / norm
           )
}

#' Functions used to describe the characteristics of estimated random
#' parameters
#'
#' @name distribution
#' 
#' @aliases distribution med rg stdev mean.rpar med.rpar stdev.rpar
#'     rg.mlogit mean.mlogit med.mlogit stdev.mlogit rg.rpar qrpar
#'     prpar drpar qrpar.rpar prpar.rpar drpar.rpar qrpar.mlogit
#'     prpar.mlogit drpar.mlogit
#' @param x a `mlogit` or a `rpar` object,
#' @param norm the variable used for normalization if any : for the
#'     `mlogit` method, this should be the name of the parameter, for
#'     the `rpar` method the absolute value of the parameter,
#' @param par the required parameter(s) for the `mlogit` methods
#'     (either the name or the position of the parameter(s). If
#'     `NULL`, all the random parameters are used.
#' @param y values for which the function has to be evaluated,
#' @param ... further arguments.
#'
#' @details `rpar` objects contain all the relevant information about
#'     the distribution of random parameters. These functions enables
#'     to obtain easily descriptive statistics, density, probability
#'     and quantiles of the distribution.
#' 
#' `mean`, `med`, `stdev` and `rg` compute respectively the mean, the
#' median, the standard deviation and the range of the random
#' parameter. `qrpar`, `prpar`, `drpar` return functions that compute
#' the quantiles, the probability and the density of the random
#' parameters (note that `sd` and `range` are not generic function in
#' `R` and that `median` is, but without `...`).
#' 
#' @return a numeric vector for `qrpar`, `drpar` and `prpar`, a
#'     numeric vector for `mean`, `stdev` and `med` and a numeric
#'     matrix for `rg`.
#' @author Yves Croissant
#' @seealso [mlogit()] for the estimation of random parameters logit
#'     models and [rpar()] for the description of `rpar` objects.
#' @keywords regression


# mean, med, rg and stdev methods for rpar and mlogit objects (sd and
# range are not generic, so create a stdev and a rg generic ; median
# is generic, but without ..., so create a med generic

#' @rdname distribution
#' @export
stdev <- function(x, ...){
    UseMethod("stdev")
}

#' @rdname distribution
#' @export
rg <- function(x, ...){
    UseMethod("rg")
}

#' @rdname distribution
#' @export
med <- function(x, ...){
    UseMethod("med")
}

#' @rdname distribution
#' @export
mean.rpar <- function(x, norm = NULL, ...){
    if (is.null(norm) & ! is.null(x$norm)) norm <- as.numeric(x$norm)
    dist <- x$dist
    m <- x$mean
    s <- abs(x$sigma)
    if (! is.null(norm)){
        s <- s2norm(s, dist, norm)
        m <- m2norm(m, dist, norm)
    }
    switch(dist,
           "n"  = m,
           "ln" = exp(m + 0.5 * s^2),
           "u"  = m,
           "t"  = m,
           "cn" = s * dnorm(- m / s) + m * (1 - pnorm(- m / s)),
           "zbu" = m,
           "zbt" = m
           )
}

#' @rdname distribution
#' @export
med.rpar <- function(x, norm = NULL, ...){
    if (is.null(norm) & ! is.null(x$norm)) norm <- as.numeric(x$norm)
    dist <- x$dist
    m <- x$mean
    s <- abs(x$sigma)
    if (! is.null(norm)){
        s <- s2norm(s, dist, norm)
        m <- m2norm(m, dist, norm) 
    }
    switch(dist,
           "n"  = m,
           "ln" = exp(m),
           "u"  = m,
           "t"  = m,
           "cn" = max(0, m),
           "zbu" = m,
           "zbt" = m
           )
}

#' @rdname distribution
#' @export
stdev.rpar <- function(x, norm = NULL, ...){
    if (is.null(norm) & ! is.null(x$norm)) norm <- as.numeric(x$norm)
    dist <- x$dist
    m <- x$mean
    s <- abs(x$sigma)
    if (! is.null(norm)){
        s <- s2norm(s, dist, norm)
        m <- m2norm(m, dist, norm)
    }
    switch(dist,
           "n"  = s,
           "ln" = sqrt(exp(s ^ 2) - 1) * exp(m + 0.5 * s ^ 2),
           "u"  = s ^ 2 / 3,
           "t"  = s,
           "cn" = sqrt( s ^ 2 * (1 - pnorm(- m / s)) + m * (s * dnorm(- m / s) + m * (1 - pnorm(- m / s))) -
                        (s * dnorm(- m / s) + m * (1 - pnorm(- m / s))) ^ 2),
           "zbu" = m ^ 2 / 3,
           "zbt" = m
           )
}


#' @rdname distribution
#' @export
rg.rpar <- function(x, norm = NULL, ...){
    if (is.null(norm) & ! is.null(x$norm)) norm <- as.numeric(x$norm)
    dist <- x$dist
    m <- x$mean
    s <- abs(x$sigma)
    if (! is.null(norm)){
        s <- s2norm(s, dist, norm)
        m <- m2norm(m, dist, norm)
    }
    result <- switch(dist,
                     "n"  = c(-Inf, +Inf),
                     "ln" = c(0, +Inf),
                     "u"  = c(m - s, m + s),
                     "t"  = c(m - s, m + s),
                     "cn" = c(0, +Inf),
                     "zbu" = c(0, 2 * m),
                     "zbt" = c(0, 2 * m)
                     )
    names(result) <- c('Min.', 'Max.')
    result
}

#' @rdname distribution
#' @export
mean.mlogit <- function(x, par = NULL, norm = NULL, ...){
    if (!is.null(norm)) norm <- abs(coef(x)[norm])
    if (is.null(par)) par <- names(x$rpar)
    sapply(x$rpar[par], function(x) mean(x, norm = norm, ...))
}

#' @rdname distribution
#' @export
med.mlogit <- function(x, par = NULL, norm = NULL, ...){
    if (!is.null(norm)) norm <- abs(coef(x)[norm])
    if (is.null(par)) par <- names(x$rpar)
    sapply(x$rpar[par], function(x) med(x, norm = norm, ...))
}

#' @rdname distribution
#' @export
stdev.mlogit <- function(x, par = NULL, norm = NULL, ...){
    if (!is.null(norm)) norm <- abs(coef(x)[norm])
    if (is.null(par)) par <- names(x$rpar)
    sapply(x$rpar[par], function(x) stdev(x, norm = norm, ...))
}

#' @rdname distribution
#' @export
rg.mlogit <- function(x, par = NULL, norm = NULL, ...){
    if (!is.null(norm)) norm <- abs(coef(x)[norm])
    if (is.null(par)) par <- names(x$rpar)
    if (length(par) > 1)
        sapply(x$rpar[par], function(x) rg(x, norm = norm, ...))
    else rg(x$rpar[[par]])
}

#' @rdname distribution
#' @export
qrpar <- function(x, ...){
    UseMethod("qrpar")
}

#' @rdname distribution
#' @export
prpar <- function(x, ...){
    UseMethod("prpar")
}

#' @rdname distribution
#' @export
drpar <- function(x, ...){
    UseMethod("drpar")
}

#' @rdname distribution
#' @export
qrpar.rpar <- function(x, norm = NULL, ...){
    dist <- x$dist
    m <- x$mean
    s <- abs(x$sigma)
    if (!is.null(norm)){
        s <- s2norm(s, dist, norm)
        m <- m2norm(m, dist, norm)
    }
    switch(dist,
           "n"  = function(x = (1:9) / 10) qnorm(x, m, s),
           "ln" = function(x = (1:9) / 10) qlnorm(x, m, s),
           "u"  = function(x = (1:9) / 10) qunif(x, m - s, m + s),
           "t"  = function(x = (1:9) / 10) (m - s + sqrt(2 * s^2 * x)) *
                                       (x <= 0.5) + (m + s - sqrt(2 * s^2 *(1 - x))) *(x > 0.5),
           "cn" = function(x=(1:9)/10) max(0, qnorm(x, m, s)),
           "zbu" = function(x = (1:9) / 10) qunif(x, 0, 2 * m),
           "zbt" = function(x = (1:9) / 10) sqrt(2 * m ^ 2 * x) *
                                            (x <= 0.5) + (2 * m - sqrt(2 * m ^ 2 *(1 - x))) * (x > 0.5)
           )
}

#' @rdname distribution
#' @export
prpar.rpar <- function(x, norm = NULL, ...){
    dist <- x$dist
    m <- x$mean
    s <- abs(x$sigma)
    if (!is.null(norm)){
        s <- s2norm(s, dist, norm)
        m <- m2norm(m, dist, norm)
    }
    switch(dist,
           "n"  = function(x) pnorm(x, m, s),
           "ln" = function(x) plnorm(x, m, s),
           "u"  = function(x) punif(x, m - s, m + s),
           "t"  = function(x) (x >= (m - s) & x < m) * (x - m + s) ^ 2 / (2 * s ^ 2) +
                              (x >= m & x <= (m + s)) * (1 - (m + s - x) ^ 2 / (2 * s ^ 2)) + (x > (m + s)) * 1 + 0,
           "cn" = function(x) pnorm(x, m, s),
           "zbu" = function(x) punif(x, 0, 2 * m),
           "zbt"  = function(x) (x >= 0 & x < m) * x ^ 2 / (2 * m ^ 2) +
                                (x >= m & x <= 2 * m) * (1 - (2 * m) ^ 2 / (2 * m ^ 2)) + (x > 2 * m) * 1 + 0

           )
}

#' @rdname distribution
#' @export
drpar.rpar <- function(x, norm = NULL, ...){
    dist <- x$dist
    m <- x$mean
    s <- abs(x$sigma)
    if (!is.null(norm)){
        s <- s2norm(s, dist, norm)
        m <- m2norm(m, dist, norm)
    }
    switch(dist,
           "n"  = function(x) dnorm(x, m, s),
           "ln" = function(x) dlnorm(x, m, s),
           "u"  = function(x) (1 / s + x * 0) * (x >= m - s & x <= m + s) + 0,
           "t"  = function(x) (x >= (m - s) & x < m) * (x - m + s) / s ^ 2 +
                              (x >= m & x <= (m + s)) * (s + m - x) / s ^ 2 + 0,
           "cn" = function(x) dnorm(x, m, s),
           "zbu" = function(x) dunif(x, 0, 2 * m),
           "zbt" = function(x) (x >= 0 & x < m) * x / m ^ 2 +
                               (x >= m & x <= 2 * m) * (2 * m - x) / m ^ 2 + 0,
           )
}

#' @rdname distribution
#' @export
qrpar.mlogit <- function(x, par = 1, y = NULL, norm = NULL, ...){
    if (is.null(rpar))
        stop("qrpar function only relevant for random parameters models")
    if (length(par) > 1)
        stop("only one parameter should be selected")
    x <- x$rpar[[par]]
    if (is.null(norm)) norm <- abs(coef(x)[norm])
    if (is.null(y)){
        qrpar(x, norm = norm, ...)
    }
    else{
        qrpar(x, norm = norm, ...)(y)
    }
}

#' @rdname distribution
#' @export
prpar.mlogit <- function(x, par = 1, y = NULL, norm = NULL, ...){
    if (is.null(rpar))
        stop("prpar function only relevant for random parameters models")
    if (length(par) > 1)
        stop("only one parameter should be selected")
    x <- x$rpar[[par]]
    if (!is.null(norm)) norm <- coef(x)[norm]
    if (is.null(y)){
        prpar(x, norm = norm, ...)
    }
    else{
        prpar(x, norm = norm, ...)(y)
    }
}

#' @rdname distribution
#' @export
drpar.mlogit <- function(x, par = 1, y = NULL, norm = NULL, ...){
    if (is.null(rpar))
        stop("drpar function only relevant for random parameters models")
    if (length(par) > 1)
        stop("only one parameter should be selected")
    x <- x$rpar[[par]]
    if (!is.null(norm)) norm <- coef(x)[norm]
    if (is.null(y)){
        drpar(x, norm = norm, ...)
    }
    else{
        drpar(x, norm = norm, ...)(y)
    }
}

