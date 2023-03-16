#' Multinomial logit model
#' 
#' Estimation by maximum likelihood of the multinomial logit model,
#' with alternative-specific and/or individual specific variables.
#' 
#' @name mlogit
#' @aliases mlogit
#' @import Formula
#' @import dfidx
#' @importFrom stats as.formula coef dlnorm dnorm formula logLik
#' @importFrom stats model.frame model.matrix model.response
#' @importFrom stats na.omit pchisq plnorm pnorm predict
#' @importFrom stats printCoefmat punif qlnorm qnorm qunif
#' @importFrom stats relevel reshape residuals rnorm runif
#' @importFrom stats terms update vcov dunif effects
#' @importFrom statmod gauss.quad
#' @importFrom MASS ginv
#' @importFrom zoo index
#' @param formula a symbolic description of the model to be estimated,
#' @param data the data: an `mlogit.data` object or an ordinary
#'     `data.frame`,
#' @param subset an optional vector specifying a subset of
#'     observations for `mlogit`,
#' @param weights an optional vector of weights,
#' @param na.action a function which indicates what should happen when
#'     the data contains `NA`s,
#' @param start a vector of starting values,
#' @param alt.subset a vector of character strings containing the
#'     subset of alternative on which the model should be estimated,
#' @param reflevel the base alternative (the one for which the
#'     coefficients of individual-specific variables are normalized to
#'     0),
#' @param nests a named list of characters vectors, each names being a
#'     nest, the corresponding vector being the set of alternatives
#'     that belong to this nest,
#' @param un.nest.el a boolean, if `TRUE`, the hypothesis of unique
#'     elasticity is imposed for nested logit models,
#' @param unscaled a boolean, if `TRUE`, the unscaled version of the
#'     nested logit model is estimated,
#' @param heterosc a boolean, if `TRUE`, the heteroscedastic logit
#'     model is estimated,
#' @param rpar a named vector whose names are the random parameters
#'     and values the distribution : `'n'` for normal, `'l'` for
#'     log-normal, `'t'` for truncated normal, `'u' ` for uniform,
#' @param probit if `TRUE`, a multinomial porbit model is estimated,
#' @param R the number of function evaluation for the gaussian
#'     quadrature method used if `heterosc = TRUE`, the number of
#'     draws of pseudo-random numbers if `rpar` is not `NULL`,
#' @param correlation only relevant if `rpar` is not `NULL`, if true,
#'     the correlation between random parameters is taken into
#'     account,
#' @param halton only relevant if `rpar` is not `NULL`, if not `NULL`,
#'     halton sequence is used instead of pseudo-random numbers. If
#'     `halton = NA`, some default values are used for the prime of
#'     the sequence (actually, the primes are used in order) and for
#'     the number of elements droped. Otherwise, `halton` should be a
#'     list with elements `prime` (the primes used) and `drop` (the
#'     number of elements droped).
#' @param random.nb only relevant if `rpar` is not `NULL`, a
#'     user-supplied matrix of random,
#' @param panel only relevant if `rpar` is not `NULL` and if the data
#'     are repeated observations of the same unit ; if `TRUE`, the
#'     mixed-logit model is estimated using panel techniques,
#' @param estimate a boolean indicating whether the model should be
#'     estimated or not: if not, the `model.frame` is returned,
#' @param seed the seed to use for random numbers (for mixed logit and
#'     probit models),
#' @param ... further arguments passed to `mlogit.data` or
#'     `mlogit.optim`.
#' 
#' @details For how to use the formula argument, see [Formula()].
#' 
#' The `data` argument may be an ordinary `data.frame`. In this case,
#' some supplementary arguments should be provided and are passed to
#' [mlogit.data()]. Note that it is not necessary to indicate the
#' choice argument as it is deduced from the formula.
#' 
#' The model is estimated using the [mlogit.optim()].
#' function.
#' 
#' The basic multinomial logit model and three important extentions of
#' this model may be estimated.
#' 
#' If `heterosc=TRUE`, the heteroscedastic logit model is estimated.
#' `J - 1` extra coefficients are estimated that represent the scale
#' parameter for `J - 1` alternatives, the scale parameter for the
#' reference alternative being normalized to 1. The probabilities
#' don't have a closed form, they are estimated using a gaussian
#' quadrature method.
#' 
#' If `nests` is not `NULL`, the nested logit model is estimated.
#' 
#' If `rpar` is not `NULL`, the random parameter model is estimated.
#' The probabilities are approximated using simulations with `R` draws
#' and halton sequences are used if `halton` is not
#' `NULL`. Pseudo-random numbers are drawns from a standard normal and
#' the relevant transformations are performed to obtain numbers drawns
#' from a normal, log-normal, censored-normal or uniform
#' distribution. If `correlation = TRUE`, the correlation between the
#' random parameters are taken into account by estimating the
#' components of the cholesky decomposition of the covariance
#' matrix. With G random parameters, without correlation G standard
#' deviations are estimated, with correlation G * (G + 1) /2
#' coefficients are estimated.

#' @return An object of class `"mlogit"`, a list with elements:
#' 
#' - coefficients: the named vector of coefficients,
#' - logLik: the value of the log-likelihood,
#' - hessian: the hessian of the log-likelihood at convergence,
#' - gradient: the gradient of the log-likelihood at convergence,
#' - call: the matched call,
#' - est.stat: some information about the estimation (time used,
#' optimisation method),
#' - freq: the frequency of choice,
#' - residuals: the residuals,
#' - fitted.values: the fitted values,
#' - formula: the formula (a `Formula` object),
#' - expanded.formula: the formula (a `formula` object),
#' - model: the model frame used,
#' - index: the index of the choice and of the alternatives.
#' 
#' @export
#' @author Yves Croissant
#' @seealso [mlogit.data()] to shape the data. [nnet::multinom()] from
#'     package `nnet` performs the estimation of the multinomial logit
#'     model with individual specific variables. [mlogit.optim()]
#'     details about the optimization function.
#' @references
#'
#' \insertRef{MCFA:73}{mlogit}
#' 
#' \insertRef{MCFA:74}{mlogit}
#' 
#' \insertRef{TRAI:09}{mlogit}
#' 
#' @keywords regression
#' @examples
#' ## Cameron and Trivedi's Microeconometrics p.493 There are two
#' ## alternative specific variables : price and catch one individual
#' ## specific variable (income) and four fishing mode : beach, pier, boat,
#' ## charter
#' 
#' data("Fishing", package = "mlogit")
#' Fish <- dfidx(Fishing, varying = 2:9, shape = "wide", choice = "mode")
#' 
#' ## a pure "conditional" model
#' summary(mlogit(mode ~ price + catch, data = Fish))
#' 
#' ## a pure "multinomial model"
#' summary(mlogit(mode ~ 0 | income, data = Fish))
#' 
#' ## which can also be estimated using multinom (package nnet)
#' summary(nnet::multinom(mode ~ income, data = Fishing))
#' 
#' ## a "mixed" model
#' m <- mlogit(mode ~ price + catch | income, data = Fish)
#' summary(m)
#' 
#' ## same model with charter as the reference level
#' m <- mlogit(mode ~ price + catch | income, data = Fish, reflevel = "charter")
#' 
#' ## same model with a subset of alternatives : charter, pier, beach
#' m <- mlogit(mode ~ price + catch | income, data = Fish,
#'             alt.subset = c("charter", "pier", "beach"))
#' 
#' ## model on unbalanced data i.e. for some observations, some
#' ## alternatives are missing
#' # a data.frame in wide format with two missing prices
#' Fishing2 <- Fishing
#' Fishing2[1, "price.pier"] <- Fishing2[3, "price.beach"] <- NA
#' mlogit(mode ~ price + catch | income, Fishing2, shape = "wide", varying = 2:9)
#' 
#' # a data.frame in long format with three missing lines
#' data("TravelMode", package = "AER")
#' Tr2 <- TravelMode[-c(2, 7, 9),]
#' mlogit(choice ~ wait + gcost | income + size, Tr2)
#' 
#' ## An heteroscedastic logit model
#' data("TravelMode", package = "AER")
#' hl <- mlogit(choice ~ wait + travel + vcost, TravelMode, heterosc = TRUE)
#' 
#' ## A nested logit model
#' TravelMode$avincome <- with(TravelMode, income * (mode == "air"))
#' TravelMode$time <- with(TravelMode, travel + wait)/60
#' TravelMode$timeair <- with(TravelMode, time * I(mode == "air"))
#' TravelMode$income <- with(TravelMode, income / 10)
#' # Hensher and Greene (2002), table 1 p.8-9 model 5
#' TravelMode$incomeother <- with(TravelMode, ifelse(mode %in% c('air', 'car'), income, 0))
#' nl <- mlogit(choice ~ gcost + wait + incomeother, TravelMode,
#'              nests = list(public = c('train', 'bus'), other = c('car','air')))
#'              
#' # same with a comon nest elasticity (model 1)
#' nl2 <- update(nl, un.nest.el = TRUE)
#' 
#' ## a probit model
#' \dontrun{
#' pr <- mlogit(choice ~ wait + travel + vcost, TravelMode, probit = TRUE)
#' }
#' 
#' ## a mixed logit model
#' \dontrun{
#' rpl <- mlogit(mode ~ price + catch | income, Fishing, varying = 2:9,
#'               rpar = c(price= 'n', catch = 'n'), correlation = TRUE,
#'               alton = NA, R = 50)
#' summary(rpl)
#' rpar(rpl)
#' cor.mlogit(rpl)
#' cov.mlogit(rpl)
#' rpar(rpl, "catch")
#' summary(rpar(rpl, "catch"))
#' }
#' 
#' # a ranked ordered model
#' data("Game", package = "mlogit")
#' g <- mlogit(ch ~ own | hours, Game, varying = 1:12, ranked = TRUE,
#'             reflevel = "PC", idnames = c("chid", "alt"))
mlogit <- function(formula, data, subset, weights, na.action, start = NULL,
                   alt.subset = NULL, reflevel = NULL,
                   nests = NULL, un.nest.el = FALSE, unscaled = FALSE,
                   heterosc = FALSE, rpar = NULL, probit = FALSE,
                   R = 40, correlation = FALSE, halton = NULL, random.nb = NULL,
                   panel = FALSE, estimate = TRUE, seed = 10, ...){
    callT <- match.call(expand.dots = TRUE)
    formula <- callT$formula <- Formula(formula)
    nframe <- length(sys.calls())
    start.time <- proc.time()
    mldata <- callT

    # 1 ######################################################
    # Check what kind of model is estimated
    ##########################################################
    mt <- mlogit.model(formula, nests, heterosc, rpar, probit)
    multinom.logit <- mt['multinom']
    nested.logit <- mt['nested']
    pair.comb.logit <- mt['rpl']
    probit <- mt['probit']
    wlogit <- mt['wlogit']
    heterosc.logit <- mt['heterosc']
    mixed.logit <- mt["mixed"]
    if (multinom.logit) callT$method <- 'nr'

    # 1 ######################################################
    # Subset the data.frame if necessary
    ##########################################################

    # Subset is an argument of mlogit and dfidx. In any case, it is
    # simpler to subset the df from mlogit, the df being either an
    # ordinary df or a dfidx ; if the data argument is an ordinary
    # data.frame, call dfidx without the subset argument


    if ("subset" %in% names(mldata)){
        sub_call <- mldata
        m <- match(c("data", "subset"), names(sub_call), 0)
        sub_call <- sub_call[c(1, m)]
        names(sub_call)[2] <- "x"
        sub_call[[1]] <- as.name("subset")
        data <- eval(sub_call, parent.frame())
    }
    # 2 ################################################################
    # Run dfidx (or mlogit.data for backward compatibility) if necessary
    ####################################################################
    # check if any of the mlogit.data argument is present
    
    common_args <- c("data", "drop.index", "shape",
                     "choice", "varying", "sep", "opposite", "ranked")
    dfidx_args <- c("idx", "as.factor", "pkg", "fancy.row.names", "idnames", "levels")
    mlogit.data_args <- c("alt.var", "chid.var", "alt.levels", "id.var", "group.var")
    match_dfidx <- match(dfidx_args, names(mldata), 0L)
    match_mlogit.data <- match(mlogit.data_args, names(mldata), 0L)
    match_common <- match(common_args, names(mldata), 0L)
    to_dfidx <- any(as.logical(match_dfidx))
    to_mlogit.data <- any(as.logical(match_mlogit.data))
    # remove the first element which matches with data
    to_common <- any(as.logical(match_common[-c(1:2)]))
    # is_data.frame is TRUE if the data frame is not a dfidx object
    is_data.frame <- ! inherits(data, "dfidx") & ! inherits(data, "mlogit.data")

    # if the data frame is already a dfidx and some dfidx arguments
    # are set, stop
    if ((to_dfidx | to_mlogit.data | to_common) & ! is_data.frame)
        stop("data is already a dfidx object")

    # if dfidx and mlogit arguments are mixed, the stop
    if (to_dfidx & to_mlogit.data)
        stop("some specificic arguments of mlogit.data and dfidx are introduced")

    # it data is an ordinary data.frame, run dfidx or mlogit.data

    if (is_data.frame){
        if (to_mlogit.data){
#            warning("the use of mlogit.data is deprecated, use dfidx' arguments instead")
            mldata[[1L]] <- as.name("mlogit.data")
            m <- match(c(common_args, mlogit.data_args),
                       names(mldata), 0L)
        }
        else{
            mldata$pkg <- "mlogit"
            mldata[[1L]] <- as.name("dfidx")
            m <- match(c(common_args, dfidx_args),
                       names(mldata), 0L)
            to_dfidx <- TRUE
        }
        if (to_mlogit.data | to_dfidx){
            # The choice argument is irrelevant, use the response
            # indicated in the formula
            mldata <- mldata[c(1L, m)]
            mldata$choice <- paste(deparse(attr(formula, "lhs")[[1]]))
            # This works only if dfidx is attached, which is not the
            # case in some packages which imports (but not depends on)
            # mlogit 
            if (to_mlogit.data) data <- eval(mldata, parent.frame())
            else{
            # Construct a list of dfidx' arguments with the default values            
                dfa <- list(data = NA, idx = NULL, drop.index = TRUE, as.factor = NULL, pkg = NULL,
                            fancy.row.names = FALSE,
                            idnames = NULL, shape = "long", choice = NULL,
                            varying = NULL, sep = ".", opposite = NULL, levels = NULL, ranked = FALSE)
                # Replace the relevant values by those indicated in mlogit
                dfa[names(mldata)[-1]] <- as.list(mldata)[- 1]
                # Run dfidx with these values
                data <- dfidx::dfidx(data = data, dfa$idx, drop.index = dfa$drop.index,
                                     as.factor = dfa$as.factor,
                                     pkg = dfa$pkg, fancy.row.names = dfa$fancy.row.names,
                                     idnames = dfa$idnames, shape = dfa$shape,
                                     choice = dfa$choice, varying = dfa$varying,
                                     sep = dfa$sep, opposite = dfa$opposite,
                                     levels = dfa$levels, ranked = dfa$ranked)
            }
        }
    }
    # if data is an ordinary dfidx object, add the dfidx_mlogit class
    # if data is a data.table object, make it a data.frame object
    if (class(data)[1] == "dfidx") {
      class(data) <- c("dfidx_mlogit", class(data))
    } else if ("data.table" %in% class(data)) {
      data <- as.data.frame(data) }
  
    # 3 ######################################################
    # compute the model.frame
    ##########################################################

    mf <- callT
    mf$data <- data
    m <- match(c("formula", "data", "subset", "na.action",
                 "weights", "alt.subset", "reflevel"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    # use the unusual order of the arguments of model.frame, first
    # data, then formula
    mf$formula <- data
    mf$data <- formula
    mf[[1L]] <- as.name("model.frame")
    # mlogit needs balanced data
    mf$balanced <- TRUE
    mf <- eval(mf, parent.frame())

    # 4 ###########################################################
    # get the dimensions of the model
    ###############################################################
    .idx <- dfidx::idx(mf)
    alt <- dfidx::idx(mf, 2)
    chid <- dfidx::idx(mf, 1)
    alt.lev <- levels(alt[drop = TRUE])
    # model.frame_dfidx returns an alt.ordering attribute which
    # contains the levels of the second index in their initial order,
    # i.e. before redefining the reference level
    if (is.null(attr(mf, "alt.ordering"))) alt.ord <- levels(alt)
    else alt.ord <- attr(mf, "alt.ordering")
    altnoNA <- rep(alt.ord, nrow(.idx) / length(alt.lev))
    altnoNA <- factor(altnoNA, levels = alt.lev, labels = alt.lev)
    J <- length(alt.lev)
    N <- length(unique(chid))
    # 5 ###########################################################
    # extract the elements of the model
    ###############################################################

    # extract the weights if necessary
    if (any(names(mf) == "(weights)")){
        #YC à voir bizarre
        weights <- matrix(mf[["(weights)"]], ncol = J, byrow = TRUE)
        weights <- as.numeric(apply(weights, 1, min, na.rm = TRUE))
        weights <- weights / mean(weights)
    }
    else weights <- 1

    # extract the individual index if it is relevant
    if (panel){
        if (! mixed.logit) stop("panel is only relevant for mixed logit models")
        id <- dfidx::idx(mf, 1, 2)
        if (is.null(id)) stop("no individual index")
        # construct a data.frame with all the unique chids and the corresponding id
        #YC ????
        id <- unique(data.frame(chid, id))$id
    }
    else id <- NULL

    # extract the X matrix for the standard deviation estimation
    # (cformula is used to save the potentially 4 part formula)

    cformula <- formula
    if (length(formula)[2] == 4){
        #YC à revoir, incompréhensible
        # sformula is used to extract XS, formula is reduced to its
        # first 3 parts
        sformula <- formula(formula, rhs = 4)
        # ca semble marcher, mais pas top
        Xs <- model.matrix(sformula, as.data.frame(mf))
#        Xs <- model.matrix(sformula, mf)
        # only one line per choice situation is kept, one has to take
        # care to NA values due to the balancing of the data performed
        # by model.frame
        Xs <- na.omit(Xs)
        to.omit <- attr(Xs, "na.action")
        Xs <- Xs[! duplicated(chid[- to.omit]), - 1, drop = FALSE]
        # incompréhensible !!!!!!
        formula <- Formula(formula(formula, rhs = 1:3))
#        formula <- mFormula(formula(as.Formula(formula), rhs = 1:3))
    }
    else Xs <- NULL
    attr(mf, "formula") <- formula
    X <- model.matrix(mf, rhs = 1:3)
    formula <- cformula
    K <- ncol(X)
    df.residual <- N - K
    colnamesX <- colnames(X)

    # extract the response and coerce it to a logical
    y <- as.logical(model.response(mf))
    # choice is a vector of length N containing the alternative chosen
    choice <- na.omit(alt[y])
    # compute the choice frequency table
    freq <- table(alt[as.logical(y)])
    #YC should do the same, much simplier
    freq <- table(choice)
    # Xl and yl are lists of length J which contains N matrix / vector
    # of covariates and response; yv is a vector that contains the
    # chosen alternative
    Xl <- vector(length = J, mode = "list")
    names(Xl) <- levels(alt)
    for (i in levels(alt))  Xl[[i]] <- X[altnoNA == i, , drop = FALSE]
    yl <- split(y, altnoNA)
    yl <- lapply(yl, function(x){x[is.na(x)] <- FALSE ; x})
    unchid <- if(is.factor(chid)) as.character(levels(chid)) else unique(chid)
    attr(yl, "chid") <- unchid
    attr(yl, "id") <- as.character(levels(id))
    # for probit the response is a vector that contains the chosen
    # alternative as a numeric ; NA values of y are replaced by FALSE
    # so that the chosen alternative is correctly returned    

    if (probit){
        y2 <- y
        y2[is.na(y2)] <- FALSE
        yv <- as.numeric(alt[y2])
        # DX is a list of covariates first differenced with respect with the
        # chosen alternative
        DX <- vector("list", length = J - 1)
        DX <- lapply(DX, function(x) return(matrix(NA, nrow = N, ncol = K)))
        for (i in 1:N){
            any <- yv[i]
            j <- 1
            newj <- 1
            for (j in 1:J){
                if (j != any){
                    DX[[newj]][i, ] <- Xl[[j]][i, ] - Xl[[any]][i, ]
                    newj <- newj + 1
                }
            }
        }
    }

    # if the nests are defined by a grouping variable or by 'rpl',
    # compute the corresponding list of nests

    if (nested.logit){
        if (is.logical(nests) && nests){
            group_var <- dfidx::idx(mf, 2, 2)
            if (is.null(group_var)) stop("no grouping variable")
            else{
                nests <- unique(data.frame(alt = as.character(alt),
                                           group = group_var))
                nests <- split(nests$alt, nests$group)
            }
        }
        if (pair.comb.logit){
            alt1 <- rep(alt.lev, c((J - 1):0))
            alt2 <- alt.lev[unlist(lapply(2:J, function(x) x:J))]
            names.nests <- paste(alt1, alt2, sep = ".")
            nests <- mapply(function(x,y) c(x,y), alt1, alt2, SIMPLIFY = FALSE)
            names(nests) <- names.nests
        }
    }

    # 6 ######################################################
    # compute the starting values
    ##########################################################

    start <- mlogit.start(formula = formula, data = data, mf = mf, start = start,
                          un.nest.el = un.nest.el, nests = nests, heterosc = heterosc,
                          rpar = rpar, probit = probit, correlation = correlation,
                          alt.subset = alt.subset, reflevel = reflevel)
    names.sup.coef <- attr(start, "names.sup.coef")

    # 7 ###################################################################
    # Estimate the model using mlogit.optim and passing the correct arguments
    #######################################################################

    # construct the call for mlogit.optim
    opt <- callT

    # if constPar is numeric, insert the relevant value in the start
    # vector and transform the constPar vector to character
    if (! is.null(opt$constPar)){
        theconstPar <- eval(opt$constPar)
        if (is.numeric(theconstPar)){
            if (is.null(names(theconstPar)))
                stop('the numeric constPar vector should be named')
            start[names(theconstPar)] <- theconstPar
            opt$constPar <- names(theconstPar)
        }
    }

    # for the probit model, the first parameter of the covariance
    # matrix is not estimated
    if (probit) if (is.null(opt$constPar)) opt$constPar <- names.sup.coef[1]

    # include the automatically computed starting values
    opt$start <- start

    # select the argument of mlogit that should be passed to
    # mlogit.optim
    m <- match(c("method", "print.level", "iterlim",
                 "start", "constPar","tol", "ftol", "steptol"),
               names(opt), 0L)
    opt <- opt[c(1L, m)]
    opt[[1]] <- as.name('mlogit.optim')
    opt$logLik <- as.name('lnl.slogit')
    opposite <- - 1
    opt[c('weights', 'opposite')] <- list(as.name('weights'), as.name('opposite'))

    # model specific arguments
    if (! probit) opt[c('X', 'y')] <- list(as.name('Xl'), as.name('yl'))
    if (wlogit) opt[c('logLik', 'Xs')] <- list(as.name('lnl.wlogit'), as.name('Xs'))
    if (mixed.logit){        
        if (is.logical(correlation)){
            if (correlation) correlation <- names(rpar) else correlation <- character(0)
        }
        if (any (! names(rpar) %in% colnamesX)) stop("unknown random parameter")
        if (any (! correlation %in% names(rpar)))
            stop("unknown random parameter in the correlation vector")
        opt[c('logLik', 'id', 'rpar', 'correlation', 'halton')] <-
            list(as.name('lnl.rlogit'), as.name('id'), as.name('rpar'),
                 as.name('correlation'), as.name('halton'))
    }
    if (wlogit | mixed.logit)
        opt[c('Xs')] <- list(as.name('Xs'))
    if (probit | mixed.logit)
        opt[c('R', 'seed')] <- list(as.name('R'), as.name('seed'))
    if (probit)
        opt[c('logLik', 'X', 'y')] <- list(as.name('lnl.mprobit'), as.name('DX'), as.name('yv'))
    if (heterosc.logit){
        rn <- gauss.quad(R, kind = "laguerre")
        opt[c('logLik', 'rn')] <- list(as.name('lnl.hlogit'), as.name('rn'))
    }
    if (nested.logit)
        opt[c('logLik', 'nests', 'un.nest.el', 'unscaled')] <-
            list(as.name('lnl.nlogit'), as.name('nests'),
                 as.name('un.nest.el'), as.name('unscaled'))

    x <- eval(opt, sys.frame(which = nframe))
    
    # compute the probabilities for all the alternatives for
    # heteroscedastic and the probit model
    if (probit | heterosc){
        opt$logLik <- opt$iterlim <- opt$method <- opt$print.level <- opt$tol <- opt$constPar <- NULL
        names(opt)[[2]] <- 'param'
        opt[[2]] <- x$coefficients
        opt$gradient <- FALSE
        opt[[1]] <- as.name('lnl.hlogit')
        if (probit) opt[[1]] <- as.name('lnl.mprobit') else opt[[1]] <- as.name('lnl.hlogit')
        probabilities <- c()
        for (k in 1:J){
            if (probit){
                they <- rep(k, N)
                opt$y <- they
                DX <- vector("list", length = J - 1)
                DX <- lapply(DX, function(x) return(matrix(NA, nrow = N, ncol = K)))
                for (i in 1:N){
                    any <- k
                    j <- 1
                    newj <- 1
                    for (j in 1:J){
                        if (j != any){
                            DX[[newj]][i, ] <- Xl[[j]][i, ] - Xl[[any]][i, ]
                            newj <- newj + 1
                        }
                    }
                }
                opt$X <- DX
            }
            if (heterosc){
                they <- vector(mode = 'list', length = J)
                they <- lapply(they, function(x) rep(FALSE, N))
                they[[k]] <- rep(TRUE, N)
                names(they) <- names(yl)
                opt$y <- they
            }
            probabilities <- cbind(probabilities,
                                   attr(eval(opt, sys.frame(which = nframe)), 'fitted'))
        }
        colnames(probabilities) <- alt.lev
        attr(x$optimum, "probabilities") <- probabilities
    }

    # 8 ###########################################################
    # put the result in form
    ###############################################################

    # some general features
    #YC n/N already defined
    n <- sum(freq)
    x$est.stat$elaps.time <- proc.time() - start.time
    logLik <- structure( - as.numeric(x$optimum),
                        df = length(x$coefficients),
                        null = sum(freq * log(freq / N)),
                        class = "logLik"
                        )
    if (mixed.logit) rpar <- make.rpar(rpar, correlation, x$coefficients, NULL) else rpar <- NULL
#    if (! nested.logit) nests <- NULL

    # if no hessian is returned, use the BHHH approximation
    if (is.null(attr(x$optimum, 'hessian'))) hessian <- - crossprod(attr(x$optimum, 'gradi'))
    else hessian <- - attr(x$optimum, 'hessian')
    fitted <- attr(x$optimum, "fitted")
    probabilities <- attr(x$optimum, "probabilities")
    linpred <- attr(x$optimum, "linpred")
    resid <- Reduce("cbind", yl) - fitted
    attr(x$coefficients, "fixed") <- attr(x$optimum, "fixed")
    attr(x$coefficients, "sup") <- names.sup.coef
    gradient <- - attr(x$optimum, "gradi")
    if (mixed.logit) indpar <- attr(x$optimum, "indpar") else indpar <- NULL
    
    # Compute the covariance matrix of the errors
    if (multinom.logit | mixed.logit){
        alt.names <- colnames(probabilities)
        J <- length(alt.names)
        Omega <- matrix(0, J, J, dimnames = list(alt.names, alt.names))
        diag(Omega) <- pi ^ 2 / 6
    }
    if (probit){
        S <- matrix(0, J - 1, J - 1)
        S[! upper.tri(S)] <- x$coefficients[- c(1:K)]
        Mi <- list()
        M <- rbind(0, diag(J - 2))
        Mi[[1]] <- diag(J - 1)
        for (i in 2:(J - 1)){
            Mi[[i]] <- cbind(M[, 0:(i - 2), drop = FALSE], -1, M[, ((i - 1):(J - 2))])
        }
        Mi[[J]] <- cbind(M, -1)
        names(Mi) <- alt.lev
        Omega <- lapply(Mi, function(x) tcrossprod(x %*% S))
        for (i in 1:J){
            colnames(Omega[[i]]) <- rownames(Omega[[i]]) <- alt.lev[-i]
        }
    }
    if (nested.logit){
        alt.names <- colnames(probabilities)
        J <- length(alt.names)
        Omega <- matrix(0, J, J, dimnames = list(alt.names, alt.names))
        for (i in 1:length(nests)){
            anEl <- x$coefficients[paste('iv', names(nests)[i], sep = ".")]
            alts <- nests[[i]]
            M <- length(alts)
            for (m in 1:M){
                for (l in 1:M){
                    Omega[alts[m], alts[l]] <- (1 - anEl ^ 2)
                }
            }
        }
        diag(Omega) <- 1
        Omega <- Omega * tcrossprod(rep(pi / sqrt(6), J))
    }
    if (heterosc.logit){
        alt.names <- colnames(probabilities)
        J <- length(alt.names)
        Omega <- matrix(0, J, J, dimnames = list(alt.names, alt.names))
        diag(Omega) <- pi ^ 2 / 6
        for (i in 2:J){
            analt <- alt.names[i]
            Omega[analt, analt] <- pi ^ 2 / 6 * x$coefficients[paste('sp', analt, sep = '.')] ^ 2
        }
    }
    if (wlogit) Omega <- NA

    mf$probabilities <- as.numeric(t(probabilities))
    if (! is.null(linpred)) mf$linpred <- as.numeric(t(linpred))
    result <- structure(
        list(
            coefficients  = x$coefficients,
            logLik        = logLik,
            gradient      = gradient,
            hessian       = hessian,
            est.stat      = x$est.stat,
            fitted.values = fitted,
            probabilities = probabilities,
            linpred       = linpred,
            indpar        = indpar,
            residuals     = resid,
            omega         = Omega,
            rpar          = rpar,
            nests         = nests,
            model         = mf,
            freq          = freq,
            formula       = formula,
            call          = callT),
        class = 'mlogit'
    ) 
    result
}

# mlogit.model checks the arguments of mlogit and returns a vector of
# named booleans which caracterize the model
mlogit.model <- function(formula, nests = NULL, heterosc = FALSE, rpar = NULL, probit = FALSE){
    wlogit <- length(formula)[2] == 4
    heterosc.logit <- heterosc
    nested.logit <- ! is.null(nests)
    pair.comb.logit <- FALSE
    if (! is.null(nests) && length(nests) == 1){
        if (is.character(nests) && nests == "pcl"){
            nested.logit <- TRUE
            pair.comb.logit <- TRUE
        }
        if (is.logical(nests)){
            nested.logit <- nests
            pair.comb.logit <- FALSE
        }
    }
    mixed.logit <- ! is.null(rpar)
    multinom.logit <- ! heterosc & ! nested.logit & is.null(rpar) & ! probit & ! wlogit
    if (heterosc.logit + nested.logit + mixed.logit + probit > 1)
        stop("only one of heterosc, rpar, nests and probit can be used")
    mm <- c(multinom.logit, nested.logit, pair.comb.logit,
            probit, wlogit, heterosc.logit, mixed.logit)
    names(mm) <- c("multinom", "nested", "rpl", "probit", "wlogit", "heterosc", "mixed")
    mm
}

# mlogit.start compute the starting values, giving the arguments of
# mlogit
mlogit.start <- function(formula, data, mf, start = NULL, un.nest.el = FALSE,
                         nests = NULL, heterosc = FALSE, rpar = NULL,
                         probit = FALSE, correlation = FALSE,
                         alt.subset = NULL, reflevel = NULL){
    nframe <- length(sys.calls())
    callT <- match.call(expand.dots = TRUE)
    # get the model type using mlogit.model
    mt <- mlogit.model(formula, nests, heterosc, rpar, probit)
    multinom.logit <- mt['multinom']
    nested.logit <- mt['nested']
    pair.comb.logit <- mt['rpl']
    probit <- mt['probit']
    wlogit <- mt['wlogit']
    heterosc.logit <- mt['heterosc']
    mixed.logit <- mt['mixed']

    alt <- dfidx::idx(mf, 2)
    
    # extract the X and eventually Xs matrix
    if (length(formula)[2] == 4){
        sformula <- formula(as.Formula(formula), rhs = 4)
        Xs <- model.matrix(sformula, as.data.frame(mf))[! duplicated(dfidx::idx(mf, 1)), - 1]
#        Xs <- model.matrix(sformula, mf)[! duplicated(index(mf)$chid), - 1]
#        formula <- mFormula(formula(as.Formula(formula), rhs = 1:3))
        formula <- Formula(formula(formula, rhs = 1:3))
    }
    else Xs <- NULL
    alt.lev <- levels(alt)
    J <- length(alt.lev)
    attr(mf, "formula") <- formula
    X <- model.matrix(mf)
    K <- ncol(X)
    colnamesX <- colnames(X)
    # give names to the supplementary coefficients and values if start
    # is null or of length K
    sup.coef <- numeric(0)
    names.sup.coef <- character(0)
    if (nested.logit){
        if (is.logical(nests) && nests){
            if (is.null(dfidx::idx(mf, 2, 2))){
                stop("no grouping variable")
            }
            else{
                nests <- unique(data.frame(alt = as.character(dfidx::idx(mf, 2)),
                                           group = dfidx::idx(mf, 2, 2)))
                nests <- split(nests$alt, nests$group)
            }
        }
        L <- length(nests)
        # set the iv coef(s) to 1 if start contains K values (normal
        # coefficients)
        if (un.nest.el){
            if (is.null(start) || length(start) == K) sup.coef <- c(iv = 1)
            names.sup.coef <- 'iv'
        }
        else{
            if (is.null(start) || length(start) == K) sup.coef <- rep(1, L)
#            names.sup.coef <- paste(names(nests), 'iv', sep = ":")
            names.sup.coef <- paste('iv', names(nests), sep = ":")            
        }
    }
    if (pair.comb.logit){
        if (un.nest.el){
            if (is.null(start)) sup.coef <- c(iv = 1)
            names.sup.coef <- 'iv'
        }
        else{
            names.sup.coef <- NULL
            for (i in 1:(J - 1)){
                names.sup.coef <- c(names.sup.coef,
                                    paste('iv', alt.lev[i],
                                          alt.lev[(i + 1):J], sep = "."))
            }
            sup.coef <- rep(1, length(names.sup.coef))
        }
    }
    if (heterosc.logit){
        if (is.null(start) || length(start) == K) sup.coef <- rep(1, J - 1)
        names.sup.coef <- paste("sp", alt.lev[- 1], sep = ".")
    }
    if (mixed.logit){
        unknowndist <- rpar[! (rpar %in% c("cn", "ln", "tn", "n", "u", "t", "zbu", "zbt"))]
        if (length(unknowndist)){
            udstr <- paste("unknown distribution",
                           paste(unique(unknowndist), collapse = ", "))
            stop(udstr)
        }        
        names.rpar <- names(rpar)
        names.rpar.sig <- names.rpar[! rpar %in% c("zbu", "zbt")]
        names.fixed <- colnamesX[! colnamesX %in% names.rpar]
        # the names of the correlated and uncorrelated random
        # parameters in the order of the X matrix
        if (is.logical(correlation)){
            if (correlation) correlation <- names(rpar) else correlation <- character(0)
        }
        if (any (! names(rpar) %in% colnamesX)) stop("unknown random parameter")
        if (any (! correlation %in% names(rpar)))
            stop("unknown random parameter in the correlation vector")

        uncorrelated <- setdiff(names(rpar), correlation)
        fixedpar <- setdiff(colnamesX, names(rpar))
        singlepar <- names(rpar)[rpar %in% c("zbu", "zbt")]
        utwopars <- intersect(names(rpar)[! rpar %in% c("zbu", "zbt")], uncorrelated)       
        correlated <-   colnamesX[sort(match(correlation,  colnamesX))]
        uncorrelated <- colnamesX[sort(match(uncorrelated, colnamesX))]
        fixedpar <- colnamesX[sort(match(fixedpar, colnamesX))]
        randompar <- colnamesX[sort(match(names(rpar), colnamesX))]
        singlepar <- colnamesX[sort(match(singlepar, colnamesX))]
        utwopars <- colnamesX[sort(match(utwopars, colnamesX))]

        Kc <- length(correlated)
        Ku <- length(uncorrelated)
        Ko <- length(singlepar)
        
        names.sup.coef <- c()
        sup.coef <- c()
        if (Ku - Ko){
            if (is.null(start) || length(start) == K)
                sup.coef <- c(sup.coef, rep(0.1, Ku - Ko))
            names.sup.coef <- paste("sd", utwopars, sep = ".")
        }
        if (Kc){
            if (is.null(start) || length(start) == K)
                sup.coef <- c(sup.coef, rep(0.1, 0.5 * Kc * (Kc + 1)))
            names.sup.coef <- c(names.sup.coef,
                                names.rpar(correlated, prefix = "chol"))
        }
        if (is.null(start) || length(start) == K) names(sup.coef) <- names.sup.coef
    }
    if (probit){
        names.sup.coef <- c()
        for (i in 2:J){
            names.sup.coef <- c(names.sup.coef,
                                paste(alt.lev[i], alt.lev[i:J], sep = "."))
        }
        if (is.null(start) || length(start) == K){
            corrMat <- matrix(0.5, J - 1, J - 1)
            diag(corrMat) <- 1
            corrMat <- chol(corrMat)
            sup.coef <- (t(corrMat)[! upper.tri(corrMat)])
            names(sup.coef) <- names.sup.coef
        }
    }
    if (wlogit){
        if (is.null(start) || length(start) == K) sup.coef <- c(sup.coef, rep(0.01, ncol(Xs)))
        names.sup.coef <- c(names.sup.coef, paste("sig", colnames(Xs), sep = "."))
    }
    ## start can be :
    ##   1. NULL, in this case estimate the multinomial logit model,
    ##   2. a vector of length K ; then add starting values for the
    ##   supplementary coefficients,
    ##   3. a full set ; then just name the coefs.

    if (is.null(start)){
        callst <- callT
        callst$formula <- formula
        start <- rep(0, K)
        names(start) <- colnamesX
        callst$start <- start
        callst$print.level <- 0
        if (! multinom.logit){
            callst[c("nests", "rpar", "constPar", "iterlim")] <- NULL
            callst[c("heterosc", "panel", "probit")] <- FALSE
            callst[c("correlation", "un.nest.el")] <- FALSE
            callst$print.level <- 0
            if (length(callst$formula)[2] == 4)
#                callst$formula <- mFormula(formula(callst$formula, rhs = 1:3))
                callst$formula <- Formula(formula(callst$formula, rhs = 1:3))
                callst[[1]] <- as.name("mlogit")
            start <- coef(eval(callst, parent.frame()))
            if (mixed.logit){
                ln <- names(rpar[rpar == "ln"])
                start[ln] <- log(start[ln])
            }
        }
    }
    if (length(start) == K){
        names(start) <- colnamesX
        if (! multinom.logit){
            names(sup.coef) <- names.sup.coef
            if (probit) start <- start * sqrt(3) / pi
            start <- c(start, sup.coef)
        }
    }
    else{
        if (! multinom.logit) names(start) <- c(colnamesX, names.sup.coef)
        else names(start) <- colnamesX
    }
    structure(start, names.sup.coef = names.sup.coef)
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
            x[[i]][is.na(x[[i]])] <- 0
            s <- s+x[[i]]
        }
    }
    s
}
