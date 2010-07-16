mlogit <- function(formula, data, subset, weights, na.action, start = NULL,
                   alt.subset = NULL, reflevel= NULL,
                   nests = NULL, un.nest.el = FALSE, unscaled = FALSE,
                   heterosc = FALSE, rpar = NULL,
                   R = 40, correlation = FALSE, halton = NULL, random.nb = NULL,
                   panel = FALSE, estimate = TRUE, ...){

  start.time <- proc.time()
  callT <- match.call(expand.dots = TRUE)
  callF <- match.call(expand.dots = FALSE)
  formula <- callF$formula <- mFormula(formula)
  nframe <- length(sys.calls())

  heterosc.logit <- heterosc
  nested.logit <- !is.null(nests)
  if (!is.null(nests) && length(nests) == 1 && nests == "pcl"){
    nested.logit <- FALSE
    pair.comb.logit <- TRUE
  }
  else pair.comb.logit <- FALSE
  
  mixed.logit <- !is.null(rpar)
  multinom.logit <- !heterosc & is.null(nests) & is.null(rpar) 

  if (heterosc.logit + nested.logit + mixed.logit > 1)
    stop("only one of heterosc, rpar and nests can be used")

  if (!multinom.logit && !is.null(callT$method) && callT$method == 'nr')
    stop("nr method only implemented for the simple multinomial logit mode")

  if (multinom.logit){
    if (is.null(callT$method)) callT$method <- 'nr'
    if (is.null(callT$print.level)) callT$print.level <- 0
  }
  
  # 1 ############################################################
  #  check whether arguments for mlogit.data are present: if so run
  #  mlogit.data
  ################################################################
  mldata <- callT
  response.name <- paste(deparse(attr(formula, "lhs")[[1]]))
  m <- match(c("data", "choice", "shape", "varying", "sep",
               "alt.var", "chid.var", "alt.levels",
               "opposite", "drop.index", "id"),
             names(mldata), 0L)
  use.mlogit.data <- sum(m[-1]) > 0
  if (use.mlogit.data){
    mldata <- mldata[c(1L, m)]
    mldata[[1L]] <- as.name("mlogit.data")
    mldata$choice <- response.name
#    data <- eval(mldata, parent.frame())
    data <- eval(mldata, sys.frame(which = nframe))
  }
  
  # 2 ###########################################################
  # model.frame (with the data provided or the one computed by
  # mlogit.data
  ###############################################################

  mf <- callT
  m <- match(c("formula", "data", "subset", "na.action", "weights"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$formula <- formula
  mf[[1L]] <- as.name("model.frame")
  if (use.mlogit.data) mf$data <- data
  # if the user called the data.frame "mldata", this conflicts with
  # the call. The following line seems to fix the bug
  mf$data <- data
  mf <- eval(mf, sys.frame(which = nframe))
  
  # change the reference level of the response if required
  if (!is.null(reflevel)){
    attr(mf, "index")[["alt"]] <-
      relevel(attr(mf, "index")[["alt"]], reflevel)
  }
  index <- attr(mf, "index")
  alt <- index[["alt"]]
  chid <- index[["chid"]]
  alt.lev <- levels(alt)

  if (panel){
    if (!mixed.logit) stop("panel is only relevant for mixed logit models")
    id <- index[["id"]]
    if (is.null(id)) stop("no individual index")
    id <- split(index[["id"]], alt)[[1]]
  }
  else id <- NULL
  # compute the relevent subset if required
  if (!is.null(alt.subset)){
    # we keep only choices that belong to the subset
    choice <- alt[model.response(mf)]
    choice <- choice %in% alt.subset
    unid <- unique(chid)
    names(choice) <- as.character(unid)
    id.kept <- choice[as.character(chid)]
    # we keep only the relevant alternatives
    alt.kept <- alt %in% alt.subset
    # the relevant subset for the data.frame and the indexes
    mf <- mf[id.kept & alt.kept, , drop = FALSE]
    alt <- alt[id.kept & alt.kept , drop = TRUE]
    chid <- chid[id.kept & alt.kept , drop = TRUE]
  }
  # balanced the data.frame i.e. insert rows with NA when an
  # alternative is not relevant
  alt.un <- unique(alt)
  chid.un <- unique(chid)
  n <- length(chid.un)
  T <- length(alt.un)
  balanced <- TRUE
  if (nrow(mf) != (n * T)){
    rownames(mf) <- paste(chid, alt, sep = ".")
    all.rn <- as.character(t(outer(chid.un, alt.un, paste, sep = ".")))
    mf <- mf[all.rn, ]
    rownames(mf) <- all.rn
    chid <- rep(chid.un, each = T)
    alt <- rep(alt.un, n)
    index <- data.frame(chid = chid, alt = alt, row.names = rownames(mf))
    balanced <- FALSE
  }
  #suppress individuals for which no choice is made
  if (FALSE){
    delete.id <- tapply(model.response(mf), chid, sum, na.rm = TRUE)
    delete.id <- names(delete.id[delete.id == 0])
    mf <- mf[!(chid %in% delete.id), ]
    index <- index[rownames(mf),]
    index[[1]] <- index[[1]][, drop=TRUE]
    index[[2]] <- alt <- index[[2]][, drop=TRUE]
    attr(mf, "index") <- index
  }
  # if estimate is FALSE, return the data.frame
  if (!estimate) return(mf)

  # 3 ###########################################################
  # extract the elements of the model
  ###############################################################

  y <- model.response(mf)
  choice <- alt[y]
  X <- model.matrix(formula, mf)
  K <- ncol(X)
  n <- nrow(X)
  df.residual <- n - K
  colnamesX <- colnames(X)
  if (any(names(mf)=="(weights)")){
    weights <- mf[["(weights)"]] <- mf[["(weights)"]]/mean(mf[["(weights)"]])
    weights <- split(weights, alt)[[1]]
  }
  else weights <- NULL
  freq <- table(alt[y])
  otime <- proc.time()

  X <- split(as.data.frame(X), alt)
  X <- lapply(X, as.matrix)
  y <- split(y, alt)
  y <- lapply(y, function(x){x[is.na(x)] <- FALSE ; x})
  
  if (mixed.logit){
    # if some random parameters are linked to factors, change the name
    # to the second level (only relevant for 2 levels factors (to be
    # tested)
    for (i in 1:length(rpar)){
      n <- names(rpar)[i]
      clvar <- class(data[[n]])
      if (inherits(clvar, "factor")){
        lvar <- levels(data[[n]])[2]
        names(rpar)[i] <- paste(n, lvar, sep="")
      }
    }
    Vara <- sort(match(names(rpar), colnames(X[[1]])))
    Varc <- (1:K)[-Vara]
    Ka <- length(Vara)
    Xa <- lapply(X, function(x) x[, Vara, drop = F])
    Xc <- lapply(X, function(x) x[, Varc, drop = F])
    # create the random numbers matrix
    if (is.null(random.nb)) random.nb <- make.random.nb(R, Ka, halton)
    colnames(random.nb) <- colnames(Xa[[1]])
  }

  # 4 ###########################################################
  # Compute the starting values
  ###############################################################

  # first give names to the supplementary coefficients and values if
  # start is null

  if (nested.logit){
    J <- length(nests)
    if (un.nest.el){
      if (is.null(start) || length(start) == K) sup.coef <- c(iv = 1)
      names.sup.coef <- 'iv'
    }
    else{
      if (is.null(start) || length(start) == K) sup.coef <- rep(1, J)
      names.sup.coef <- paste("iv", names(nests), sep = ".")
    }
  }
  if (pair.comb.logit){
    unalt <- levels(alt)
    J <- length(unalt)
    if (un.nest.el){
      if (is.null(start)) sup.coef <- c(iv = 1)
      names.sup.coef <- 'iv'
    }
    else{
      names.sup.coef <- NULL
      for (i in 1:(J-1)){
        names.sup.coef <- c(names.sup.coef, paste('iv', unalt[i], unalt[(i+1):J], sep = "_"))
      }
      sup.coef <- rep(1, length(names.sup.coef))
    }
  }
  
  if (heterosc.logit){
    unalt <- levels(alt)
    J <- length(unalt)
    if (is.null(start) || length(start) == K) sup.coef <- rep(1, J-1)
    names.sup.coef <- paste("sp", unalt[-1], sep = ".")
  }
  if (mixed.logit){
    nmean <- length(c(Varc, Vara))
    nvar <- length(Vara)
    if (!correlation){
      if (is.null(start) || length(start) == K) sup.coef <- rep(.1, nvar)
      names.sup.coef <- paste("sd", colnames(Xa[[1]]), sep = ".")
    }
    else{
      if (is.null(start) || length(start) == K) sup.coef <- rep(.1, 0.5 * nvar * (nvar + 1))
      names.sup.coef <- c()
      Ka <- length(rpar)
      for (i in 1:Ka){
        names.sup.coef <- c(names.sup.coef,
                            paste(names(rpar)[i], names(rpar)[i:Ka], sep = "."))
      }
    }
  }
  
## if no starting values are provided, provide a 0 vector and, if the
## model is not the multinomial logit, estimate the model
##   start can be :
##   1. NULL, in this case estimate the multinomial logit model,
##   2. a vector of lengthK ; then add starting values for the
##   supplementary coefficients,
##   3. a full set ; then just name the coefs.

  if (is.null(start)){
    callst <- callT
    start <- rep(0, K)
    names(start) <- colnames(X[[1]])
    callst$start <- start
    callst$print.level <- 0
    if (!multinom.logit){
      callst$nests <- NULL
      callst$heterosc <- FALSE
      callst$rpar <- NULL
      callst$panel <- FALSE
      callst$constPar <- NULL
      callst$iterlim <- NULL
      start <- coef(eval(callst, sys.frame(which=nframe)))
      if (mixed.logit){
        ln <- names(rpar[rpar == "ln"])
        start[ln] <- log(start[ln])
      }
    }
  }
  if (length(start) == K){
    names(start) <- colnames(X[[1]])
    if (!multinom.logit){
      names(sup.coef) <- names.sup.coef
      start <- c(start, sup.coef)
    }
  }
  else{
    if (!multinom.logit) names(start) <- c(colnames(X[[1]]), names.sup.coef)
    else names(start) <- colnames(X[[1]])
  }

  # 5 ###################################################################
  # Estimate the model using mlogit.nlm and passing the correct arguments
  #######################################################################

  opt <- callT
  opt$start <- start
  m <- match(c("method", "print.level", "iterlim",
               "start", "constPar","tol", "ftol", "steptol"),
             names(opt), 0L)
  opt <- opt[c(1L, m)]
  opt[[1]] <- as.name('mlogit.optim')
  opt$logLik <- as.name('lnl.mlogits')
  if (!mixed.logit) opt[c('X', 'y')] <- list(as.name('X'), as.name('y'))
  if (mixed.logit){
    opt$logLik <- as.name('lnl.rlogit')
    opt[c('Xa', 'Xc', 'y', 'Varc', 'Vara', 'random.nb', 'id', 'rpar', 'correlation')] <-
    list(as.name('Xa'), as.name('Xc'), as.name('y'), as.name('Varc'), as.name('Vara'),
         as.name('random.nb'), as.name('id'), as.name('rpar'), as.name('correlation'))
  }
  if (heterosc.logit){
    opt$logLik <- as.name('lnl.hlogit')
    rn <- gauss.quad(R, kind = "laguerre")
    opt[c('rn', 'choice')] <- list(as.name('rn'), as.name('choice'))
  }
  if (nested.logit){
    opt$logLik <- as.name('lnl.nlogit')
    opt$nests <- as.name('nests')
    opt$un.nest.el <- as.name('un.nest.el')
    opt$unscaled <- as.name('unscaled')
  }
  if (pair.comb.logit){
    opt$logLik <- as.name('lnl.nlogit')
    alt.lev <- levels(alt)
    J <- length(alt.lev)
    alt1 <- rep(alt.lev, c((J-1):0))
    alt2 <- alt.lev[unlist(lapply(2:J, function(x) x:J))]
    names.nests <- paste(alt1, alt2, sep = "")
    lnests <- mapply(function(x,y) c(x,y), alt1, alt2, SIMPLIFY = FALSE)
    names(lnests) <- names.nests
    opt$nests <- lnests
    opt$unscaled <- as.name('unscaled')
    opt$un.nest.el <- as.name('un.nest.el')
  }
  if (!is.null(weights)) opt$weights <- as.name('weights')
  opt$opposite <- TRUE

  x <- eval(opt, sys.frame(which = nframe))

  # 6 ###########################################################
  # put the result in form
  ###############################################################

  n <- sum(freq)
  logLik <- - as.numeric(x$optimum)
  attr(logLik,"df") <- length(x$coefficients)
  attr(logLik, 'null') <- sum(freq*log(freq/n))
  class(logLik) <- "logLik"
  fit <- as.matrix(data.frame(attr(x$optimum, 'probabilities')))
  
  if (heterosc){
    opt$gradient <- FALSE
    opt$logLik <- opt$iterlim <- opt$method <- opt$print.level <- opt$tol <- NULL
    opt[[1]] <- as.name('lnl.hlogit')
    names(opt)[[2]] <- 'param'
    fit <- c()
    for (j in 1:J){
      they <- vector(mode='list', length= J)
      they <- lapply(they, function(x) rep(FALSE, n))
      they[[j]] <- rep(TRUE, n)
      opt$y <- they
      fit <- cbind(fit,
                   attr(eval(opt, sys.frame(which = nframe)), 'probabilities'))
    }
    colnames(fit) <- alt.lev
  }
  if (!nested.logit)  resid <- as.matrix(data.frame(y)) - fit
  else resid <- NULL

    
  # if no hessian is returned, use the BHHH approximation
  if (is.null(attr(x$optimum, 'hessian')))
    hessian <- - crossprod(attr(x$optimum, 'gradi'))
  else hessian <- - attr(x$optimum, 'hessian')
  x$est.stat$elaps.time <- proc.time() - start.time

  if (mixed.logit)   rpar <- make.rpar(rpar, correlation, x$coefficients, NULL)
  else rpar <- NULL

  thecoef <- x$coefficients
  attr(thecoef, "fixed") <- attr(x$optimum, "fixed")
  structure(
            list(
                 coefficients  = thecoef,
                 logLik        = logLik,
                 gradient      = - attr(x$optimum, 'gradient'),
                 gradi         = - attr(x$optimum, 'gradi'),
                 hessian       = hessian,
                 est.stat      = x$est.stat,
                 fitted.values = fit,
                 residuals     = resid,
                 rpar          = rpar,
                 model         = mf,
                 freq          = freq,
                 formula       = formula,
                 call          = callT),
            class = 'mlogit'
            )
}

