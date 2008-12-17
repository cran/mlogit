mlogit <- function(formula, data, subset, weights = NULL, na.action, alt.subset = NULL, reflevel= NULL, ...){
  
  # first check whether arguments for mlogit.data are present: if so run mlogit.data
  mf <- match.call(expand.dots = TRUE)
  m <- match(c("data","choice","shape","varying","sep","alt.var","id.var","alt.levels","opposite"),
             names(mf), 0L)
  use.mlogit.data <- sum(m[-1]) > 0
  nframe <- length(sys.calls())
  if (use.mlogit.data){
    mf <- mf[c(1L, m)]
    mf[[1L]] <- as.name("mlogit.data")
    data <- eval(mf, parent.frame())
    new.data.name <- "mydata"
    assign(new.data.name,data,env=sys.frame(which=nframe))
  }
  
  class(formula) <- c("logitform","formula")
  expanded.formula <- expand.formula(formula,data)
  y.name <- as.character(formula[[2]])
  alt.name <- names(data)[2]
  chid.name <- names(data)[1]
  start.time <- proc.time()
  cl <- match.call()
  cl$formula <- formula

  # model.frame (with the data provided or the one computed by mlogit.data
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action","weights"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf$formula <- do.call("logitform",list(mf$formula))
  if (use.mlogit.data){
    mf$data <- as.name(new.data.name)
    mf <- eval(mf, sys.frame(which=nframe))
  }
  else{
    mf <- eval(mf, parent.frame())
  }

  if (!is.null(reflevel)){
    mf[[alt.name]] <- relevel(mf[[alt.name]],reflevel)
  }

  if (!is.null(alt.subset)){
    choice <- mf[[2]][mf[[y.name]]]
    choice <- choice %in% alt.subset
    id <- unique(mf[[1]])
    names(choice) <- as.character(id)
    id.kept <- choice[as.character(mf[[1]])]
    alt.kept <- mf[[2]] %in% alt.subset
    mf <- mf[id.kept & alt.kept,]
    mf[[2]] <- mf[[2]][,drop=TRUE]
  }
  else{
    choice <- mf[[2]]
  }

  
  # balanced the data.frame i.e. insert rows with NA when an
  # alternative is not relevant

##   alt.un <- unique(data[[alt.name]])
##   chid.un <- unique(data[[chid.name]])
##   n <- length(chid.un)
##   T <- length(alt.un)
##   rownames(mf) <- paste(data[[chid.name]],data[[alt.name]],sep=".")
##   all.rn <- as.character(t(outer(chid.un,alt.un,paste,sep=".")))
##   mf <- mf[all.rn,]
##   rownames(mf) <- all.rn
##   mf[[1]] <- rep(chid.un,each=T)
##   mf[[2]] <- rep(alt.un,n)

  alt.un <- unique(mf[[alt.name]])
  chid.un <- unique(mf[[chid.name]])
  n <- length(chid.un)
  T <- length(alt.un)
  rownames(mf) <- paste(mf[[chid.name]],mf[[alt.name]],sep=".")
  all.rn <- as.character(t(outer(chid.un,alt.un,paste,sep=".")))
  mf <- mf[all.rn,]
  rownames(mf) <- all.rn
  mf[[1]] <- rep(chid.un,each=T)
  mf[[2]] <- rep(alt.un,n)

  
  #suppress individuals for which no choice is made
  y <- mf[[y.name]]
  delete.id <- tapply(y,mf[[1]],sum,na.rm=TRUE)
  delete.id <- names(delete.id[delete.id==0])
  mf <- mf[!(mf[[1]]%in%delete.id),]
  mf[[1]] <- mf[[1]][,drop=TRUE]

  
  y <- mf[[y.name]]
  X <- model.matrix(formula,mf)

  namesX <- colnames(X)
  if (any(names(mf)=="(weights)")) mf[["(weights)"]] <- mf[["(weights)"]]/mean(mf[["(weights)"]])

  alt <- mf[[alt.name]]
#  chid <- as.character(mf[[chid.name]])
  chid <- mf[[chid.name]]
  mframe <- mf

  freq <- table(alt[y])

  X <- split(as.data.frame(X),alt)
  X <- lapply(X,as.matrix)
  y <- split(y,alt)


#  n <- length(y[[1]])
#  K <- ncol(X[[1]])

  rownames.X <- split(chid,alt)

  y <- lapply(y,function(x){x[is.na(x)] <- FALSE;x})
  X <- lapply(X,function(x){x[is.na(x)] <- 0;x})
  
  f <- function(param) mlogit.lnl(param, X, y, weights = NULL)
  g <- function(param) mlogit.grad(param, X, y, weights = NULL)
  h <- function(param) mlogit.hess(param, X, y, weights = NULL)

  mf <- match.call(expand.dots = TRUE)
  m <- match(c("method","print.level","iterlim","start","constPar","activePar"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("maxLik")
  mf[c('logLik', 'grad', 'hess')] = c(as.name('f'),as.name('g'),as.name('h'))
  if (is.null(mf$start)) mf$start <- rep(0,ncol(X[[1]]))
  result <- eval(mf,sys.frame(which=nframe))

#  choice <- data[[choice.name]]
  eff <- table(alt[choice])
  n <- sum(eff)
  logLik0 <- sum(eff*log(eff/n))
  
#  result <- maxLik(f,g,h,start=rep(0,ncol(X[[1]])))
  gradient <- g(result$estimate)
  coef <- result$estimate
  P <- mlogit.P(coef,X)
  fitted.values <- as.matrix(as.data.frame(P))
  residuals <- as.matrix(as.data.frame(y))-fitted.values
  logLik <- result$maximum
  attr(logLik,"df") <- length(coef)
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
  result <- list(coefficients = coef, logLik = logLik, logLik0 = logLik0,
                 hessian = hessian, gradient = gradient,
                 call = cl, est.stat = est.stat, freq = freq,
                 residuals = residuals, fitted.values = fitted.values,
                 formula = formula, expanded.formula=expanded.formula, model= mframe, index=mf[1:2])
  class(result) <- "mlogit"
  result
}

mlogit.P <- function(param,X){
  eXb <- lapply(X,function(x) exp(crossprod(t(x),param)))
#  eXb <- lapply(eXb,function(x){x[is.na(x)] <- 0;x})
  seXb <- suml(eXb)
  P <- lapply(eXb,function(x){ v <- x/seXb; as.vector(v)})
  P
}
  
mlogit.lnl <- function(param, X, y, weights = NULL){
  if (is.null(weights)) weights <- 1
  P <- mlogit.P(param,X)
#  y <- lapply(y,function(x){x[is.na(x)] <- FALSE;x})
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

expand.formula <- function(formula, data){
  alt.name <- names(data)[2]
  ind.spec <- NULL
  rhs <- formula[[3]]
  saveformula <- formula
  if(length(rhs)>1 && rhs[[1]]=="|"){
#    ind.spec <- paste("alt:(",deparse(rhs[[3]]),")",sep="")
    ind.spec <- paste(alt.name,":(",deparse(rhs[[3]]),")",sep="")
    alt.spec <- paste(deparse(rhs[[2]]))
    y <- as.character(formula[[2]])
    alt.spec <- paste(alt.spec,"+",ind.spec)
    formula <- as.formula(paste(y," ~ ",alt.spec,sep=""))
  }
  else{
    class(formula) <- "formula"
  }
  formula
}

suml <- function(x){
  n <- length(x)
  if (!is.null(dim(x[[1]]))){
    d <- dim(x[[1]])
    s <- matrix(0,d[1],d[2])
    for (i in 1:n){
      s <- s+x[[i]]
#      s <- sum(c(s,x[[i]]),na.rm=na.rm)
    }
  }
  else{
    s <- rep(0,length(x[[n]]))
    for (i in 1:n){
      s <- s+x[[i]]
#      s <- sum(c(s,x[[i]]),na.rm=na.rm)
    }
  }
  s
}

mlogit.data <- function(data, choice, shape = c("wide","long"), varying = NULL, sep = ".",
                         alt.var = NULL, id.var = NULL, alt.levels = NULL, opposite = NULL, ...){
  if (shape=="long"){
    if (is.null(id.var)){
      chid.name <- "chid"
      chid.is.variable <- FALSE
    }
    else{
      chid.name <- id.var
      chid.is.variable <- ifelse(is.null(data[[id.var]]),FALSE,TRUE)
    }
    choice.name <- choice
    choice <- data[[choice]]

    if (is.null(alt.var) && is.null(alt.levels)) stop("at least one of alt.var and alt.levels should be filled")
    
    if (!is.null(alt.levels)){
      J <- length(alt.levels)
      n <- nrow(data)/J
      alt <- factor(rep(alt.levels,n),levels=alt.levels)
      if (!is.null(alt.var) && !is.null(data[[alt.var]])){
        warning(paste("variable",alt.var,"exists and will be replaced"))
        alt.is.variable <- TRUE
      }
      else{
        alt.is.variable <- FALSE
      }
      alt.name <- ifelse(is.null(alt.var),"alt",alt.var)
    }
    else{
      alt.name <- alt.var
      alt.is.variable <- TRUE
      if (!is.factor(data[[alt.name]])) data[[alt.name]] <- factor(data[[alt.name]])
      alt.levels <- levels(data[[alt.name]])
      J <- length(alt.levels)
      alt <- data[[alt.name]]
    }
    n <- nrow(data)/J
    if (!chid.is.variable) choiceid <- rep(1:n,each=J) else choiceid <- data[[chid.name]]

    if (!is.logical(data[[choice.name]])){
      if (is.factor(choice) && 'yes' %in% levels(choice)) data[[choice.name]] <- data[[choice.name]] == 'yes'
      if (is.numeric(choice)) data[[choice.name]] <- data[[choice.name]] != 0
    }
    data <- cbind(ch=choiceid,alt=alt,data)
    if (alt.is.variable){
      i <- which(names(data)==alt.name)
      data <- data[-i]
    }
    if (chid.is.variable){
      i <- which(names(data)==chid.name)
      data <- data[-i]
    }
    names(data)[1:2] <- c(chid.name,alt.name)
    row.names(data) <- paste(data[[chid.name]],data[[alt.name]],sep=".")
  }


  
  if (shape == "wide"){
    class(data[[choice]]) <- 'factor'
    if (is.null(alt.var)) alt <- "alt" else alt <- alt.var
    if (is.null(id.var)) chid <- "chid" else chid <- id.var
    if (!is.null(varying)){
      data <- reshape(data,varying = varying, direction = "long", sep = sep, timevar = alt, idvar = chid,  ...)
    }
    else{
      id.names <- as.numeric(rownames(data))
      nb.id <- length(id.names)
      data[[chid]] <- id.names
      lev.ch <- levels(data[[choice]])
      data <- data.frame(lapply(data,rep,length(lev.ch)))
      data[[alt]] <- rep(lev.ch,each=nb.id)
      row.names(data) <- paste(data[[chid]],data[[alt]],sep=".")
#      rownames(data) <- paste(outer(data[[chid]],data[[alt]],sep="."))
    }
    data <- data[order(data[[chid]],data[[alt]]),]
    id <- data[[chid]]
    time <- data[[alt]]
    idpos <- which(names(data)==chid)
    timepos <- which(names(data)==alt)
    data <- data[,-c(idpos,timepos)]
    id <- as.factor(id)
    time <- as.factor(time)
    data <- cbind(id,time,data)
    names(data)[1:2] <- c(chid,alt)
    if (!is.null(alt.levels)){
      levels(data[[choice]]) <- alt.levels
      levels(data[[alt]]) <- alt.levels
      row.names(data) <- paste(data[[chid]],data[[alt]],sep=".")
    }
    data[[choice]] <- data[[choice]]==data[[alt]]
  }

  if (!is.null(opposite)){
    for (i in opposite){
      data[[i]] <- -data[[i]]
    }
  }
  data[[1]] <- factor(data[[1]])
  data[[2]] <- factor(data[[2]])
  data
}
