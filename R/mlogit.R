mlogit <- function(formula, data, subset, weights = NULL, na.action, alt.subset = NULL, reflevel= NULL, ...){
  class(formula) <- c("logitform","formula")
  expanded.formula <- expand.formula(formula,data)
  y.name <- as.character(formula[[2]])
  alt.name <- names(data)[2]
  chid.name <- names(data)[1]
  if (!is.null(alt.subset)){
    choice <- data[[2]][data[[y.name]]]
    choice <- choice %in% alt.subset
    id <- unique(data[[1]])
    names(choice) <- as.character(id)
    id.kept <- choice[as.character(data[[1]])]
    alt.kept <- data[[2]] %in% alt.subset
    data <- data[id.kept & alt.kept,]
#    data <- subset(data,id.kept & alt.kept)
    data[[2]] <- data[[2]][,drop=TRUE]
  }
  if (!is.null(reflevel)){
    data[[alt.name]] <- relevel(data[[alt.name]],reflevel)
  }
  start.time <- proc.time()
  cl <- match.call()
  cl$formula <- formula
  # model.frame
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action","weights"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf$formula <- do.call("logitform",list(mf$formula))
  mf <- eval(mf, parent.frame())
  
  # balanced the data.frame i.e. insert rows with NA when an
  # alternative is not relevant
  alt.un <- unique(data[[alt.name]])
  chid.un <- unique(data[[chid.name]])
  n <- length(chid.un)
  T <- length(alt.un)
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
  
  result <- maxLik(f,g,h,start=rep(0,ncol(X[[1]])), ...)
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
  result <- list(coefficients = coef, logLik = logLik,
                 hessian = hessian, gradient = gradient,
                 call = cl, est.stat = est.stat, freq = freq,
                 residuals = residuals, fitted.values = fitted.values,
                 formula = formula, expanded.formula=expanded.formula, model= mf, index=data[1:2])
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

mlogit.data <- function(x, choice, shape = c("wide","long"), varying = NULL, sep = ".",
                         alt.var = NULL, id.var = NULL, alt.levels = NULL, opposite = NULL, ...){
  if (shape=="long"){
    if (is.null(id.var)){
      chid.name <- "chid"
      chid.is.variable <- FALSE
    }
    else{
      chid.name <- id.var
      chid.is.variable <- ifelse(is.null(x[[id.var]]),FALSE,TRUE)
    }
    choice.name <- choice
    choice <- x[[choice]]
    if (is.null(alt.var) && is.null(alt.levels)) stop("at least one of alt.var and alt.levels should be filled")
    
    if (!is.null(alt.levels)){
      J <- length(alt.levels)
      n <- nrow(x)/J
      alt <- factor(rep(alt.levels,n),levels=alt.levels)
      if (!is.null(alt.var) && !is.null(x[[alt.var]])){
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
      if (!is.factor(x[[alt.name]])) x[[alt.name]] <- factor(x[[alt.name]])
      alt.levels <- levels(x[[alt.name]])
      J <- length(alt.levels)
      alt <- x[[alt.name]]
    }
    n <- nrow(x)/J
    if (!chid.is.variable) choiceid <- rep(1:n,each=J) else choiceid <- x[[chid.name]]
    if (!is.logical(x[[choice.name]])) x[[choice.name]] <- as.logical(x[[choice.name]])
    x <- cbind(ch=choiceid,alt=alt,x)
    if (alt.is.variable){
      i <- which(names(x)==alt.name)
      x <- x[-i]
    }
    if (chid.is.variable){
      i <- which(names(x)==chid.name)
      x <- x[-i]
    }
    names(x)[1:2] <- c(chid.name,alt.name)
    row.names(x) <- paste(x[[chid.name]],x[[alt.name]],sep=".")
  }


  
  if (shape == "wide"){
    if (is.null(alt.var)) alt <- "alt" else alt <- alt.var
    if (is.null(id.var)) chid <- "chid" else chid <- id.var
    x <- reshape(x,varying = varying, direction = "long", sep = sep, timevar = alt, idvar = chid,  ...)
    x <- x[order(x[[chid]],x[[alt]]),]
    id <- x[[chid]]
    time <- x[[alt]]
    idpos <- which(names(x)==chid)
    timepos <- which(names(x)==alt)
    x <- x[,-c(idpos,timepos)]
    id <- as.factor(id)
    time <- as.factor(time)
    x <- cbind(id,time,x)
    names(x)[1:2] <- c(chid,alt)
    if (!is.null(alt.levels)){
      levels(x[[choice]]) <- alt.levels
      levels(x[[alt]]) <- alt.levels
      row.names(x) <- paste(x[[chid]],x[[alt]],sep=".")
    }
    x[[choice]] <- x[[choice]]==x[[alt]]
  }

  if (!is.null(opposite)){
    for (i in opposite){
      x[[i]] <- -x[[i]]
    }
  }
  x[[1]] <- factor(x[[1]])
  x[[2]] <- factor(x[[2]])
  x
}
