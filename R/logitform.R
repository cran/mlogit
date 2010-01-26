has.intercept <- function(object, ...) {
  UseMethod("has.intercept")
}

has.intercept.default <- function(object, ...) {
  has.intercept(formula(object), ...)
}

has.intercept.formula <- function(object, ...) {
  attr(terms(object), "intercept") == 1L
}

has.intercept.Formula <- function(object, rhs = NULL, ...) {
  ## NOTE: return a logical vector of the necessary length
  ## (which might be > 1)
  if(is.null(rhs)) rhs <- 1:length(attr(object, "rhs"))
  sapply(rhs, function(x) has.intercept(formula(object, lhs = 0, rhs = x)))
}



## pFormula:
## methods : formula, model.frame, model.matrix, pmodel.response

logitform <- function(object){
  UseMethod("logitform")
}

is.logitform <- function(object)
  inherits(object, "logitform")

logitform.formula <- function(object){
  if (!inherits(object, "Formula")) object <- Formula(object)
  class(object) <- c("logitform", "Formula", "formula")
  object
}

expand.logitform <- function(x, data){
  ind.spec <- NULL
  rhs <- x[[3]]
  saveformula <- x
  if(length(rhs)>1 && rhs[[1]]=="|"){
    ind.spec <- paste("alt",":(",deparse(rhs[[3]]),")",sep="")
    alt.spec <- paste(deparse(rhs[[2]]))
    y <- as.character(x[[2]])
    alt.spec <- paste(alt.spec,"+",ind.spec)
    x <- as.formula(paste(y," ~ ",alt.spec,sep=""))
  }
  else{
    class(x) <- "formula"
  }
  x
}

logitform <- function(object) {
  stopifnot(inherits(object, "formula"))
  object <- Formula(object)
  class(object) <- c("logitform", class(object))
  object
}

as.Formula.logitform <- function(x, ...){
  class(x) <- class(x)[-1]
  x
}

model.frame.logitform <- function(formula, data, ..., lhs = NULL, rhs = NULL){
  if (is.null(rhs)) rhs <- 1:(length(formula)[2])
  if (is.null(lhs)) lhs <- ifelse(length(formula)[1]>0, 1, 0)
  index <- attr(data, "index")
  mf <- model.frame(as.Formula(formula), as.data.frame(data), ..., rhs = rhs)
  index <- index[rownames(mf),]
  index <- data.frame(lapply(index, function(x) x[drop = TRUE]), row.names = rownames(index))
  structure(mf,
            index = index,
            class = c("mlogit.data", class(mf)))
}

model.matrix.logitform <- function(object, data, ...){
  K <- length(data)
  omitlines <- attr(na.omit(data), "na.action")
  index <- attr(data, "index")
  alt <- index[["alt"]]
  chid <- index[["chid"]]
  has.int <- has.intercept(object)
  formula <- expand.logitform(object)
  if (has.int){
    upform <- as.formula(paste(".~", "alt", "+.+1"))
    formula <- update(formula, upform)
  }
  else{
    formula <- update(formula, .~.+1)
  }
  data$alt <- alt
  X <- model.matrix(formula, data)[, -1, drop = F]
  lev1 <- levels(alt)[1]
  motif <- paste("alt", lev1, ":", sep = "")
  varkeep <- regexpr(motif, colnames(X))<0
  X <- X[,varkeep]
  X[omitlines, ] <- NA
  X
}

has.intercept.logitform <- function(object, ...){
  attr(object, "class") <- "Formula"
  has.int <- has.intercept(object)
  ifelse(length(has.int) > 1, has.int[2], has.int[1])
}
