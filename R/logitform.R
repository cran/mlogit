##-------------------------------
## methods for logitform objects |
##-------------------------------
##    * terms.logitform          |
##    * update.logitform         |
##    * is.logitform             |
##    * logitform.formula        |
##    * logitform.list           |
##    * model.matrix.logitform   |
##    * model.frame.logitform    |
##-------------------------------

terms.logitform <- function(x,...){
  rhs <- x[[3]]
  saveformula <- x
  if (length(rhs==1)){
    alt.spec <- paste(deparse(rhs))
    ind.spec <- NULL
  }
  if(length(rhs)>1 && rhs[[1]]=="|"){
    ind.spec <- deparse(rhs[[3]])
    alt.spec <- paste(deparse(rhs[[2]]))
  }
    y <- as.character(x[[2]])
  x <- list(y=y,ind.spec=ind.spec,alt.spec=alt.spec)
  x
}

update.logitform <- function(object, new,...){
  old <- object
  if (!is.logitform(old)) old <- logitform(old)
  if (!is.logitform(new)) new <- logitform(new)
  to <- terms(old)
  tn <- terms(new)
  if (tn$y!=".") to$y=tn$y
  if (tn$alt.spec!="."){
    old.alt.spec <- formula(paste("~",to$alt.spec))
    new.alt.spec <- formula(paste("~",tn$alt.spec))
    alt.spec <- update(old.alt.spec,new.alt.spec)
    to$alt.spec <- paste(deparse(alt.spec[[2]]))
  }
  if (!is.null(tn$ind.spec) && tn$ind.spec!="."){
    new.ind.spec <- formula(paste("~",tn$ind.spec))
    if (is.null(to$ind.spec)){
      ind.spec <- new.ind.spec
    }
    else{
      old.ind.spec <- formula(paste("~",to$ind.spec))
      ind.spec <- update(old.ind.spec,new.ind.spec)
    }
    to$ind.spec <- paste(deparse(ind.spec[[2]]))
    if (to$ind.spec=="1") to$ind.spec <- NULL
  }
#  if (to$alt.spec=="1") to$alt.spec <- "0"
  logitform(to)
  
}

is.logitform <- function(object){
  ifelse(class(object)[1]=="logitform",TRUE,FALSE)
}

logitform <- function(object){
  UseMethod("logitform")
}

logitform.formula <- function(object){
  class(object) <- c("logitform","formula")
  object
}

logitform.list <- function(object){
  if (is.null(object$y)) y <- "" else y <- object$y
  if (is.null(object$alt.spec)) alt.spec <- 0 else alt.spec <- object$alt.spec
  if (is.null(object$ind.spec)){
    x <- paste(y,"~",alt.spec)
  }
  else{
    ind.spec <- object$ind.spec
    x <- paste(y,"~",alt.spec,"|",ind.spec)
  }
  x <- formula(x)
  logitform(x)
}

model.matrix.logitform <- function(object, data, ...){
  alt.name <- names(data)[2]
  formula <- expand.formula(object, data)
  if (attr(terms(formula),"intercept")==1){
    upform <- as.formula(paste(".~",alt.name,"+.+1"))
    formula <- update(formula,upform)
  }
  else{
    formula <- update(formula,.~.+1)
  }
  mf <- model.frame(formula,data)
  X <- model.matrix(formula,data)[,-1,drop=F]
  lev1 <- levels(data[[alt.name]])[1]
  motif <- paste(alt.name,lev1,":",sep="")
  varkeep <- regexpr(motif,colnames(X))<0
  X <- X[,varkeep]
  X
}

model.frame.logitform <- function(formula, data, ...){
  indexes <- data[,c(1:2)]
  formula <- expand.formula(formula, data)
  class(formula) <- "formula"
  mf <- model.frame(formula,data, ...)
  terms.mf <- attr(mf,"terms")
  selected.rows <- intersect(rownames(mf),rownames(indexes))
  mf <- cbind(indexes[selected.rows,],mf[selected.rows,,drop = FALSE])
  attr(mf,"terms") <- terms.mf
  mf
}
  
