mlogit.data <- function(data, choice, shape = c("wide","long"), varying = NULL,
                        sep = ".", alt.var = NULL, chid.var = NULL, 
                        alt.levels = NULL, id.var = NULL, opposite = NULL,
                        drop.index = FALSE, ...){
  # chid.name, alt.name : the name of the index variables
  # chid, alt : the index variables
  if (shape == "long"){
    if (is.null(chid.var)){
      chid.name <- "chid"
      chid.is.variable <- FALSE
    }
    else{
      chid.name <- chid.var
      chid.is.variable <- ifelse(is.null(data[[chid.var]]), FALSE, TRUE)
    }
    choice.name <- choice
    choice <- data[[choice]]

    if (is.null(alt.var) && is.null(alt.levels))
      stop("at least one of alt.var and alt.levels should be filled")
    
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
    if (!chid.is.variable) chid <- rep(1:n, each = J) else chid <- data[[chid.name]]
    if (!is.logical(data[[choice.name]])){
      if (is.factor(choice) && 'yes' %in% levels(choice))
        data[[choice.name]] <- data[[choice.name]] == 'yes'
      if (is.numeric(choice)) data[[choice.name]] <- data[[choice.name]] != 0
    }
    # remplacer id par chid à gauche

    chid <- as.factor(chid)
    alt <- as.factor(alt)
    row.names(data) <- paste(chid, alt, sep = ".")
  }
  
  if (shape == "wide"){
    if (is.ordered(data[[choice]])) class(data[[choice]]) <- "factor"
    else data[[choice]] <- as.factor(data[[choice]])
    # this doesn't work for ordered factors which remains ordered
    if (is.null(alt.var)) alt.name <- "alt" else alt.name <- alt.var
    if (is.null(chid.var)) chid.name <- "chid" else chid.name <- chid.var
    if (!is.null(varying)){
      data <- reshape(data, varying = varying, direction = "long", sep = sep,
                      timevar = alt.name, idvar = chid.name,  ...)
    }
    else{
      id.names <- as.numeric(rownames(data))
      nb.id <- length(id.names)
      data[[chid.name]] <- id.names
      lev.ch <- levels(data[[choice]])
      data <- data.frame(lapply(data, rep, length(lev.ch)))
      data[[alt.name]] <- rep(lev.ch, each = nb.id)
      row.names(data) <- paste(data[[chid.name]], data[[alt.name]], sep = ".")
    }
    data <- data[order(data[[chid.name]], data[[alt.name]]), ]
    # remplacer id par chid à gauche
    chid <- as.factor(data[[chid.name]])
    alt <- as.factor(data[[alt.name]])
    if (!is.null(alt.levels)){
      levels(data[[choice]]) <- alt.levels
      levels(alt) <- alt.levels
      row.names(data) <- paste(chid, alt, sep = ".")
    }
    data[[choice]] <- data[[choice]] == alt
  }
  chidpos <- which(names(data) == chid.name)
  altpos <- which(names(data) == alt.name)
  if (!is.null(id.var)){
    idpos <- which(names(data) == id.var)
    id.var <- as.factor(data[[id.var]])
  }
  if (drop.index){
    if (!is.null(id.var)) data <- data[, -c(chidpos, altpos, idpos)]
    else data <- data[, -c(chidpos, altpos)] 
  }
  
  if (!is.null(opposite)){
    for (i in opposite){
      data[[i]] <- -data[[i]]
    }
  }
  index <- data.frame(chid = chid, alt = alt)
  if (!is.null(id.var)) index <- cbind(index, id = id.var)
  rownames(index) <- rownames(data)
  attr(data, "index") <- index
  attr(data, "class") <- c("mlogit.data", "data.frame")
  data
}


"[.mlogit.data" <- function(x, i, j, drop = TRUE){
  mydata <- `[.data.frame`(x, i, j, drop = drop)
  index <- "[.data.frame"(attr(x, "index"), i,)
  index <- data.frame(lapply(index, function(x) x[drop = TRUE]), row.names = rownames(mydata))
  
  if (is.null(dim(mydata))){
    structure(mydata,
              index = index,
              class = c("pseries", class(mydata))
              )
  }
  else{
    structure(mydata,
              index = index,
              class = c("mlogit.data", "data.frame"))
  }
}

print.mlogit.data <- function(x, ...){
  attr(x, "index") <- NULL
  class(x) <- "data.frame"
  print(x, ...)
}

"[[.mlogit.data" <- function(x, y){
  index <- attr(x, "index")
  attr(x, "index") <- NULL
  class(x) <- "data.frame"
  result <- x[[y]]
  if (!is.null(result)){
    result <- structure(result,
                        class = c("pseries", class(x[[y]])),
                        index = index,
                        names = row.names(x)
                        )
  }
  result
}  

"$.mlogit.data" <- function(x,y){
  "[["(x, paste(as.name(y)))
}

print.pseries <- function(x, ...){
  attr(x, "index") <- NULL
  attr(x, "class") <- attr(x, "class")[-1]
  if (length(attr(x, "class")) == 1
      && class(x) %in% c("character", "logical", "numeric"))
    attr(x, "class") <- NULL
  print(x, ...)
}
