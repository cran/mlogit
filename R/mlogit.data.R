mlogit.data <- function(data, choice, shape = c("wide","long"), varying = NULL,
                        sep = ".", alt.var = NULL, chid.var = NULL, 
                        alt.levels = NULL, id = NULL, opposite = NULL, drop.index = FALSE, ...){
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
    if (!chid.is.variable) choiceid <- rep(1:n, each = J) else choiceid <- data[[chid.name]]

    if (!is.logical(data[[choice.name]])){
      if (is.factor(choice) && 'yes' %in% levels(choice))
        data[[choice.name]] <- data[[choice.name]] == 'yes'
      if (is.numeric(choice)) data[[choice.name]] <- data[[choice.name]] != 0
    }
    if (alt.is.variable){
      i <- which(names(data) == alt.name)
      if (drop.index) data <- data[-i]     
    }
    if (chid.is.variable){
      i <- which(names(data) == chid.name)
      if (drop.index) data <- data[-i]
    }
    # remplacer id par chid à gauche
    chid <- as.factor(choiceid)
    alt <- as.factor(alt)
    row.names(data) <- paste(choiceid, alt, sep = ".")
  }
  
  if (shape == "wide"){
    data[[choice]] <- as.factor(data[[choice]])
    if (is.null(alt.var)) alt <- "alt" else alt <- alt.var
    if (is.null(chid.var)) chid <- "chid" else chid <- chid.var
    if (!is.null(varying)){
      data <- reshape(data, varying = varying, direction = "long", sep = sep,
                      timevar = alt, idvar = chid,  ...)
    }
    else{
      id.names <- as.numeric(rownames(data))
      nb.id <- length(id.names)
      data[[chid]] <- id.names
      lev.ch <- levels(data[[choice]])
      data <- data.frame(lapply(data,rep,length(lev.ch)))
      data[[alt]] <- rep(lev.ch, each = nb.id)
      row.names(data) <- paste(data[[chid]], data[[alt]], sep = ".")
    }
    data <- data[order(data[[chid]], data[[alt]]), ]
    chidpos <- which(names(data) == chid)
    altpos <- which(names(data) == alt)
    # remplacer id par chid à gauche
    chid <- as.factor(data[[chid]])
    alt <- as.factor(data[[alt]])
    if (!is.null(id)){
      idpos <- which(names(data) == id)
      id <- as.factor(data[[id]])
      data <- data[, -c(chidpos, altpos, idpos)]
    }
    else  data <- data[, -c(chidpos, altpos)]
    if (!is.null(alt.levels)){
      levels(data[[choice]]) <- alt.levels
      levels(alt) <- alt.levels
      row.names(data) <- paste(chid, alt, sep = ".")
    }
    data[[choice]] <- data[[choice]] == alt
  }

  if (!is.null(opposite)){
    for (i in opposite){
      data[[i]] <- -data[[i]]
    }
  }
  index <- data.frame(chid = chid, alt = alt)
  if (!is.null(id)) index <- cbind(index, id = id)
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
