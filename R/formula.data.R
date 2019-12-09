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

has.intercept.mFormula <- function(object, ...){
    attr(object, "class") <- "Formula"
    has.int <- has.intercept(object)
    ifelse(length(has.int) > 1, has.int[2], has.int[1])
}

#' Model formula for logit models
#' 
#' Two kinds of variables are used in logit models: alternative specific and
#' individual specific variables. `mFormula` provides a relevant class to
#' deal with this specificity and suitable methods to extract the elements of
#' the model.
#' 
#' @name mFormula
#' @aliases mFormula is.mFormula mFormula.formula
#'     model.matrix.mFormula model.frame.mFormula
#' @import Formula
#' @param object for the `mFormula` function, a formula, for the
#'     `update` and `model.matrix` methods, a `mFormula` object,
#' @param formula a `mFormula` object,
#' @param data a `data.frame`,
#' @param lhs see `Formula`,
#' @param rhs see `Formula`,
#' @param alt.subset a vector of subset of alternatives one want to
#'     select,
#' @param reflevel the alternative selected to be the reference
#'     alternative,
#' @param ... further arguments.
#' @details Let `J` being the number of alternatives.  The formula may
#'     include alternative-specific and individual specific
#'     variables. For the latter, `J - 1` coefficients are estimated for
#'     each variable. For the former, only one (generic) coefficient
#'     or `J` different coefficient may be estimated.
#' 
#' A `mFormula` is a formula for which the right hand side may contain
#' three parts: the first one contains the alternative specific
#' variables with generic coefficient, *i.e.* a unique coefficient for
#' all the alternatives ; the second one contains the individual
#' specific variables for which one coefficient is estimated for all
#' the alternatives except one of them ; the third one contains the
#' alternative specific variables with alternative specific
#' coefficients.  The different parts are separeted by a `|` sign. If
#' a standard formula is writen, it is assumed that there are only
#' alternative specific variables with generic coefficients.
#' 
#' The intercept is necessarely alternative specific (a generic
#' intercept is not identified because only utility differences are
#' relevant). Therefore, it deals with the second part of the
#' formula. As it is usual in `R`, the default behaviour is to include
#' an intercept. A model without an intercept may be specified by
#' including `+ 0` or `- 1` in the second right-hand side part of the
#' formula. `+ 0` or `- 1` in the first and in the third part of the
#' formula are simply ignored.
#' 
#' Specific methods are provided to build correctly the model matrix
#' and to update the formula. The `mFormula` function is not intended
#' to be use directly. While using the [mlogit()] function, the first
#' argument is automaticaly coerced to a `mFormula` object.
#' @return an object of class `mFormula`.
#' @export
#' @author Yves Croissant
#' @keywords models
#' @examples
#' 
#' data("Fishing", package = "mlogit")
#' Fish <- mlogit.data(Fishing, varying = c(2:9), shape = "wide", choice =
#' "mode")
#' 
#' # a formula with to alternative specific variables (price and
#' # catch) and an intercept
#' f1 <- mFormula(mode ~ price + catch)
#' head(model.matrix(f1, Fish), 2)
#' 
#' # same, with an individual specific variable (income)
#' f2 <- mFormula(mode ~ price + catch | income)
#' head(model.matrix(f2, Fish), 2)
#' 
#' # same, without an intercept
#' f3 <- mFormula(mode ~ price + catch | income + 0)
#' head(model.matrix(f3, Fish), 2)
#' 
#' # same as f2, but now, coefficients of catch are alternative
#' # specific
#' f4 <- mFormula(mode ~ price | income | catch)
#' head(model.matrix(f4, Fish), 2)
#' 
mFormula <- function(object){
    UseMethod("mFormula")
}

#' @rdname mFormula
#' @method mFormula formula
#' @export
mFormula.formula <- function(object){
    if (!inherits(object, "Formula")) object <- Formula(object)
    class(object) <- c("mFormula", "Formula", "formula")
    object
}

mFormula.default <- function(object){
    stopifnot(inherits(object, "formula"))
    if (!inherits(object, "Formula")) object <- Formula(object)
    if (!inherits(object, "mFormula"))
        class(object) <- c("mFormula", class(object))
    object
}

#' @rdname mFormula
#' @export
is.mFormula <- function(object){
    inherits(object, "mFormula")
}

#' @rdname mFormula
#' @method model.frame mFormula
#' @export
model.frame.mFormula <- function(formula, data, ..., lhs = NULL, rhs = NULL, alt.subset = NULL, reflevel = NULL){
    if (is.null(rhs)) rhs <- 1:(length(formula)[2])
    if (is.null(lhs)) lhs <- ifelse(length(formula)[1] > 0, 1, 0)
    index <- attr(data, "index")
    mf <- model.frame(as.Formula(formula), as.data.frame(data), ..., rhs = rhs)
    if (inherits(data, "mlogit.data")){
        # obvious condition, model.frame.mFormula is intended for
        # mlogit.data and not ordinary data.frame, but the mpbart
        # package use it in this way
        index <- index[rownames(mf), ]
        oindex <- index
        # a/ coerce the response to a logical if necessary
        y <- model.response(mf)
        if (! is.logical(y)){
            if (is.factor(y)){
                if (length(levels(y)) != 2)
                    stop("the number of levels for the choice variable should equal two")
                y <- y == levels(y)[2]
            }
            if (is.numeric(y)) y <- y != 0
        }
        mf[[1]] <- y
        # b/ change the reference level of the response if required
        if (! is.null(reflevel)){
            attr(index, "alt.ordering") <- levels(index$alt)
            index$alt <- relevel(index$alt, reflevel)
        }
        # c/ compute the relevent subset if required
        if (! is.null(alt.subset)){
            # we keep only choices that belong to the subset
            choice <- index$alt[model.response(mf)]
            choice <- choice %in% alt.subset
            unid <- unique(index$chid)
            names(choice) <- as.character(unid)
            id.kept <- choice[as.character(index$chid)]
            # we keep only the relevant alternatives
            alt.kept <- index$alt %in% alt.subset
            # the relevant subset for the data.frame and the indexes
            mf <- mf[id.kept & alt.kept, , drop = FALSE]
            index <- index[id.kept & alt.kept, , drop = FALSE]
            # drop the non-selected levels
            index$alt <- index$alt[drop = TRUE]
            if (! is.null(attr(index, "alt.ordering"))){
                attr(index, "alt.ordering") <-
                    intersect(attr(index, "alt.ordering"), levels(index$alt))
            }
        }
        # d/ balance the data.frame i.e. insert rows with NA when an
        # alternative is not relevant
#        alt.un <- unique(index$alt)
        alt.un <- levels(index$alt)
        # be very careful, if some observations are missing unique may
        # be different from levels
        chid.un <- unique(index$chid)
        alt.lev <- levels(index$alt)
        J <- length(alt.lev)
        n <- length(chid.un)
        T <- length(alt.un)
        if (nrow(mf) != (n * T)){
            rownames(mf) <- paste(index$chid, index$alt, sep = ".")
            all.rn <- as.character(t(outer(chid.un, alt.un, paste, sep = ".")))
            mf <- mf[all.rn, ]
            rownames(mf) <- all.rn
            chid <- rep(chid.un, each = T)
            alt <- factor(rep(alt.un, n), levels = alt.un, labels = alt.un)
            # take care here to make a factor with the levels in the
            # correct order
            alt[is.na(mf[[1]])] <- NA
            index <- data.frame(chid = chid, alt = alt, row.names = rownames(mf))
            if (! is.null(oindex$group)){
                ra <- oindex[c("alt", "group")][! duplicated(oindex$alt), ]
                gps <- ra$group
                names(gps) <- ra$alt
                index$group <- gps[index$alt]
            }
            if (! is.null(oindex$id)){
                ra <- oindex[c("chid", "id")][! duplicated(oindex$chid), ]
                ids <- ra$id
                names(ids) <- ra$chid
# YC 2019/12/09 : factor indexation leads to a position and not a character indexing                
#                index$id <- ids[index$chid]
                index$id <- ids[as.character(index$chid)]
                index$id <- index$id[drop = TRUE]
            }
        }
    }
#    index <- data.frame(lapply(index, function(x) x[drop = TRUE]), row.names = rownames(index))
    index$chid <- index$chid[drop = TRUE]
    structure(mf,
              index = index,
              formula = formula,
              choice = names(mf)[1],
              class = c("mlogit.data", class(mf)))
}

#' @rdname mFormula
#' @method model.matrix mFormula
#' @export
model.matrix.mFormula <- function(object, data, ...){
    if (length(object)[2] == 4) object <- formula(object, rhs = 1:3)
    
    K <- length(data)
    omitlines <- attr(na.omit(data), "na.action")
    index <- attr(data, "index")
    alt <- index[["alt"]]
    chid <- index[["chid"]]
    data$alt <- alt
    resp.name <- as.character(attr(object, "lhs")[[1]])
    # keep track of the existence of an intercept
    has.int <- has.intercept(object)
    if (has.int) intercept.char <- "alt" else intercept.char <- NULL
  
    ## for ind.spec : remove any 0 or 1 or -1 in the formula and get the
    ## list of the variables
    if (length(object)[2] > 1){
        ind.spec <- formula(object, rhs = 2, lhs = 0)
        if (!has.int) ind.spec <- update(ind.spec, ~ . + 1)
        ind.spec <- update(ind.spec, ~ .)
        ind.spec.char <- as.character(ind.spec)[2]
        if (ind.spec.char == "1") ind.spec.char <- ind.spec.var <- NULL
        else{
            # ind.spec.var <- attr(terms(ind.spec), "term.labels") the
            # following lines extract the effects and not the variable
            # names, useful for factors
            ind.spec.var <- colnames(model.matrix(update(ind.spec, ~.+1), data))[-1]
            ind.spec.char <- paste("(", ind.spec.char, "):alt", sep="")
        }
    }
    else ind.spec <- ind.spec.char <- ind.spec.var <- NULL

    # alternative specific variables
    alt.spec <- formula(object, rhs = 1, lhs = 0)
    alt.spec <- update(update(alt.spec, ~ . + 1), ~ .)
    alt.spec.char <- as.character(alt.spec)[2]
    if (alt.spec.char == "1") als.spec <- alt.spec.char <- NULL

    # specific coefficient for alternative specific variables
    if (length(object)[2] == 3){
        coef.spec <- formula(object, rhs = 3, lhs = 0)
        coef.spec <- update(update(coef.spec, ~ . + 1), ~ .)
        coef.spec.char <- as.character(coef.spec)[2]
        if (!is.null(coef.spec.char)) coef.spec.char <- paste("(", coef.spec.char, "):alt", sep="")
    }
    else coef.spec <- coef.spec.char <- NULL
    form.char <- paste(c(intercept.char, alt.spec.char,
                         ind.spec.char, coef.spec.char),
                       collapse = "+")
    formula <- as.formula(paste(resp.name, " ~ ", form.char))
    X <- model.matrix(formula, data)[, -1, drop = F]
    lev1 <- levels(alt)[1]
    lev1 <- paste("alt", lev1, sep = "")
    toremove <- unlist(lapply(as.list(ind.spec.var), function(x) paste(lev1, x, sep = ":")))
    revtoremove <- unlist(lapply(as.list(ind.spec.var), function(x) paste(x, lev1, sep = ":")))
    toremove <- colnames(X) %in% c(toremove, revtoremove)
    X <- X[, ! toremove, drop = FALSE]
    # the following lines suppress the mentions to 'alt' in the names of
    # the effects and add a mention to '(intercept)'
    namesX <- colnames(X)
    for (i in 1:length(namesX)) namesX[i] <- sub('alt', '', namesX[i])
    z <- match(levels(alt), namesX)
    namesX[na.omit(z)] <- paste(levels(alt)[!is.na(z)], '(intercept)', sep=":")
    colnames(X) <- namesX
    qrX <- qr(na.omit(X))
    # remove the linear dependant columns
    # X <- X[, qrX$pivot[1:qrX$rank], drop = FALSE]
    attr(X, "index") <- index
    X
}

as.Formula.mFormula <- function(x, ...){
    class(x) <- class(x)[-1]
    x
}

#' data.frame for logit model
#' 
#' shape a `data.frame` in a suitable form for the use of the
#' `mlogit` function.
#' 
#' @name mlogit.data
#' @aliases mlogit.data [[<-.mlogit.data $<-.mlogit.data print.pseries
#'     index mean.mlogit.data formula.mlogit.data print.mlogit.data
#' @param data a `data.frame`,
#' @param x,object a `mlogit.data` or a `pseries` object,
#' @param choice the variable indicating the choice made: it can be
#'     either a logical vector, a numerical vector with 0 where the
#'     alternative is not chosen, a factor with level 'yes' when the
#'     alternative is chosen
#' @param shape the shape of the `data.frame`: whether `long` if each
#'     row is an alternative or `wide` if each row is an observation,
#' @param varying the indexes of the variables that are alternative
#'     specific,
#' @param sep the seperator of the variable name and the alternative
#'     name (only relevant for a `wide` `data.frame`),
#' @param alt.var the name of the variable that contains the
#'     alternative index (for a `long` `data.frame` only) or the name
#'     under which the alternative index will be stored (the default
#'     name is `alt`),
#' @param chid.var the name of the variable that contains the choice
#'     index or the name under which the choice index will be stored,
#' @param alt.levels the name of the alternatives: if null, for a
#'     `wide` data.frame, they are guessed from the variable names and
#'     the choice variable (both should be the same), for a `long`
#'     `data.frame`, they are guessed from the `alt.var` argument,
#' @param id.var the name of the variable that contains the individual
#'     index if any,
#' @param group.var the name of the variable that contains the group
#'     index if any,
#' @param opposite returns the opposite of the specified variables,
#' @param drop.index should the index variables be dropped from the
#'     `data.frame`,
#' @param ranked a logical value which is true if the response is a
#'     rank,
#' @param subset a logical expression which defines the subset of
#'     observations to be selected,
#' @param i the rows to extract,
#' @param j the columns to extract,
#' @param y the column of the `data.frame` to extract or to replace,
#' @param value the replacement value,
#' @param drop a boolean, equal to `FALSE` if one wants that a
#'     `data.frame` is always returned,
#' @param ... further arguments passed to `reshape`.
#' @return A `mlogit.data` object, which is a `data.frame` in `long`
#'     format, *i.e.* one line for each alternative.  It has a `index`
#'     attribute, which is a `data.frame` that contains the index of
#'     the choice made (`chid`), the index of the alternative (`alt`)
#'     and, if any, the index of the individual (`id`) and of the
#'     alternative groups (`group`).  The choice variable is a boolean
#'     which indicates the choice made. This function use
#'     [stats::reshape()] if the `data.frame` is in `wide` format.
#' @export
#' @author Yves Croissant
#' @seealso [stats::reshape()]
#' @keywords attribute
#' @examples
#' # ModeChoice is a long data.frame 
#' 
#' data("TravelMode", package = "AER")
#' TM <- mlogit.data(TravelMode, choice = "choice", shape = "long",
#'                  alt.levels = c("air", "train", "bus", "car"))
#' 
#' # Same but the alt variable called mode is provided
#' 
#' TM <- mlogit.data(TravelMode ,choice = "choice", shape = "long",
#'                   alt.var = "mode")
#' 
#' # Same but the chid variable called individual is provided
#' 
#' TM <- mlogit.data(TravelMode, choice = "choice",
#'                   shape = "long", id.var = "individual",
#'                   alt.levels = c("air", "train", "bus", "car"))
#' 
#' # Same but with two own provided variables
#' 
#' TM <- mlogit.data(TravelMode, choice = "choice", shape = "long",
#'                  id.var = "individual", alt.var = "mode")
#' 
#' #  Same but with two own provided variables which are deleted from the
#' #  data.frame
#' 
#' TM <- mlogit.data(TravelMode, choice = "choice", shape = "long",
#'                  id.var = "individual", alt.var = "mode", drop.index = TRUE)
#' 
#' #  Train is a wide data.frame with columns 'choiceid' is the choice
#' #  index, the alternatives are named "ch1" and "ch2", the opposite of
#' #  the variables is returned
#' 
#' data("Train", package = "mlogit")
#' Train <- mlogit.data(Train, choice = "choice", shape = "wide",
#'                      varying = 4:11, alt.levels = c("A", "B"), sep = "_",
#'                      opposite = c("price", "time", "change", "comfort"))
#' 
#' data("HC", package = "mlogit")
#' HC <- mlogit.data(HC, choice = "depvar", varying=c(2:8, 10:16), shape="wide")
#' 
#' # Game is a data.frame in wide format for which the response is a
#' #  ranking variable
#' 
#' data("Game", package = "mlogit")
#' G <- mlogit.data(Game, shape="wide", varying = 1:12, alt.var = 'platform',
#'                  drop.index = TRUE, choice="ch", ranked =TRUE)
#' 
#' # Game2 contains the same data, but in long format 
#' data("Game2", package = "mlogit")
#' G2 <- mlogit.data(Game2,  shape='long', choice="ch", alt.var = 'platform', ranked = TRUE)
mlogit.data <- function(data, choice = NULL, shape = c("long", "wide"), varying = NULL,
                        sep = ".", alt.var = NULL, chid.var = NULL, 
                        alt.levels = NULL, id.var = NULL, group.var = NULL,
                        opposite = NULL, drop.index = FALSE, ranked = FALSE,
                        subset = NULL, ...){
    # chid.name, alt.name : the name of the index variables
    # chid, alt : the index variables

    # if a subset argument is prodided, subset the original data frame
    cldata <- match.call(expand.dots = TRUE)
    if (match("subset", names(cldata), 0)){
        m <- match(c("data", "subset"), names(cldata), 0)
        cldata <- cldata[c(1, m)]
        names(cldata)[2] <- "x"
        cldata[[1]] <- as.name("subset")
        data <- eval(cldata, parent.frame())
    }
    
    # the long shape is now the default
    shape <- match.arg(shape)
    
    if (shape == "long"){
        if (is.null(chid.var)){
            chid.name <- "chid"
            chid.is.variable <- FALSE
        }
        else{
            chid.name <- chid.var
            chid.is.variable <- ifelse(is.null(data[[chid.var]]), FALSE, TRUE)
        }
        if (is.null(alt.var) && is.null(alt.levels))
            stop("at least one of alt.var and alt.levels should be filled")
        
        if (! is.null(alt.levels)){
            J <- length(alt.levels)
            n <- nrow(data) / J
            alt <- factor(rep(alt.levels, n), levels = alt.levels)
            if (!is.null(alt.var) && !is.null(data[[alt.var]])){
                warning(paste("variable", alt.var, "exists and will be replaced"))
                alt.is.variable <- TRUE
            }
            else alt.is.variable <- FALSE
            alt.name <- ifelse(is.null(alt.var), "alt", alt.var)
        }
        else{
            alt.name <- alt.var
            alt.is.variable <- TRUE
            if (!is.factor(data[[alt.name]])) data[[alt.name]] <- factor(data[[alt.name]])
            alt.levels <- levels(data[[alt.name]])
            J <- length(alt.levels)
            alt <- data[[alt.name]]

        }
        n <- nrow(data) / J
        if (! chid.is.variable) chid <- rep(1:n, each = J) else chid <- data[[chid.name]]
        if (! ranked & ! is.null(choice)){
            choice.name <- choice
            choice <- data[[choice]]
            # coerce the choice variable to a logical if necessary
            if (! is.logical(choice)){
                if (is.factor(choice)){
                    if (length(levels(choice)) != 2)
                        stop("the number of levels for the choice variable should equal two")
                    choice <- choice == levels(choice)[2]
                }
                if (is.numeric(choice)) choice <- choice != 0
            }
            data[[choice.name]] <- choice
        }
        chid <- as.factor(chid)
        alt <- as.factor(alt)
        row.names(data) <- paste(chid, alt, sep = ".")
    }
    
    if (shape == "wide"){
        if (! ranked){
            if (is.null(choice)) stop("the choice argument is mandatory for wide-shaped data.frame")
            choice.name <- choice
            if (is.ordered(data[[choice]])) class(data[[choice]]) <- "factor"
            else data[[choice]] <- as.factor(data[[choice]])
        }
        # this doesn't work for ordered factors which remains ordered
        if (is.null(alt.var)) alt.name <- "alt" else alt.name <- alt.var
        if (is.null(chid.var)) chid.name <- "chid" else chid.name <- chid.var
        if (! is.null(varying)){
            data <- reshape(data, varying = varying, direction = "long", sep = sep,
                            timevar = alt.name, idvar = chid.name,  ...)
        }
        else{
            if (ranked)
                stop("for ranked data in wide format, the varying argument is mandatory and should contain at least the choice variable")
            id.names <- as.numeric(rownames(data))
            nb.id <- length(id.names)
            data[[chid.name]] <- id.names
            lev.ch <- levels(data[[choice]])
            data <- data.frame(lapply(data, rep, length(lev.ch)))
            data[[alt.name]] <- rep(lev.ch, each = nb.id)
            row.names(data) <- paste(data[[chid.name]], data[[alt.name]], sep = ".")
        }
        data <- data[order(data[[chid.name]], data[[alt.name]]), ]
        chid <- as.factor(data[[chid.name]])
        alt <- as.factor(data[[alt.name]])
        if (! is.null(alt.levels)){
            levels(data[[choice]]) <- alt.levels
            levels(alt) <- alt.levels
            row.names(data) <- paste(chid, alt, sep = ".")
        }       
        if (! ranked){
            #YC 20180106 never chosen alternatives are added to the levels of the choice variable
            data[[choice]] <- factor(data[[choice]], levels = levels(alt))
            data[[choice]] <- data[[choice]] == alt
        }
        else{
            if (is.null(data[[choice]])) stop("the choice variable doesn't exist")
        }
    }

    chidpos <- which(names(data) == chid.name)
    altpos <- which(names(data) == alt.name)
    if (! is.null(id.var)){
        idpos <- which(names(data) == id.var)
        id.var <- as.factor(data[[id.var]])
    }
    if (! is.null(group.var)){
        grouppos <- which(names(data) == group.var)
        group.var <- as.factor(data[[group.var]])
    }
    
    if (drop.index){
        if (! is.null(id.var)) data <- data[, -c(chidpos, altpos, idpos)]
        else data <- data[, -c(chidpos, altpos)] 
    }
    
    if (!is.null(opposite)){
        for (i in opposite){
            data[[i]] <- - data[[i]]
        }
    }
    index <- data.frame(chid = chid, alt = alt)
    if (! is.null(id.var)) index <- cbind(index, id = id.var)
    if (! is.null(group.var)) index <- cbind(index, group = group.var)
    rownames(index) <- rownames(data)
    attr(data, "index") <- index
    attr(data, "class") <- c("mlogit.data", "data.frame")
    if (ranked) data <- mlogit2rank(data, choicename = choice)
    if (! ranked & ! is.null(choice)) attr(data, "choice") <- choice.name
    data
}


formula.mlogit.data <- function(x, ...) attr(x, "formula")

#' @rdname mlogit.data
#' @export
print.mlogit.data <- function(x, ...){
  attr(x, "index") <- NULL
  class(x) <- "data.frame"
  print(x, ...)
}

#' @rdname mlogit.data
#' @export
index.mlogit.data <- function(x, ...){
  attr(x, "index")
}

#' @rdname mlogit.data
#' @export
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

#' @rdname mlogit.data
#' @export
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

#' @rdname mlogit.data
#' @export
"$.mlogit.data" <- function(x,y){
  "[["(x, paste(as.name(y)))
}

#' @rdname mlogit.data
#' @export
"$<-.mlogit.data" <- function(object, y, value){
  # object : le data.frame
  # y : la variable
  # value : la nouvelle valeur
  object[[y]] <- value
  object
}

#' @rdname mlogit.data
#' @export
"[[<-.mlogit.data" <- function(object, y, value){
  if (class(value)[1] == "pseries"){
    class(value) <- class(value)[-1]
    attr(value, "index") <- NULL
  }
  object <- "[[<-.data.frame"(object, y, value = value)
  object
}

#' @rdname mlogit.data
#' @method mean mlogit.data
#' @export
mean.mlogit.data <- function(x, ...){
    alt <- index(x)$alt
    J <- length(levels(alt))
    result <- data.frame(lapply(x,
                                function(x){
                                    if (is.numeric(x)) result <- as.numeric(tapply(x, alt, mean))
                                    else{
                                        if (is.logical(x)){
                                            z <- tapply(x, alt, sum)
                                            result <- z == max(z)
                                        }
                                        if(is.character(x)){
                                            x <- factor(x, levels = unique(x))
                                        }
                                        if (is.factor(x)){
                                            result <- factor(names(which.max(table(x))), levels = levels(x))
                                        }
                                    }
                                    result
                                }
                                )
                         )
    attr(result, "index") <- data.frame(alt = factor(levels(alt), levels =  levels(alt)), chid = rep(1, J))
    rownames(result) <- rownames(attr(result, "index")) <- paste(rep(1, J), levels(alt), sep = ".")
    class(result) <- c("mlogit.data", "data.frame")
    result
}

mlogit2rank <- function(x, choicename, ...){
  choicepos <- match(choicename, names(x))
  id <- attr(x, "index")$chid
  lev.id <- levels(id)
  theid <- as.numeric(id)
  oalt <-  attr(x, "index")$alt
  lev.alt <- levels(oalt)
  choice <- x[[choicename]]
  J <- length(unique(choice))
  d <- data.frame()
  chid <- c()
  alt <- c()
  id <- c()
  k <- 0
  for (i in unique(theid)){
    aid <- which(theid == i)
    adata <- x[aid, - choicepos]
    achoice <- choice[aid]
    aalt <- oalt[aid]
    remAlts <- rep(TRUE, J)
    alogchoice <- achoice == 1
    d <- rbind(d, cbind(adata, alogchoice))
    Z <- sum(remAlts)
    k <- k + 1
    chid <- c(chid, rep(k, Z))
    id <- c(id, rep(i, Z))
    alt <- c(alt, aalt)
    for (j in 1:(J - 2)){
      k <- k + 1
      min.index <- achoice == j
      remAlts[min.index] <- FALSE
      Z <- sum(remAlts)
      chid <- c(chid, rep(k, Z))
      alt <- c(alt, aalt[remAlts])
      id <- c(id, rep(i, Z))
      alogchoice <- achoice[remAlts] == j + 1
      d <- rbind(d, cbind(adata[remAlts,], alogchoice))
    }
  }
  colnames(d)[length(d)] <- choicename
  alt <- factor(alt, labels = lev.alt)
  index <- data.frame(chid = chid, alt = alt, id = id)
  rownames(d) <- rownames(index) <- paste(chid, as.character(alt), sep = ".")
  structure(d, index = index, class = c('mlogit.data', 'data.frame'))
}
