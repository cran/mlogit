#' Some deprecated functions, especially `mlogit.data`, `index` and
#' `mFormula`
#' 
#' `mlogit.data` is deprecated, use [dfidx::dfidx()] instead,
#' `mFormula` is replaced by [Formula::Formula()] and [zoo::index()]
#' by `idx`.
#' 
#' @name mlogit-deprecated
#' @aliases mlogit.data mFormula
#' @param x,object a `formula`, a `dfidx` or a `mlogit` object,
#' @param data a `data.frame`,
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
#' @param drop a boolean, equal to `FALSE` if one wants that a
#'     `data.frame` is always returned,
#' @param ... further arguments passed to `reshape`.
#' @return `mlogit.data` now returns a `dfidx` object, `mFormula`
#'     simply calls [Formula::Formula()] and returns a `Formula`
#'     object.
#' @author Yves Croissant
#' @seealso [stats::reshape()]
#' @keywords attribute
NULL

#' @rdname mlogit-deprecated
#' @export
mlogit.data <- function (data, choice = NULL, shape = c("long", "wide"), varying = NULL, 
    sep = ".", alt.var = NULL, chid.var = NULL, alt.levels = NULL, 
    id.var = NULL, group.var = NULL, opposite = NULL, drop.index = FALSE, 
    ranked = FALSE, subset = NULL, ...) {
    cldata <- match.call(expand.dots = TRUE)

    if ("subset" %in% names(cldata)){
        sub_call <- cldata
        m <- match(c("data", "subset"), names(sub_call), 0)
        sub_call <- sub_call[c(1, m)]
        names(sub_call)[2] <- "x"
        sub_call[[1]] <- as.name("subset")
        data <- eval(sub_call, parent.frame())
    }

    if (is.null(alt.var) & is.null(chid.var)) idx <- NULL
    if (is.null(alt.var) & ! is.null(chid.var)) idx <- list(chid.var)
    if (! is.null(alt.var) & is.null(chid.var)) idx <- list(NA, alt.var)                # cas compliqué
    if (! is.null(alt.var) & ! is.null(chid.var)) idx <- list(chid.var, alt.var)
    if (! is.null(group.var)) idx[[2]] <- c(idx[[2]], group.var)
    if (! is.null(id.var)){
        if (is.null(idx)) idx <- list(c(NA, id.var))
        else idx[[1]] <- c(idx[[1]], id.var)
    }
    cldata$idx <- idx
    cldata$pkg <- "mlogit"
    cldata$idnames <- c("chid", "alt")
    cldata$drop.index <- drop.index
    levpos <- match("alt.levels", names(cldata))
    if (! is.na(levpos)) names(cldata)[levpos] <- "levels"
    m <- match(c("data", "idx", "choice", "shape", "varying", "sep", "levels", "opposite", "idnames",
                 "pkg", "drop.index", "subset", "ranked"), names(cldata), 0)
    cldata <- cldata[c(1, m)]
    cldata$pkg <- "mlogit"
    cldata[[1]] <- as.name("dfidx")
    cldata[[1]] <- as.name("dfidx")

    # This works only if dfidx is attached, which is not the case in
    # some packages which imports (but not depends on) mlogit data <-
    # data <- eval(mldata, parent.frame())

    # Construct a list of dfidx' arguments with the default values            
    dfa <- list(data = NA, idx = NULL, drop.index = TRUE, as.factor = NULL, pkg = NULL,
                fancy.row.names = FALSE, subset = NULL,
                idnames = NULL, shape = "long", choice = NULL,
                varying = NULL, sep = ".", opposite = NULL, levels = NULL, ranked = FALSE)
    # Replace the relevant values by those indicated in mlogit
    dfa[names(cldata)[-1]] <- as.list(cldata)[- 1]
    # If mlogit.data is called from another function say choice
    # argument is a name in the call and should be replaced by its
    # value
    for (i in 2:length(dfa))
        if (is.name(dfa[[i]]))
            dfa[[i]] <- eval(dfa[[i]], parent.frame())
 
    # Run dfidx with these values

    data <- dfidx::dfidx(data = data, dfa$idx, drop.index = dfa$drop.index, as.factor = dfa$as.factor,
                         pkg = "mlogit", fancy.row.names = dfa$fancy.row.names,
                         idnames = dfa$idnames, shape = dfa$shape,
                         choice = dfa$choice, varying = dfa$varying,
                         sep = dfa$sep, opposite = dfa$opposite,
                  levels = dfa$levels, ranked = dfa$ranked)
    # add mlogit.data for backward compatibility with gmnl
    class(data) <- c(class(data), "mlogit.data")
    data
}

#' @rdname mlogit-deprecated
#' @export
mFormula <- function(object){
    UseMethod("mFormula")
}

#' @rdname mlogit-deprecated
#' @export
mFormula.formula <- function(object){
    if (! inherits(object, "Formula")) object <- Formula(object)
    if (! inherits(object, "mFormula")) class(object) <- c("mFormula", class(object))
    object
}

#' @rdname mlogit-deprecated
#' @export
mFormula.default <- function(object){
    stopifnot(inherits(object, "formula"))
    if (! inherits(object, "Formula")) object <- Formula(object)
    if (! inherits(object, "mFormula")) class(object) <- c("mFormula", class(object))
    object
}

#' @rdname mlogit-deprecated
#' @export
model.matrix.mFormula <- function(object, data, ...){
    if (inherits(data, "dfidx")){
        mf <- model.frame(formula = data, data = object, ...)
        model.matrix(mf, ...)
    }
    else{
        class(object) <- c("Formula", "formula")
        model.matrix(object, data, ...)
    }   
}


#' @rdname mlogit-deprecated
#' @export
is.mFormula <- function(object){
    inherits(object, "mFormula")
}


#' @importFrom zoo index
#' @export
zoo::index

#' @rdname mlogit-deprecated
#' @export
index.dfidx <- function(x, ...){
  idx(x, ...)
}

#' @rdname mlogit-deprecated
#' @export
index.mlogit <- function(x, ...){
  idx(x, ...)
}
