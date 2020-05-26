#' Indicates whether the formula contains an intercept
#' 
#' This is a generic which provide convenient methods for
#' formula/Formula object and for specific fitted models
#' 
#' @name has.intercept
#' @param object the object
#' @param rhs for the Formula method the rhs for which one wants to
#'     know if there is an intercept may be specified
#' @param ... further arguments
#' @export
#' @author Yves Croissant
#' @keywords attribute
has.intercept <- function(object, ...) {
    UseMethod("has.intercept")
}


#' @rdname has.intercept
#' @export
has.intercept.default <- function(object, ...) {
    has.intercept(formula(object), ...)
}

#' @rdname has.intercept
#' @export
has.intercept.formula <- function(object, ...) {
    attr(terms(object), "intercept") == 1L
}

#' @rdname has.intercept
#' @export
has.intercept.Formula <- function(object, rhs = NULL, ...) {
    ## NOTE: return a logical vector of the necessary length
    ## (which might be > 1)
    if(is.null(rhs)) rhs <- 1:length(attr(object, "rhs"))
    sapply(rhs, function(x) has.intercept(formula(object, lhs = 0, rhs = x)))
}

#' @rdname has.intercept
#' @export
has.intercept.mlogit <- function(object, ...){
    .formula <- object$formula
    ifelse(length(.formula)[2] == 1,
           has.intercept(.formula),
           has.intercept(.formula)[2])
}

#' Compute the model matrix for RUM
#' 
#' specific stuff compared to the model.matrix.dfidx method which
#' simply applies the Formula method
#' 
#' @name model.matrix.dfidx_mlogit
#' @param object the object
#' @param ...,lhs,rhs,dot see the `Formula` method
#' @export
#' @author Yves Croissant
#' @keywords attribute
model.matrix.dfidx_mlogit <- function(object, ..., lhs = NULL, rhs = 1, dot = "separate"){

    spliteffects <- TRUE
    
    if (is.null(attr(object, "formula"))) stop("the argument is an ordinary dfidx object")
    .formula <- attr(object, "formula")
    
    # check the presence of an intercept (in the second part or by
    # default in the first)

    .int <- ifelse(length(.formula)[2] == 1,
                   has.intercept(.formula),
                   has.intercept(.formula)[2])
    altname <- idx_name(object, 2)
    .altlevels <- levels(idx(object, 2))
    .reflevel <- .altlevels[1]
    toremove <- c()
    object <- unfold_idx(object)
    
    # For get consistent effets' names, always introduce the alt index
    # as a covariate.

    first_part <- formula(.formula, rhs = 1)
    first_part <- update(first_part, . ~ . - 1)
    first_part <- update(first_part, as.formula(paste(". ~ ", altname, " + .", sep = "")))

    # If there are no intercept, remove the J intercepts columns,
    # otherwise only the first (reflevel) one

    if (! .int) toremove <- c(toremove, paste(altname, .altlevels, sep = ""))
    else toremove <- c(toremove, paste(altname, .reflevel, sep = ""))
    if (length(.formula)[2] >= 2){
        second_part <- formula(.formula, rhs = 2, lhs = 0)
        second_part <- update(second_part, ~ . + 1)
        chid_spec_cov <- attr(terms(second_part), "term.labels")
        if (length(chid_spec_cov)){
            
            # get the effects and not the variable names in case of factors
            
            termlab <- colnames(model.matrix(second_part, object))[- 1]
            toremove <- c(toremove,
                          paste(paste(altname, .reflevel, sep = ""),
                                termlab, sep = ":")
                          )
            second_part <- update(second_part, as.formula(paste("~  (.): ", altname, sep = "")))
        }
        else second_part <- ~ 0
    }
    if (length(.formula)[2] >= 3){
        third_part <- formula(.formula, rhs = 3, lhs = 0)
        third_part <- update(third_part, as.formula(paste("~ (.): ", altname, sep = "")))
        # remove the intercept if any (bug with the vignette of the mnlogit package")
        if (has.intercept(third_part)) third_part <- update(third_part, . ~ . - 1)
    }

    # paste the 1/3 parts together to get a Formula object
    
    if (length(.formula)[2] == 1) .formula <- Formula(first_part)
    if (length(.formula)[2] == 2) .formula <- as.Formula(first_part, second_part)
    if (length(.formula)[2] == 3) .formula <- as.Formula(first_part, second_part, third_part)

    # get the model matrix and remove the relevant columns
    
    nrhs <- min(length(.formula)[2], 3)
    X <- model.matrix(.formula, object, ..., lhs = lhs, rhs = 1:nrhs, dot = dot)
    X <- X[, - c(match(toremove, colnames(X))), drop = FALSE]

    # customize the effects' names

    cnamesX <- colnames(X)
    
    # customize the intercepts analt:(Intercept)

    intnames <- paste(altname, .altlevels, sep = "")
    z <- match(intnames, cnamesX)
    cnamesX[na.omit(z)] <- paste(intnames[! is.na(z)], "(Intercept)", sep = ":")
#    cnamesX[na.omit(z)] <- paste("(Intercept)", .altlevels[! is.na(z)], sep = ":")

    # get a list of effects' names, split the interacted effects    
    # split the order of the two terms separated by the ':' operator
    # if the first element is part of intnames

    if (spliteffects){
        cnamesX <- strsplit(cnamesX, ":")
        cnamesX <- sapply(cnamesX, function(x)
            ifelse(length(x) > 1,
            ifelse(x[1] %in% intnames,
                   paste(x[2], x[1], sep = ":"),
                   paste(x[1], x[2], sep = ":")),
            x))
    }

    # replace for all effects' names namelevel by level
    
    for (i in seq_len(length(intnames))) cnamesX <- gsub(intnames[i], .altlevels[i], cnamesX)
    colnames(X) <- cnamesX
    X
}

