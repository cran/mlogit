lnl.mlogits <- function(param, X, y, weights = NULL, gradient = FALSE,
                        hessian = FALSE, opposite = TRUE, sumlnl = TRUE,
                        direction = 0, initial.value = NULL, steptol = 1E-01){
  opposite <- ifelse(opposite, -1, +1)
  balanced <- FALSE
  if (is.null(weights)) weights <- 1
  step <- 2
  repeat{
    step <- step / 2
    if (step < steptol) break
    eXb <- lapply(X, function(x) exp(crossprod(t(x), param + step * direction)))
#    seXb <- Reduce("+", eXb)
    # Reduce doesn't manage correctly the NAs, so we switch back to suml
    seXb <- suml(eXb)
    P <- lapply(eXb, function(x){v <- x/seXb; v[is.na(v)] <- 0; as.vector(v)})
    Pch <- Reduce("+", mapply("*", P, y, SIMPLIFY = FALSE))
    lnl <- sum(opposite * weights * log(Pch))
    if (is.null(initial.value) || lnl <= initial.value) break
  }
  if (gradient | hessian) PX <- suml(mapply("*", X, P, SIMPLIFY = FALSE))
  if (gradient){
    Xch <- suml(mapply("*", X, y, SIMPLIFY = FALSE))
    gradi <-  opposite * weights * (Xch - PX)
    attr(lnl, "gradi") <- gradi
    attr(lnl, "gradient") <- if (is.matrix(gradi)) apply(gradi,2,sum) else sum(gradi)
  }
  if (hessian){
    XmPX <- lapply(X, function(x){g <- x - PX; g[is.na(g)] <- 0; g})
    hessian <-   - suml( mapply(function(x, y) crossprod(x*y, y),
                                P, XmPX, SIMPLIFY = FALSE))
    attr(lnl, "hessian") <- opposite * hessian
  }
  if (step < steptol) lnl <- NULL
  else{
    attr(lnl, "probabilities") <- P
    attr(lnl, "step") <- step
  }
  lnl
}

lnl.rlogit <- function(param, y, Xa, Xc,
                       weights = NULL, gradient = FALSE, hessian = FALSE,
                       opposite = TRUE, sumlnl = TRUE,
                       direction = rep(0, length(param)), initial.value = NULL,
                       steptol = 1E-1,
                       Varc, Vara, random.nb,
                       id, rpar, correlation){

  if (is.null(weights)) weights <- rep(1, length(y[[1]]))
  opposite <- ifelse(opposite, -1, +1)
  step <- 2
  repeat{
    step <- step / 2
    if (step < steptol) break
    betac <- param[Varc] + step * direction[Varc]
    mua <- param[Vara] + step * direction[Vara]
    Kc <- length(Varc)
    Ka <- length(Vara)
    n <- length(y[[1]])
    K <- Kc + Ka
    R <- nrow(random.nb)
    siga <- param[-c(1:K)] + step * direction[-c(1:K)]
    # seems redondant for uncorrelated models and false for correlated ones
    if (!correlation)
      names(mua) <- names(siga) <- colnames(Xa[[1]])
    b <- make.beta(mua, siga, rpar, random.nb, correlation)
    betaa <- b$betaa
    A <- lapply(Xc, function(x) as.vector(crossprod(t(as.matrix(x)), betac)))
    B <- lapply(Xa, function(x) tcrossprod(as.matrix(x), betaa))
    AB <- mapply(function(x, y) exp(x + y), A, B, SIMPLIFY=FALSE)
    S <- suml(AB)
    P <- lapply(AB, function(x) x/S)
    probabilities <- as.matrix(data.frame(lapply(P,
                                                 function(x) apply(x, 1, mean))))
    Pch <- suml(mapply("*", P, y, SIMPLIFY=FALSE))
    if (!is.null(id)){
      Pch <- apply(Pch, 2, tapply, id, prod)
    }
    pm <- apply(Pch, 1, mean)
    if (!is.null(id)) lnl <- opposite * sum(weights[!duplicated(id)]*log(pm))
    else lnl <- opposite * sum(weights*log(pm))
    if (is.null(initial.value) || lnl <= initial.value) break
  }
  if (gradient){
    xac <- suml(mapply("*", Xa, y, SIMPLIFY=FALSE))
    xcc <- suml(mapply("*", Xc, y, SIMPLIFY=FALSE))
    if (!is.null(id)) Pch <- Pch[as.character(id), ]
    if (correlation){
      names.cor <- c()
      for (i in 1:Ka){
        names.cor <- c(names.cor, paste(names(rpar)[i], names(rpar)[i:Ka], sep=":"))
      }
      vecX <- c()
      for (i in 1:Ka){
        vecX <- c(vecX, i:Ka)
      }
      Xas <- lapply(Xa,  function(x) x[, vecX])
      xac <- suml(mapply("*", Xa, y, SIMPLIFY=FALSE))
      xacs <- suml(mapply("*", Xas, y, SIMPLIFY=FALSE))
      colnames(xacs) <- names(param)[-(1:(Ka+Kc))]
    }
    else{
      xacs <- xac
      Xas <- Xa
    }
    if (!is.null(id)){
      Pch <- Pch[as.character(id), ]
      pm <- apply(Pch, 1, mean)
    }

    PCP <- lapply(P, function(x) Pch*x)
    PCPs <- lapply(PCP, function(x) apply(x, 1, sum))
    grad.cst <- xcc-suml(mapply("*", Xc, PCPs, SIMPLIFY=FALSE))/(R*pm)
    grad.mu <- (tcrossprod(Pch, t(b$betaa.mu))*xac -
                suml(mapply(function(x, y) x*tcrossprod(y, t(b$betaa.mu)),
                            Xa, PCP, SIMPLIFY=FALSE)))/(R*pm)
    grad.sd <- (tcrossprod(Pch, t(b$betaa.sigma))*xacs -
                suml(mapply(function(x, y) x*tcrossprod(y, t(b$betaa.sigma)),
                            Xas, PCP, SIMPLIFY=FALSE)))/(R*pm)
    gradi <- matrix(NA, n, K + ncol(grad.sd))
    gradi[, Varc] <- - grad.cst
    gradi[, Vara] <- - grad.mu
    gradi[, -c(1:K)] <- - grad.sd
    if (!is.null(weights)) gradi <- weights * gradi
    colnames(gradi) <- names(param)
    
    attr(lnl, "gradi") <-  opposite * gradi
    attr(lnl, "gradient") <- apply(gradi, 2, sum)
    attr(lnl, "step") <- step
  }
  if (step < steptol) lnl <- NULL
  else{
    attr(lnl, "probabilities") <- probabilities
    attr(lnl, "step") <- step
  }
  lnl
}
    
lnl.nlogit <- function(param, X, y, weights = NULL, gradient = FALSE, hessian = FALSE,
                       opposite = TRUE, sumlnl = TRUE, nests, un.nest.el = FALSE,
                       unscaled = FALSE, direction = rep(0, length(param)),
                       initial.value = NULL, steptol = 1E-01){
  if (un.nest.el){
    lambda <- param[length(param)]
    param <- c(param, rep(lambda, length(nests)-1))
  }
  thealts <- names(X)
  posalts <- lapply(thealts, function(x) which(unlist(nests) %in% x))
  opposite <- ifelse(opposite, -1, +1)
  if (is.null(weights)) weights <- 1
  K <- ncol(X[[1]])
  n <- nrow(X[[1]])
  J <- length(nests)
  X <- lapply(nests, function(x) X[x])
  Y <- lapply(nests, function(x) y[x])
  Yn <- lapply(Y, function(x) Reduce("+", x))
  step <- 2
  repeat{
    step <- step / 2
    if (step < steptol) break
    beta <- param[1:K] + step * direction[1:K]
    lambda <- param[-c(1:K)] + step * direction[-c(1:K)]
    names(lambda) <- names(nests)
    V <- lapply(X, function(x)
                lapply(x, function(y)
                       as.numeric(crossprod(t(y), beta))
                       )
                )
    if (!unscaled)
      W <- mapply(function(v, l) lapply(v, function(x) x / l), V, lambda, SIMPLIFY = FALSE)
    else W <- V
    A <- lapply(W, function(x) lapply(x, exp))
    N <- lapply(A, function(x) Reduce("+", x))
    Pjl <- mapply(function(a, n) lapply(a, function(x) x/n), A, N, SIMPLIFY = FALSE)
    Pl <- mapply(function(n, l) n^l, N, lambda, SIMPLIFY=FALSE)
    D <- Reduce("+", Pl)
    Pl <- lapply(Pl, function(x) x/D)
    P <- mapply(function(pjl, pl) lapply(pjl, function(x) x * pl), Pjl, Pl, SIMPLIFY = FALSE)
    L <- mapply(function(p, y) mapply("*", p, y, SIMPLIFY = FALSE), P, Y, SIMPLIFY = FALSE)
    L <- Reduce("+", lapply(L, function(x) Reduce("+", x)))
    lnl <- sum(opposite * weights * log(L))
    if (is.null(initial.value) || lnl <= initial.value) break
  }
  if (gradient){
    ### For overlaping nests
    Pj <- unlist(P, recursive = FALSE)
    Pj <- lapply(posalts, function(x) Reduce("+", Pj[x]))
    names(Pj) <- thealts
    proba <- Pj
    Pj <- lapply(nests, function(x) Pj[x])
    Pond <- mapply(function(p, pj) mapply("/", p, pj, SIMPLIFY = FALSE), P, Pj, SIMPLIFY = FALSE)
    ###
    Xb <- mapply(function(x, pjl)
                 Reduce("+",mapply("*", x, pjl, SIMPLIFY = FALSE)),
                 X, Pjl, SIMPLIFY = FALSE)
    Vb <- mapply(function(v, pjl)
                 Reduce("+",mapply("*", v, pjl, SIMPLIFY = FALSE)),
                 V, Pjl, SIMPLIFY = FALSE)
    if (!unscaled)
      Xbtot <- Reduce("+", mapply("*", Pl, Xb, SIMPLIFY = FALSE))
    else
      Xbtot <- Reduce("+", mapply(function(pl, xb, l) pl * xb * l,
                                  Pl, Xb, lambda, SIMPLIFY = FALSE))

    if (!unscaled)
      Gb <- mapply(function(x, xb, l)
                   lapply(x, function(z) (z+(l-1)*xb)/l), X, Xb, lambda, SIMPLIFY = FALSE)
    else
      Gb <- mapply(function(x, xb, l)
                   lapply(x, function(z) (z+(l-1)*xb)), X, Xb, lambda, SIMPLIFY = FALSE)
    
    Gb <-  mapply(function(gb, y)
                  mapply("*", gb, y, SIMPLIFY = FALSE),
                  Gb, Y, SIMPLIFY = FALSE)
    Gb <-  mapply(function(gb, pond)
                  mapply("*", gb, pond, SIMPLIFY = FALSE),
                  Gb, Pond, SIMPLIFY = FALSE)
    Gb <- Reduce("+", lapply(Gb, function(x) Reduce("+", x))) - Xbtot

    if (!unscaled){
      Gl1 <- mapply(function(v, n, vb, l)
                    lapply(v, function(z) -(z - l^2 * log(n) + (l-1) * vb)/l^2),
                    V, N, Vb, lambda, SIMPLIFY = FALSE)
      Gl1 <- mapply(function(gl1, y) mapply("*", gl1, y, SIMPLIFY = FALSE),
                    Gl1, Y, SIMPLIFY = FALSE)
      Gl1 <- mapply(function(gl1, pond)
                    mapply("*", gl1, pond, SIMPLIFY = FALSE),
                    Gl1, Pond, SIMPLIFY = FALSE)
      Gl1 <- lapply(Gl1, function(x) Reduce("+", x))
      Gl2 <- mapply(function(vb, n, l, pl) - pl * (l^2 * log(n)- l * vb)/ l^2,
                    Vb, N, lambda, Pl, SIMPLIFY = FALSE)
    }
    else{
      Gl1 <- mapply(function(n, y)
                    lapply(y, function(x) x * log(n)),
                    N, Y, SIMPLIFY = FALSE)
      Gl1 <- mapply(function(gl1, pond)
                    mapply("*", gl1, pond, SIMPLIFY = FALSE),
                    Gl1, Pond, SIMPLIFY = FALSE)
      Gl1 <- lapply(Gl1, function(x) Reduce("+", x))
      Gl2 <- mapply(function(n, pl) - pl * log(n),
                    N, Pl, SIMPLIFY = FALSE)
    }      
    Gl <- mapply("+", Gl1, Gl2)
    if (un.nest.el) Gl <- apply(Gl, 1, sum)
    gradi <- opposite * weights * cbind(Gb, Gl)
    attr(lnl, "gradi") <- gradi
    attr(lnl, "gradient") <- apply(gradi, 2, sum)
  }
  if (step < steptol) lnl <- NULL
  else{
    attr(lnl, "probabilities") <- L
    attr(lnl, "step") <- step
  }
  lnl
}

lnl.hlogit <- function(param, X, y, weights = NULL, gradient = FALSE,
                       hessian = FALSE, opposite = TRUE, sumlnl = TRUE,
                       direction = rep(0, length(param)), initial.value = NULL, steptol = 1E-01,
                       rn, choice){
  otime <- proc.time()
  
  opposite <- ifelse(opposite, -1, +1)
  if (is.null(weights)) weights <- 1
  balanced <- TRUE
  u <- rn$nodes
  w <- rn$weights
  K <- ncol(X[[1]])
  n <- nrow(X[[1]])
  J <- length(X)

  step <- 2
  repeat{
    step <- step / 2
    if (step < steptol) break
    beta <- param[1:K] + step * direction[1:K]
    theta <- c(1, param[-c(1:K)]) + step * c(0, direction[-c(1:K)])
    V <- lapply(X, function(x) as.numeric(crossprod(t(x), beta)))
    Vi <- Reduce("+", (mapply("*", V, y, SIMPLIFY = FALSE)))
    DVi <- lapply(V, function(x) Vi - x)
    names(theta) <- levels(choice)
    thetai <- theta[choice]
    alpha <- mapply(function(dvi, th)
                    exp( - (dvi - thetai %o% log(u)) / th),
                    DVi, theta, SIMPLIFY = FALSE)
    A <- t( t(Reduce("+", alpha)))
    G <- exp(- A)
    P <- apply(t(t(G) * w * exp(u)), 1, sum)
    lnl <- sum (opposite * weights * log(P))
    if (is.null(initial.value) || lnl <= initial.value) break
  }
  if (gradient){
    Xi <- Reduce("+", mapply("*", X, y, SIMPLIFY = FALSE))
    DX <- lapply(X, function(x) x - Xi)
    DXt <- mapply("/", DX, theta, SIMPLIFY = FALSE)
    Gb <- lapply(alpha, function(a) apply(t(t(a * G) * w * exp(u)), 1, sum))
    Gb <- mapply(function(a, dxt) a * dxt,
                 Gb, DXt, SIMPLIFY = FALSE)
    Gb <- - Reduce("+", Gb)
    
    Gtj <- mapply(function(a, th) a * log(a) / th,
                  alpha, theta, SIMPLIFY = FALSE)

    Gtl <- mapply(function(a, th) - t( t(a / th) * log(u)),
                  alpha, theta, SIMPLIFY = FALSE)
    Gtl <- Reduce("+", Gtl)
    Gt <- mapply(function(gtj, ay) gtj  + Gtl * ay,
                 Gtj, y, SIMPLIFY = FALSE)
    
    Gt <- sapply(Gt, function(x) apply(t(t(x * G) * exp(u) * w ), 1, sum))

    gradi <- opposite * cbind(Gb, Gt[, -1]) / P
    attr(lnl, "gradi") <- gradi
    attr(lnl, "gradient") <- if (is.matrix(gradi)) apply(gradi, 2, sum) else sum(gradi)
  }
  if (step < steptol) lnl <- NULL
  else{
    attr(lnl, "probabilities") <- P
    attr(lnl, "step") <- step
  }
  lnl
}

