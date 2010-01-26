lnl.mlogits <- function(param, X, y, weights = NULL, gradient = FALSE,
                        hessian = FALSE, opposite = TRUE, sumlnl = TRUE,
                        direction = 0, initial.value = NULL){

  opposite <- ifelse(opposite, -1, +1)
  if (is.null(weights)) weights <- 1
  balanced <- TRUE
  
  # compute the probabilities and the likelihood for step = 1
  step <- 1
  eXb <- lapply(X, function(x) exp(crossprod(t(x), param + direction)))
  seXb <- suml(eXb)
  P <- lapply(eXb, function(x){v <- x/seXb; v[is.na(v)] <- 0; as.vector(v)})
  Pch <- suml(mapply("*", P, y, SIMPLIFY = FALSE))
  lnl <- sum(opposite * weights * log(Pch))
  
  if (!is.null(initial.value)){
    while(lnl > initial.value){
      step <- step / 2
      eXb <- lapply(X, function(x) exp(crossprod(t(x), param + step * direction)))
      seXb <- suml(eXb)
      P <- lapply(eXb, function(x){v <- x/seXb; v[is.na(v)] <- 0; as.vector(v)})
      Pch <- suml(mapply("*", P, y, SIMPLIFY = FALSE))
      lnl <- sum(opposite * weights * log(Pch))
    }
  }

  if (gradient | hessian) PX <- suml(mapply("*", X, P, SIMPLIFY = FALSE))
  if (gradient){
    Xch <- suml(mapply("*", X, y, SIMPLIFY = FALSE))
    gradi <-  opposite * weights * (Xch - PX)
    attr(lnl, "gradi") <- gradi
    attr(lnl, "gradient") <- if (is.matrix(gradi)) apply(gradi,2,sum) else sum(gradi)
  }
  if (hessian){
    if (balanced)
      XmPX <- lapply(X, function(x) x - PX)
    else
      XmPX <- lapply(X, function(x){g <- x - PX; g[is.na(g)] <- 0; g})
    hessian <-   - suml( mapply(function(x, y) crossprod(x*y, y),
                               P, XmPX, SIMPLIFY = FALSE))
    attr(lnl, "hessian") <- opposite * hessian
  }
  attr(lnl, "probabilities") <- P
  attr(lnl, "step") <- step
  lnl
}

lnl.nlogit <- function(param, X, y, weights = NULL, gradient = FALSE, hessian = FALSE,
                       opposite = TRUE, sumlnl = TRUE, nests,
                       direction = rep(0, length(param)), initial.value = NULL){
  opposite <- ifelse(opposite, -1, +1)
  if (is.null(weights)) weights <- 1
  K <- ncol(X[[1]])
  n <- nrow(X[[1]])
  J <- length(nests)
  nestv <- rep(names(nests), sapply(nests, length))
  names(nestv) <- unlist(nests)
  nestv <- nestv[names(X)]
  Y <- as.matrix(as.data.frame(y))

  step <- 1
  beta <- param[1:K] + direction[1:K]
  lambda <- param[-c(1:K)] + direction[-c(1:K)]
  names(lambda) <- names(nests)
  V <- sapply(X, function(x) as.numeric(crossprod(t(x),beta)))
  A <- exp(t(t(V)/lambda[nestv]))
  N <- c()
  ynest <- c()
  for (i in 1:J){
    anestname <- names(nests)[i]
    anest <- nests[[i]]
    N <- cbind(N, apply(A[, anest, drop = FALSE], 1, sum))
    ynest <- cbind(ynest, apply(Y[, anest, drop = FALSE], 1, sum))
    colnames(N)[ncol(N)] <- colnames(ynest)[ncol(ynest)] <- anestname
  }
  Denom <- apply(t( t(N)^lambda ), 1, sum)
  P <- A * t(t(N[, nestv, drop = FALSE])^(lambda[nestv]-1)) / Denom
  lnl <- sum(opposite * weights * log(apply(P * Y, 1, sum)))
  if (!is.null(initial.value)){
    while(lnl > initial.value){
      step <- step / 2
      beta <- param[1:K] + step * direction[1:K]
      lambda <- param[-c(1:K)] + step * direction[-c(1:K)]
      names(lambda) <- names(nests)
      V <- sapply(X, function(x) as.numeric(crossprod(t(x),beta)))
      A <- exp(t(t(V)/lambda[nestv]))
      N <- c()
      ynest <- c()
      for (i in 1:J){
        anestname <- names(nests)[i]
        anest <- nests[[i]]
        N <- cbind(N, apply(A[, anest, drop = FALSE], 1, sum))
        ynest <- cbind(ynest, apply(Y[, anest, drop = FALSE], 1, sum))
        colnames(N)[ncol(N)] <- colnames(ynest)[ncol(ynest)] <- anestname
      }
      Denom <- apply(t( t(N)^lambda ), 1, sum)
      P <- A * t(t(N[, nestv, drop = FALSE])^(lambda[nestv]-1)) / Denom
      lnl <- sum(opposite * weights * log(apply(P * Y, 1, sum)))
    }
  }

  if (gradient){
    lambdanestl <- as.list(lambda)
    lambdal <- as.list(lambda[nestv])
    Al <- as.list(data.frame(A))
    Vl <- as.list(data.frame(V))
    Nl <- as.list(data.frame(N))
    Gbeta <- mapply(function(x, y)  x/y  , X, lambdal, SIMPLIFY = FALSE)
    Glamb <- mapply(function(x, y) -x/y^2, Vl, lambdal, SIMPLIFY = FALSE)
    Abeta <- mapply("*", Gbeta, Al, SIMPLIFY = FALSE)
    Alamb <- mapply("*", Glamb, Al, SIMPLIFY = FALSE)
    
    Nbeta <- list()
    Nlamb <- list()
    ynest <- list()
    for (i in 1:J){
      anestname <- names(nests)[i]
      anest <- nests[[i]]
      Nbeta[[i]] <- suml(Abeta[anest])
      Nlamb[[i]] <- suml(Alamb[anest])
      ynest[[i]] <- suml(y[anest])
    }
    names(Nbeta) <- names(Nlamb) <- names(ynest) <- names(nests)
    Dbeta <- suml(mapply(function(x, y, z) x*y^(x-1)*z,
                         lambdanestl, Nl, Nbeta,
                         SIMPLIFY = FALSE))
    Dlamb <- mapply(function(x, y, z) y^x*(log(y)+x/y*z),
                    lambdanestl, Nl, Nlamb,
                    SIMPLIFY = TRUE)
    Nbeta <- mapply(function(x, y, z) (x-1)/y*z,
                 lambdanestl, Nl, Nbeta,
                 SIMPLIFY = FALSE)
    Nlamb <- mapply(function(x, y, z) log(y)+(x-1)/y*z,
                 lambdanestl, Nl, Nlamb,
                 SIMPLIFY = FALSE)
    Gbetai <- suml(mapply("*", Gbeta, y, SIMPLIFY = FALSE))
    Nbetai <- suml(mapply("*", Nbeta, ynest, SIMPLIFY = FALSE))
    Nlambi <- mapply("*", Nlamb, ynest, SIMPLIFY = TRUE) 
    Glambi <- mapply("*", Glamb, y, SIMPLIFY = TRUE)   

    Glambi2 <- c()
    for (i in 1:J){
      anestname <- names(nests)[i]
      anest <- nests[[i]]
      Glambi2 <- cbind(Glambi2, apply(Glambi[, anest, drop = FALSE], 1, sum))
      colnames(Glambi2)[ncol(Glambi2)] <- anestname
    }
    gradlambi <- Glambi2+Nlambi-Dlamb/Denom
    gradbetai <- Gbetai + Nbetai - Dbeta/Denom
    gradi <- opposite * weights * cbind(gradbetai, gradlambi)
    attr(lnl, "gradi") <- gradi
    attr(lnl, "gradient") <- apply(gradi, 2, sum)
    attr(lnl, "step") <- step
  }
  attr(lnl, "probabilities") <- P
  lnl
}

lnl.hlogit <- function(param, X, y, weights = NULL, gradient = FALSE,
                       hessian = FALSE, opposite = TRUE, sumlnl = TRUE,
                       direction = rep(0, length(param)), initial.value = NULL, rn, choice){
  opposite <- ifelse(opposite, -1, +1)
  if (is.null(weights)) weights <- 1
  balanced <- TRUE
  u <- rn$nodes
  w <- rn$weights
  K <- ncol(X[[1]])
  n <- nrow(X[[1]])
  J <- length(X)

  step <- 1
  beta <- param[1:K] + direction[1:K]
  theta <- c(1,param[-c(1:K)]) + c(0, direction[-c(1:K)])
  V <- lapply(X,function(x) as.numeric(crossprod(t(x),beta)))
  Vi <- suml(mapply("*",V,y,SIMPLIFY=FALSE))
  DVi <- sapply(V,function(x) Vi-x)
  names(theta) <- levels(choice)
  thetai <- theta[choice]
  DVi[DVi==0] <- NA
  DVi <- lapply(u,function(x) t(-t(DVi-thetai*log(x))/theta) )
  alpha <- lapply(DVi,exp)
  A <- lapply(alpha,function(x) apply(x,1,function(x){sum(x,na.rm=TRUE)}))
  G <- lapply(A,function(x) exp(-x))
  P <- suml(mapply("*",w,G,SIMPLIFY=FALSE))
  lnl <- sum (opposite * weights * log(P))
  if (!is.null(initial.value)){
    while(lnl > initial.value){
      step <- step / 2
      beta <- param[1:K] + step * direction[1:K]
      theta <- c(1,param[-c(1:K)]) + step * c(0, direction[-c(1:K)])
      V <- lapply(X,function(x) as.numeric(crossprod(t(x),beta)))
      Vi <- suml(mapply("*",V,y,SIMPLIFY=FALSE))
      DVi <- sapply(V,function(x) Vi-x)
      names(theta) <- levels(choice)
      thetai <- theta[choice]
      DVi[DVi==0] <- NA
      DVi <- lapply(u,function(x) t(-t(DVi-thetai*log(x))/theta) )
      alpha <- lapply(DVi,exp)
      A <- lapply(alpha,function(x) apply(x,1,function(x){sum(x,na.rm=TRUE)}))
      G <- lapply(A,function(x) exp(-x))
      P <- suml(mapply("*",w,G,SIMPLIFY=FALSE))
      lnl <- sum (opposite * weights * log(P))
    }
  }
  if (gradient){
    Xi <- suml(mapply("*",X,y,SIMPLIFY=FALSE))
    DX <- lapply(X,function(x) Xi-x)
    DX <- array(unlist(DX),dim=c(n,K,J))
    DX <- aperm(DX,c(1,3,2))
    alphaDtheta <- lapply(alpha, function(x) t( t(x)/theta))
    alphaDthetaR <- lapply(alphaDtheta,function(x) array(rep(x,K),dim=c(n,J,K)))
    Dbeta <- lapply(alphaDthetaR,function(x) -x*DX)
    Dbeta <- lapply(Dbeta,function(x) apply(x,c(1,3),function(x) sum(x,na.rm=TRUE)))
    ym <- matrix(unlist(y),ncol=J)
    alphaDthetaTDVi <- mapply(function(x,y) -log(x)*y,alpha,alphaDtheta,SIMPLIFY=FALSE)
    alphaDthetaTDVi <- lapply(alphaDthetaTDVi,function(x){x[is.na(x)]=0;x})
    SalphaDtheta <- lapply(alphaDtheta,function(x)
                           matrix(rep(apply(x,1,function(x) sum(x,na.rm=TRUE)),J),ncol=J))
    SalphaDthetalnu <- mapply(function(x,y) x*log(y),SalphaDtheta,u,SIMPLIFY=FALSE)
    Dtheta <- mapply("+",
                     lapply(alphaDthetaTDVi,function(x) x*(!ym)),
                     lapply(SalphaDthetalnu,function(x) x*ym),
                     SIMPLIFY=FALSE)
    Dtheta <- lapply(Dtheta,function(x) x[,-1])
    DD <- mapply(cbind,Dbeta,Dtheta,SIMPLIFY=FALSE)
#    DD <- lapply(DD,function(x){colnames(x) <- names(param);x})
    DD <- suml(mapply(function(x,y,z) x*y*z,G,DD,w,SIMPLIFY=FALSE))
    colnames(DD) <- names(param)
    gradi <- - opposite * (DD/P)
    attr(lnl, "gradi") <- gradi
    attr(lnl, "gradient") <- if (is.matrix(gradi)) apply(gradi,2,sum) else sum(gradi)
  }
  attr(lnl, "probabilities") <- P
  attr(lnl, "step") <- step
  lnl
}

lnl.rlogit <- function(param, y, Xa, Xc,
                       weights = NULL, gradient = FALSE, hessian = FALSE,
                       opposite = TRUE, sumlnl = TRUE,
                       direction = rep(0, length(param)), initial.value = NULL,
                       Varc, Vara, random.nb,
                       id, rpar, correlation){

  if (is.null(weights)) weights <- rep(1, length(y[[1]]))
  opposite <- ifelse(opposite, -1, +1)
  step <- 1
  betac <- param[Varc] + direction[Varc]
  mua <- param[Vara] + direction[Vara]
  Kc <- length(Varc)
  Ka <- length(Vara)
  n <- length(y[[1]])
  K <- Kc + Ka
  R <- nrow(random.nb)
  siga <- param[-c(1:K)] + direction[-c(1:K)]
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

  if (!is.null(initial.value)){
    while(lnl > initial.value){
      step <- step / 2
      betac <- param[Varc] + step * direction[Varc]
      mua <- param[Vara] + step * direction[Vara]
      siga <- param[-c(1:K)] + step * direction[-c(1:K)]
      names(mua) <- names(siga) <- colnames(Xa[[1]])
      b <- make.beta(mua, siga, rpar, random.nb, correlation)
      betaa <- b$betaa
      A <- lapply(Xc, function(x) as.vector(crossprod(t(as.matrix(x)), betac)))
      B <- lapply(Xa, function(x) tcrossprod(as.matrix(x), betaa))
      AB <- mapply(function(x, y) exp(x + y), A, B, SIMPLIFY = FALSE)
      S <- suml(AB)
      P <- lapply(AB, function(x) x/S)
      Pch <- suml(mapply("*", P, y, SIMPLIFY=FALSE))
      if (!is.null(id)){
        Pch <- apply(Pch, 2, tapply, id, prod)
      }
      pm <- apply(Pch, 1, mean)
      if (!is.null(id)) lnl <- opposite * sum(weights[!duplicated(id)] * log(pm))
      else lnl <- opposite * sum(weights * log(pm))
    }
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
  attr(lnl, "probabilities") <- probabilities
  lnl
}
    
  
