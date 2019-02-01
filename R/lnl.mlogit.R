lnl.slogit <- function(param, X, y, weights = NULL, gradient = FALSE,
                       hessian = FALSE, opposite, direction = rep(0, length(param)),
                       initial.value = NULL,stptol = 1E-01){
    step <- 2
    repeat{
        step <- step / 2
        if (step < stptol) break
        Xb <- lapply(X, function(x) crossprod(t(x), param + step * direction))
        eXb <- lapply(Xb, exp)
        seXb <- suml(eXb)
        P <- lapply(eXb, function(x){v <- x / seXb; v[is.na(v)] <- 0; as.vector(v)})
        Pch <- Reduce("+", mapply("*", P, y, SIMPLIFY = FALSE))
        names(Pch) <- attr(y, "chid")
        lnl <- sum(opposite * weights * log(Pch))
        if (is.null(initial.value) || lnl <= initial.value) break
    }
    if (gradient | hessian) PX <- suml(mapply("*", X, P, SIMPLIFY = FALSE))
    if (gradient){
        Xch <- suml(mapply("*", X, y, SIMPLIFY = FALSE))
        gradi <-  opposite * weights * (Xch - PX)
        attr(lnl, "gradi") <- gradi
        attr(lnl, "gradient") <- if (is.matrix(gradi)) apply(gradi, 2, sum) else sum(gradi)
    }
    if (hessian){
        XmPX <- lapply(X, function(x){g <- x - PX; g[is.na(g)] <- 0; g})
        hessian <-   - suml( mapply(function(x, y) crossprod(x * y, y),
                                    P, XmPX, SIMPLIFY = FALSE))
        attr(lnl, "hessian") <- opposite * hessian
    }
    if (step < stptol) lnl <- NULL
    else{
        Xb <- Reduce("cbind", Xb)
        P <- Reduce("cbind", P)
        dimnames(P) <- dimnames(Xb) <- list(attr(y, "chid"), names(y))
        attr(lnl, "probabilities") <- P
        attr(lnl, "linpred") <- Xb
        attr(lnl, "fitted") <- Pch
        attr(lnl, "step") <- step
    }
    lnl
}


lnl.wlogit <- function(param, X, y, Xs, weights = NULL, gradient = FALSE,
                       hessian = FALSE, opposite, direction = rep(0, length(param)),
                       initial.value = NULL,stptol = 1E-01){
    balanced <- FALSE
    step <- 2
    repeat{
        step <- step / 2
        if (step < stptol) break
        Ks <- ncol(Xs)
        K <- ncol(X[[1]])
        beta <- param[1:K] + step * direction[1:K]
        lambda <- param[(K + 1):(K + Ks)] + direction[(K + 1):(K + Ks)]
        sigma <- 1 + as.numeric(Xs %*% lambda)
        fpsigma <- Xs
        Xb <- lapply(X, function(x) crossprod(t(x / sigma), beta))
        eXb <- lapply(Xb, exp)
        seXb <- suml(eXb)
        P <- lapply(eXb, function(x){v <- x / seXb; v[is.na(v)] <- 0; as.vector(v)})
        Pch <- Reduce("+", mapply("*", P, y, SIMPLIFY = FALSE))
        names(Pch) <- attr(y, "chid")
        lnl <- sum(opposite * weights * log(Pch))
        if (is.null(initial.value) || lnl <= initial.value) break
    }

    if (gradient){
        PX <- suml(mapply("*", X, P, SIMPLIFY = FALSE))
        Xch <- suml(mapply("*", X, y, SIMPLIFY = FALSE))
        gradi <-  opposite * weights * (Xch - PX)
        gradsi <- sapply(as.data.frame(fpsigma),
                         function(x) apply(x * t(t(gradi) * beta), 1, sum)) / (- sigma ^ 2)
        gradi <- gradi / sigma
        gradi <- cbind(gradi, gradsi)
        attr(lnl, "gradi") <- gradi
        attr(lnl, "gradient") <- if (is.matrix(gradi)) apply(gradi, 2, sum) else sum(gradi)
    }
    if (step < stptol) lnl <- NULL
    else{
        P <- Reduce("cbind", P)
        Xb <- log(Reduce("cbind", eXb))
        colnames(P) <- colnames(Xb) <- names(y)
        attr(lnl, "linpred") <- Xb
        attr(lnl, "probabilities") <- P
        attr(lnl, "fitted") <- Pch
        attr(lnl, "step") <- step
    }
    lnl
}

lnl.nlogit <- function(param, X, y, weights = NULL, gradient = FALSE,
                       hessian = FALSE, opposite = TRUE, initial.value = NULL,
                       direction = rep(0, length(param)), stptol = 1E-01,
                       nests, un.nest.el = FALSE, unscaled = FALSE){

    if (un.nest.el){
        lambda <- param[length(param)]
        param <- c(param, rep(lambda, length(nests) - 1))
    }
    thealts <- names(X)
    posalts <- lapply(thealts, function(x) which(unlist(nests) %in% x))
    K <- ncol(X[[1]])
    n <- nrow(X[[1]])
    J <- length(nests)
    Xn <- lapply(nests, function(x) X[x])
    Y <- lapply(nests, function(x) y[x])
    Yn <- lapply(Y, suml)
    step <- 2
    repeat{
        step <- step / 2
        if (step < stptol) break
        beta <- param[1:K] + step * direction[1:K]
        lambda <- param[-c(1:K)] + step * direction[-c(1:K)]
        names(lambda) <- names(nests)
        Xb <- lapply(X, function(x) as.numeric(crossprod(t(x), beta)))
        Xb <- Reduce("cbind", Xb)
        V <- lapply(Xn, function(x)
            lapply(x, function(y)
                as.numeric(crossprod(t(y), beta))
                )
            )
        if (! unscaled) W <- mapply(function(v, l) lapply(v, function(x) x / l), V, lambda, SIMPLIFY = FALSE)
        else W <- V
        A <- lapply(W, function(x) lapply(x, exp))
        N <- lapply(A, suml)
        Pjl <- mapply(function(a, n) lapply(a, function(x) x / n), A, N, SIMPLIFY = FALSE)
        Pl <- mapply(function(n, l) n^l, N, lambda, SIMPLIFY = FALSE)
        D <- suml(Pl)
        Pl <- lapply(Pl, function(x) x / D)
        P <- mapply(function(pjl, pl) lapply(pjl, function(x) x * pl), Pjl, Pl, SIMPLIFY = FALSE)
        Pch <- mapply(function(p, y) mapply("*", p, y, SIMPLIFY = FALSE), P, Y, SIMPLIFY = FALSE)
        Pch <- suml(lapply(Pch, suml))
        lnl <- sum(opposite * weights * log(Pch))
        if (is.null(initial.value) || lnl <= initial.value) break
    }
    if (gradient){
        ### For overlaping nests
        Pj <- unlist(P, recursive = FALSE)
        Pj <- lapply(posalts, function(x) suml(Pj[x]))
        names(Pj) <- thealts
        proba <- Pj
        Pj <- lapply(nests, function(x) Pj[x])
        Pond <- mapply(function(p, pj) mapply("/", p, pj, SIMPLIFY = FALSE), P, Pj, SIMPLIFY = FALSE)
        Xb <- mapply(function(x, pjl)
            suml(mapply("*", x, pjl, SIMPLIFY = FALSE)),
            Xn, Pjl, SIMPLIFY = FALSE)
        Vb <- mapply(function(v, pjl)
            suml(mapply("*", v, pjl, SIMPLIFY = FALSE)),
            V, Pjl, SIMPLIFY = FALSE)
        
        if (!unscaled)
            Xbtot <- suml(mapply("*", Pl, Xb, SIMPLIFY = FALSE))
        else
            Xbtot <- suml(mapply(function(pl, xb, l) pl * xb * l,
                                 Pl, Xb, lambda, SIMPLIFY = FALSE))
        
        if (!unscaled)
            Gb <- mapply(function(x, xb, l)
                lapply(x, function(z) (z+(l-1) * xb)/l), Xn, Xb, lambda, SIMPLIFY = FALSE)
        else
            Gb <- mapply(function(x, xb, l)
                lapply(x, function(z) (z+(l-1) * xb)), Xn, Xb, lambda, SIMPLIFY = FALSE)
        
        Gb <-  mapply(function(gb, y)
            mapply("*", gb, y, SIMPLIFY = FALSE),
            Gb, Y, SIMPLIFY = FALSE)
        Gb <-  mapply(function(gb, pond)
            mapply("*", gb, pond, SIMPLIFY = FALSE),
            Gb, Pond, SIMPLIFY = FALSE)
        Gb <- suml(lapply(Gb, suml)) - Xbtot
        
        if (!unscaled){
            Gl1 <- mapply(function(v, n, vb, l)
                lapply(v, function(z) -(z - l^2 * log(n) + (l-1) * vb)/l^2),
                V, N, Vb, lambda, SIMPLIFY = FALSE)
            Gl1 <- mapply(function(gl1, y) mapply("*", gl1, y, SIMPLIFY = FALSE),
                          Gl1, Y, SIMPLIFY = FALSE)
            Gl1 <- mapply(function(gl1, pond)
                mapply("*", gl1, pond, SIMPLIFY = FALSE),
                Gl1, Pond, SIMPLIFY = FALSE)
            mylog <- function(x){
                nullx <- abs(x) < 1E-20
                x[nullx] <- 0
                x[!nullx] <- log(x[!nullx])
                x
            }
            Gl1 <- lapply(Gl1, suml)
            Gl2 <- mapply(function(vb, n, l, pl) - pl * (l^2 * mylog(n)- l * vb)/ l^2,
                          Vb, N, lambda, Pl, SIMPLIFY = FALSE)
            # log is replaced by mylog which replace log(0) by 0.

        }
        else{
            Gl1 <- mapply(function(n, y)
                lapply(y, function(x) x * log(n)),
                N, Y, SIMPLIFY = FALSE)
            Gl1 <- mapply(function(gl1, pond)
                mapply("*", gl1, pond, SIMPLIFY = FALSE),
                Gl1, Pond, SIMPLIFY = FALSE)
            Gl1 <- lapply(Gl1, suml)
            Gl2 <- mapply(function(n, pl) - pl * log(n),
                          N, Pl, SIMPLIFY = FALSE)
        }
        Gl <- mapply("+", Gl1, Gl2)
        if (un.nest.el) Gl <- apply(Gl, 1, sum)
        gradi <- opposite * weights * cbind(Gb, Gl)
        attr(lnl, "gradi") <- gradi
        attr(lnl, "gradient") <- apply(gradi, 2, sum)
    }
    if (step < stptol) lnl <- NULL
    else{
        P <- unlist(P, recursive = FALSE)
        P <- sapply(posalts, function(x) suml(P[x]))
#        colnames(P) <- thealts
        Xb <- Reduce("cbind",
                     lapply(X, function(x) as.numeric(crossprod(t(x), beta))))
        dimnames(P) <- dimnames(Xb) <- list(attr(y, "chid"), names(y))
        attr(lnl, "probabilities") <- P
        attr(lnl, "linpred") <- Xb
        attr(lnl, "fitted") <- Pch
        attr(lnl, "step") <- step
        
    }
    lnl
} 

lnl.hlogit <- function(param, X, y, weights = NULL,
                       gradient = FALSE, hessian = FALSE, opposite = TRUE,
                       direction = rep(0, length(param)), initial.value = NULL, stptol = 1E-01,
                       rn){
    choice <- Reduce("cbind", y)
    choice <- factor(apply(choice * col(choice), 1, sum), labels = names(y), levels = 1:length(y))
    names(choice) <- NULL
    balanced <- TRUE
    u <- rn$nodes
    w <- rn$weights
    K <- ncol(X[[1]])
    n <- nrow(X[[1]])
    J <- length(X)
    
    step <- 2
    repeat{
        step <- step / 2
        if (step < stptol) break
        beta <- param[1:K] + step * direction[1:K]
        theta <- c(1, param[-c(1:K)]) + step * c(0, direction[-c(1:K)])
        V <- lapply(X, function(x) as.numeric(crossprod(t(x), beta)))
        Vi <- suml(mapply("*", V, y, SIMPLIFY = FALSE))
        DVi <- lapply(V, function(x) Vi - x)
        names(theta) <- levels(choice)
        thetai <- theta[choice]
        alpha <- mapply(function(dvi, th){
            exp( - (dvi - thetai %o% log(u)) / th)},
            DVi, theta, SIMPLIFY = FALSE)
        A <- suml(alpha)
        G <- exp(- A)
        P <- apply(t(t(G) * w * exp(u)), 1, sum)
        lnl <- sum (opposite * weights * log(P))
        if (is.null(initial.value) || lnl <= initial.value) break
    }
    if (gradient){
        Xi <- suml(mapply("*", X, y, SIMPLIFY = FALSE))
        DX <- lapply(X, function(x) x - Xi)
        DXt <- mapply("/", DX, theta, SIMPLIFY = FALSE)
        Gb <- lapply(alpha, function(a) apply(t(t(a * G) * w * exp(u)), 1, sum))
        Gb <- mapply(function(a, dxt) a * dxt,
                     Gb, DXt, SIMPLIFY = FALSE)
        Gb <- - suml(Gb)
        Gtj <- mapply(function(a, th) a * log(a) / th,
                      alpha, theta, SIMPLIFY = FALSE)
        Gtl <- mapply(function(a, th) - t( t(a / th) * log(u)),
                      alpha, theta, SIMPLIFY = FALSE)
        Gtl <- suml(Gtl)
        Gtj <- lapply(Gtj, function(x){x[is.na(x)] <- 0;x})
        Gt <- mapply(function(gtj, ay) gtj  + Gtl * ay,
                     Gtj, y, SIMPLIFY = FALSE)
        Gt <- sapply(Gt, function(x) apply(t(t(x * G) * exp(u) * w ), 1, sum))
        gradi <- opposite * cbind(Gb, Gt[, -1]) / P
        attr(lnl, "gradi") <- gradi
        attr(lnl, "gradient") <- if (is.matrix(gradi)) apply(gradi, 2, sum) else sum(gradi)
    }
    if (step < stptol) lnl <- NULL
    else{
        attr(lnl, "fitted") <- P
        attr(lnl, "step") <- step
        attr(lnl, "linpred") <- Reduce("cbind", V)
    }
    lnl
}

lnl.mprobit <- function(param, y, X, weights = NULL,
                        gradient = FALSE, hessian = FALSE, opposite = TRUE,
                        direction = rep(0, length(param)), initial.value = NULL, stptol = 1E-1,
                        R, seed){

    names(direction) <- names(param)
    mills <- function(x) exp(dnorm(x, log = TRUE) - pnorm(x, log.p = TRUE))
    
    K <- ncol(X[[1]])
    J <- length(X) + 1
    n <- length(y)
    step <- 2
    
    # Mi is a list containing the linear transformation of the utility
    # differences with respect to the first alternative to utility
    # differences with respect to any alternative
    Mi <- list()
    M <- rbind(0, diag(J - 2))
    Mi[[1]] <- diag(J-1)
    for (i in 2:(J-1)){
        Mi[[i]] <- cbind(M[, 0:(i-2), drop = FALSE], -1, M[, ((i-1):(J-2))])
    }
    Mi[[J]] <- cbind(M, -1)
    
    repeat{
        step <- step / 2
        if (step < stptol) break
        beta <- param[1:K] + step * direction[1:K]
        corrCoef <- param[- c(1:K)] + step * direction[- c(1:K)]
        DV <- sapply(X, function(x) crossprod(t(x), beta))
        if (! is.matrix(DV)) DV <- matrix(DV, nrow = 1)
        # Cholesky matrix and covariance matrix of U-U_1
        S <- matrix(0, J - 1, J - 1)
        S[!upper.tri(S)] <- corrCoef
        omega <- S %*% t(S)
        # Si is a list containing the cholesky decomposition of the
        # covariance matrix of utility differences
        Si <- lapply(Mi, function(x) t(chol(x %*% omega %*% t(x))))
        set.seed(seed)
        
        A <- vector("list", length = J - 1)
        ETA <- vector("list", length = J - 2)
        for (id in 1:n){
            ay <- y[id]
            dv <- DV[id, ]
            if (! length(na.omit(dv)) == 0){
            # Jn is the nb of alternative for individual n
                Jn <- length(na.omit(dv)) + 1
                # Compute the random numbers only if the number of alt is at
                # least 3
                if (Jn >= 3) RN <- matrix(runif( (Jn - 2) * R), R, Jn - 2)
                # Compute the E matrix when some alternatives are not available
                if (Jn < J){
                    # insert in the utility diff a 0 for the chosen alternative at
                    # the right place
                    mydv <- DV[id, ]
                    if (ay == J) theTail <- c()
                    else theTail <- dv[ay:(J - 1)]
                    mydv <- c(mydv[0:(ay - 1)], 0, theTail)
                    # check for missing alternatives and 
                    na.alt <- which(is.na(mydv))
                    na.alt[ay < na.alt] <- na.alt[ay < na.alt] - 1
                    E <- (diag(J - 1))[- na.alt,]
                    Min <- E %*% Mi[[ay]]
                    si <- t(chol(Min %*% omega %*% t(Min)))
                    dv <- na.omit(dv)
                }
                else si <- Si[[ay]]
                A[[1]] <- rbind(A[[1]], rep(- dv[1] / si[1, 1], R))
                if (Jn >= 3){
                    eta <- qnorm(RN[, 1, drop = FALSE] * pnorm(A[[1]][id,]))
                    ETA[[1]] <- rbind(ETA[[1]], t(eta))
                    for (l in 2:(Jn - 1)){
                        A[[l]] <- rbind(A[[l]], t(- (dv[l] + eta[, 1:(l-1), drop = FALSE] %*% si[l, 1:(l-1)]) / si[l, l]))
                        if (l != (Jn - 1)){
                            etai <- qnorm(RN[, l] * pnorm(A[[l]][id, ]))
                            ETA[[l]] <- rbind(ETA[[l]], etai)
                            eta <- cbind(eta, etai)
                        }
                    }
                }
                if (Jn < J){
                    for (l in Jn:(J-1)) A[[l]] <- rbind(A[[l]], rep(1, R))
                    for (l in (Jn - 1):(J - 2)) ETA[[l]] <- rbind(ETA[[l]], rep(NA, R))
                }
            }
            else{
                for (l in 1:(J - 1)) A[[l]] <- rbind(A[[l]], rep(NA, R))
                for (l in 1:(J - 2)) ETA[[l]] <- rbind(ETA[[l]], rep(NA, R))
            }
        }
        PR <- lapply(A, pnorm)
        probai <- Reduce("*", PR)
        P <- apply(probai, 1, mean)
        lnl <- opposite * sum(log(P))
        if (is.null(initial.value) || lnl <= initial.value) break
    }
    if (gradient){
        set.seed(seed)
        # Two matrices VK and VL maps s to vec(S) and vec(S')
        M <- matrix(NA, J - 1, J - 1)
        M[! upper.tri(M)] <- 1:(J * (J - 1) / 2)
        m <- c(M)
        mp <- c(t(M))
        Id <- diag(J * (J - 1) / 2)
        VK <- VL <- matrix(0, (J - 1) ^ 2, J * (J - 1) / 2)
        for (i in 1:length(m)){
            if (!is.na(m[i]))  VL[i, ] <- Id[m[i],]
            if (!is.na(mp[i])) VK[i, ] <- Id[mp[i],]
        }
        VK <- t(VK)
        VL <- t(VL)
        DB <- c()
        DS <- c()
        
        pos <- matrix(0, J - 1, J- 1)
        pos[!upper.tri(pos)] <- 1:(J * (J - 1) / 2)
        for (id in 1:n){
            ay <- y[id]
            dv <- DV[id, ]
            # Jn is the nb of alternative for individual n
            Jn <- length(na.omit(dv)) + 1
            # Compute the random numbers only if the number of alt is at
            # least 3
            if (Jn >= 3) RN <- matrix(runif( (Jn - 2) * R), R, Jn - 2)
            # Compute the E matrix when some alternatives are not available
            if (Jn < J){
                # insert in the utility diff a 0 for the chosen alternative at
                # the right place
                mydv <- DV[id, ]
                if (ay == J) theTail <- c()
                else theTail <- dv[ay:(J - 1)]
                mydv <- c(mydv[0:(ay - 1)], 0, theTail)
                # check for missing alternatives and 
                na.alt <- which(is.na(mydv))
                na.alt.rm.rows <- na.alt
                na.alt.rm.rows[ay < na.alt] <- na.alt[ay < na.alt] - 1
                E <- (diag(J - 1))[- na.alt.rm.rows,]
                Min <- E %*% Mi[[ay]]
                si <- t(chol(Min %*% omega %*% t(Min)));
                odv <- dv
                dv <- na.omit(dv)
                # Two matrices VKn and VLn maps sn to vec(Sn) and vec(Sn')
                M <- matrix(NA, Jn - 1, Jn - 1)
                M[! upper.tri(M)] <- 1:(Jn * (Jn - 1) / 2)
                m <- c(M)
                mp <- c(t(M))
                Id <- diag(Jn * (Jn - 1) / 2)
                VKn <- VLn <- matrix(0, (Jn - 1) ^ 2, Jn * (Jn - 1) / 2)
                for (i in 1:length(m)){
                    if (!is.na(m[i]))  VLn[i, ] <- Id[m[i],]
                    if (!is.na(mp[i])) VKn[i, ] <- Id[mp[i],]
                }
                VKn <- t(VKn)
                VLn <- t(VLn)
                posn <- matrix(0, Jn - 1, Jn- 1)
                posn[!upper.tri(posn)] <- 1:(Jn * (Jn - 1) / 2)
            }
            else{
                si <- Si[[ay]]
                Min <- Mi[[ay]]
                posn <- pos
                VKn <- VK
                VLn <- VL
                na.alt <- an.alt.rm.rows <- c()
            }
            JacB <- ((Min %*% S) %x% Min) %*% t(VL) + (Min %x% (Min %*% S)) %*% t(VK)
            JacA <- (si %x% diag(Jn - 1)) %*% t(VLn) + (diag(Jn - 1) %x% si) %*% t(VKn)
            Jac <- t(JacB) %*% t(ginv(JacA))
            Gr <- c()
            theAs <- vector("list", length = Jn * (Jn - 1) / 2)
            for (u in 1:(Jn * (Jn - 1) / 2)) theAs[[u]] <- matrix(0, R, Jn - 1)
            Xi <- lapply(X, function(x) matrix(x[id, ], R, K, byrow = TRUE))
            if (Jn < J) Xi <- Xi[- na.alt.rm.rows]
            Abeta <- vector("list", length = Jn - 1)
            Abeta[[1]] <- - Xi[[1]] / si[1, 1]
            Dbeta <- Abeta[[1]] * mills(A[[1]][id, ])
            if (Jn >= 3){
                for (j in 2:(Jn - 1)){
                    Abeta[[j]] <- - Xi[[j]] / si[j, j]
                    for (k in 1:(j-1))
                        Abeta[[j]] <- Abeta[[j]] - si[j, k] / si[j, j] * RN[, k] * dnorm(A[[k]][id, ]) /
                            dnorm(ETA[[k]][id, ]) * Abeta[[k]]
                    Dbeta <- Dbeta + mills(A[[j]][id, ]) * Abeta[[j]]
                }
            }
            Dbeta <- apply(Dbeta * probai[id,], 2, mean) / P[id]
            DB <- rbind(DB, Dbeta)
            s <- c(t(si)[!lower.tri(si)])
            As <- theAs
            # si i = j = l (-> k)
            for (k in 1:(Jn - 1)) As[[posn[k, k]]][, k] <- - A[[k]][id, ] / si[k, k]
           # si j < (i = l)
            if (Jn >2){
                for (i in 2:(Jn - 1))
                    for (j in 1:(i - 1)) As[[posn[i, j]]][, i] <- - ETA[[j]][id, ] / si[i, i]
            }
            # i = 1 => j = 1     => 1 < l => l = 2,3   => l = (i+1):(J-1)
            # i = 2 => j = 1,2   => 2 < l => l = 3     => l = (i+1):(J-1)
            # i = 3 => j = 1,2,3 => 3 < l => l = rien
            # pour Jn = 3; i = 1:1; j = 1:1, l= 2:2
            if (Jn >= 3){
                for (i in 1:(Jn - 2)){
                    for (j in 1:i){
                        for (l in (i+1):(Jn-1)){
                            for (h in 1:(l-1)){
                                As[[posn[i, j]]][, l] <- As[[posn[i, j]]][, l] -
                                    RN[, h] * si[l, h] / si[l, l] * dnorm(A[[h]][id, ]) / dnorm(ETA[[h]][id,]) * As[[pos[i,j]]][, h]
                            }
                        }
                    }
                }
            }
            lambda <- sapply(A[1:(Jn - 1)], function(x) mills(x[id, ]))
            Ds <- sapply(As, function(x) x * lambda)
            Dss <- Ds[1:R, , drop = FALSE]
            if (Jn > 2)
                for (l in 2:(Jn - 1)) Dss <- Dss + Ds[(R * (l - 1) + 1:R), ]
            Ds <- Dss
            Ds <- apply(Ds * probai[id,], 2, mean) / P[id]
            Ds <- as.numeric(Jac %*% Ds)
            DS <- rbind(DS, Ds)
        }
        Gr <- cbind(DB, DS)
        colnames(Gr) <- c(names(beta), names(corrCoef))
        attr(lnl, "gradi") <- opposite * Gr
        attr(lnl, "gradient") <- opposite * apply(Gr, 2, sum)
    }
    if (step < stptol) lnl <- NULL
    else{
        attr(lnl, "fitted") <- P
        attr(lnl, "step") <- step
    }
    lnl
}

