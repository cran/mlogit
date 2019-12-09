names.rpar <- function(rpar, prefix = NULL, diag = NULL, unique = FALSE){
    K <- length(rpar)
    nms <- vector(mode = "character", length = K * (K + 1) / 2)
    pos <- 0
    for (i in 1:K){
        for (j in 1:i){
            pos <- pos + 1
            if (is.null(prefix)) nms[pos] <- paste(rpar[j], ":", rpar[i], sep = "")
            else{
                if (is.null(diag)) nms[pos] <- paste(prefix, ".", rpar[j], ":", rpar[i], sep = "")
                else{
                    ifelse(i == j,
                           nms[pos] <- paste(diag, ".",
                                             ifelse(unique,
                                                    rpar[i],
                                                    paste(rpar[j], ":", rpar[i], sep = "")), sep = ""),
                           nms[pos] <- paste(prefix, ".", rpar[j], ":", rpar[i], sep = "")
                           )
                }
            }
        }
    }
    nms
}


make.beta <- function(mua, siga, rpar, random.nb, correlation){
    uncorrelated <- setdiff(names(rpar), correlation)
    correlated <-   names(mua)[sort(match(correlation,  names(mua)))]
    uncorrelated <- names(mua)[sort(match(uncorrelated, names(mua)))]
    singlepar <- names(rpar)[rpar %in% c("zbu", "zbt")]
    singlepar <- names(mua)[sort(match(singlepar, names(mua)))]
    utwopars <- setdiff(uncorrelated, singlepar)

    if (is.logical(correlation) && ! correlation) correlation <- character(0)
    nr <- names(mua)
    rpar <- rpar[nr]

    censored <-    nr[rpar == "cn"]
    lognormal <-   nr[rpar == "ln"]
    truncated <-   nr[rpar == "tn"]
    normal  <-     nr[rpar ==  "n"]
    uniform  <-    nr[rpar ==  "u"]
    triangular  <- nr[rpar ==  "t"]

    zbuniform <- nr[rpar == "zbu"]
    zbtriangular <- nr[rpar == "zbt"]

    Ko <- sum(rpar %in% c("zbu", "zbt"))
        
    R <- nrow(random.nb)
    Ku <- length(uncorrelated)
    Kc <- length(correlated)
    Ka <- Kc + Ku

    betaa <- matrix(NA, R, Ka)
    colnames(betaa) <- nr
    betaa.mu <- betaa

    if (Ku){
        random.nbu <- random.nb[, 1:Ku, drop = FALSE]
        if (Ku - Ko){
            betaa.sigmau <- matrix(NA, R, Ku - Ko)       
            colnames(betaa.sigmau) <- paste("sd.", utwopars, sep = "")
        }
        colnames(random.nbu) <- paste("sd.", uncorrelated, sep = "")


        sel <- intersect(uniform, uncorrelated)
        sd.sel <- paste("sd.", sel, sep = "")
        if (length(sel)){
            etauni <- pnorm(random.nbu[, sd.sel, drop = FALSE])
            betaa[, sel] <- t(mua[sel] - siga[sd.sel] + 2 * t(etauni) * siga[sd.sel])
            betaa.mu[, sel] <- 1
            betaa.sigmau[, sd.sel] <- 2 * etauni - 1
        }

        sel <- intersect(zbuniform, uncorrelated)
        if (length(sel)){
            etauni <- pnorm(random.nbu[, paste("sd.", sel, sep = ""), drop = FALSE])
            betaa[, sel] <- 2 * etauni * mua[sel]
            betaa.mu[, sel] <- 2 * etauni
        }

        sel <- intersect(triangular, uncorrelated)
        sd.sel <- paste("sd.", sel, sep = "")
        if (length(sel)){
            eta <- pnorm(random.nbu[, sd.sel, drop = FALSE])
            betaa.mu[, sel] <- 1
            betaa.sigmau[, sd.sel] <- (eta < 0.5) * (sqrt(2 * eta) - 1) +
                (eta > 0.5) * (1 - sqrt(2 * (1 - eta)))
            betaa[, sel] <- t(mua[sel] + siga[sd.sel] * t(betaa.sigmau[, sd.sel]))
        }

        sel <- intersect(zbtriangular, uncorrelated)
        sd.sel <- paste("sd.", sel, sep = "")
        if (length(sel)){
            eta <- pnorm(random.nbu[, sd.sel, drop = FALSE])
            betaa.mu[, sel] <- (eta < 0.5) * sqrt(2 * eta) +
                (eta > 0.5) * (2 - sqrt(2 * (1 - eta)))
            betaa[, sel] <- t(t(betaa.mu[, sel]) * mua[sel])
        }
        
        sel <- intersect(censored, uncorrelated)
        sd.sel <- paste("sd.", sel, sep = "")
        if (length(sel)){
            betaa[, sel] <- pmax(t(mua[sel] + siga[sd.sel] *
                                   t(random.nbu[, sd.sel, drop = FALSE])), 0)
            betaa.mu[, sel] <- as.numeric(betaa[, sel] > 0)
            betaa.sigmau[, sd.sel] <- betaa.mu[, sel] * random.nbu[, sd.sel]
        }

        sel <- intersect(lognormal, uncorrelated)
        sd.sel <- paste("sd.", sel, sep = "")
        if (length(sel)){
            betaa[, sel] <- exp(t(mua[sel] + siga[sd.sel] * t(random.nbu[, sd.sel, drop = FALSE])))
            betaa.mu[, sel] <- betaa[, sel, drop = FALSE]
            betaa.sigmau[, sd.sel] <- betaa.mu[, sel] * random.nbu[, sd.sel]
        }
        sel <- intersect(normal, uncorrelated)
        sd.sel <- paste("sd.", sel, sep = "")
        if (length(sel)){
            betaa[, sel] <- t(mua[sel] + siga[sd.sel] * t(random.nbu[, sd.sel, drop = FALSE]))
            betaa.mu[, sel] <- 1
            betaa.sigmau[, sd.sel] <- random.nbu[, sd.sel, drop = FALSE]
        }        
    }
    if (Kc){
        random.nbc <- random.nb[, (Ku + 1):(Ku + Kc), drop = FALSE]
        names.corr.coef <- names.rpar(correlated, prefix = "chol")
        CC <- ltm(siga[(Ku - Ko + 1):length(siga)], to = "ltm")
        sigeta <- random.nbc %*% t(CC)
        colnames(sigeta) <- correlated
        betaa[, correlated] <- t(mua[correlated] + t(sigeta))
        betaa.mu[, correlated] <- matrix(1, R, length(correlated))
        betaa.sigmac <- random.nbc[, Reduce("c", lapply(1:Kc, function(i) 1:i))]
        colnames(betaa.sigmac) <- names.corr.coef
        for (i in 1:Kc){
#            sigi <- i + cumsum(c(0, (Kc - 1):1))[1:i]
            sigi <- (i * (i - 1) / 2 + 1):(i * (i + 1) / 2)
#            print(sigi)
            if (rpar[i] == "cn"){
                betaa[, i] <- pmax(betaa[, i], 0)
                betaa.mu[, i] <- as.numeric(betaa[, i] > 0)
                betaa.sigmac[, sigi] <- as.numeric(betaa[, i] > 0) * betaa.sigmac[, sigi]
            }
            if (rpar[i] == "ln"){
                betaa[, i] <- exp(betaa[, i])
                betaa.mu[, i] <- betaa[, i]
                betaa.sigmac[, sigi] <- betaa[, i] * betaa.sigmac[, sigi]
            }
        }
    }
    if (Kc &   (Ku - Ko)) betaa.sigma <- cbind(betaa.sigmau, betaa.sigmac)
    if (Kc & ! (Ku - Ko)) betaa.sigma <- betaa.sigmac
    if (! Kc & (Ku - Ko)) betaa.sigma <- betaa.sigmau
    if (! Kc & ! (Ku - Ko)) betaa.sigma <- NULL

    list(betaa = betaa, betaa.mu = betaa.mu, betaa.sigma = betaa.sigma)
}

gnrpoints <- function(low, up, n = 100){
    low + (up - low) * (0:n) / n
}

halton <- function(prime = 3, length = 100, drop = 10){
    halt <- 0
    t <- 0
    while(length(halt) < length + drop){
        t <- t + 1
        halt <- c(halt, rep(halt, prime - 1) + rep(seq(1, prime - 1, 1) / prime ^ t, each = length(halt)))
    }
    halt[(drop + 1):(length + drop)]
}

make.random.nb <- function(R, Ka, halton){
    # Create the matrix of random numbers
    if (! is.null(halton)){
        length.halton <- rep(R,Ka)
        prime <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
                   47, 53, 59, 61, 71, 73, 79, 83, 89, 97, 101, 103,
                   107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167,
                   173, 179, 181, 191, 193, 197, 199)
        drop.halton <- rep(100,Ka)
        if (! is.na(halton) && ! is.null(halton$prime)){
            if (length(halton$prime) != Ka){
                stop("wrong number of prime numbers indicated")
            }
            else{
                prime <- halton$prime
            }
            if (! is.na(halton) && ! is.null(halton$drop)){
                if (! length(halton$drop) %in% c(1, Ka)) stop("wrong number of drop indicated")
                if (length(halton$drop) == 1){
                    drop.halton <- rep(halton$drop, Ka)
                }
                else{
                    drop.halton <- halton$drop
                }
            }
        }
        random.nb <- numeric(0)
        i <- 0
        for (i in 1:Ka){
            random.nb <- cbind(random.nb, qnorm(halton(prime[i], R, drop.halton[i])))
        }
    }
    else{
        random.nb <- matrix(rnorm(R * Ka), ncol = Ka, nrow = R)
    }
    random.nb
}

makeC <- function(x){
    # create the lower triangular C matrix
    K <- (- 1 + sqrt(1 + 8 * length(x))) / 2
    mat <- matrix(0, K, K)
    mat[lower.tri(mat, diag = TRUE)] <- x
    mat
}

make.rpar <- function(rpar, correlation, estimate, norm){
    K <- length(rpar)
    Kc <- length(correlation)
    Ku <- K - Kc
    nr <- names(rpar)
    rpar <- lapply(rpar, function(x) list(dist = x))
    if (Kc){
        Ktot <- length(estimate)
        index.corr <- (Ktot - 0.5 * Kc * (Kc + 1) + 1):Ktot
        v <- estimate[index.corr]
        v <- tcrossprod(ltm(v, to = "ltm"))
        colnames(v) <- rownames(v) <- correlation
        sc <- sqrt(diag(v))
        names(sc) <- correlation
    }
    else sc <- c()
    for (i in (1:K)){
        m <- estimate[nr[i]]
        if (! nr[i] %in% correlation){
            s <- estimate[paste("sd.", nr[i], sep = "")]
        }
        else{
            s <- sc[nr[i]]
        }
        names(m) <- names(s) <- NULL
        rpar[[i]]$mean <- m
        rpar[[i]]$sigma <- s
        rpar[[i]]$name <- nr[[i]]
        if (! is.null(norm)){
            vn <- estimate[norm]
            names(vn) <- NULL
            rpar[[i]]$norm <- vn
        }
    }
    z <- lapply(rpar, function(x){attr(x, "class") = "rpar" ; x})
    if (length(correlation)) attr(z, 'covariance') <- v
    z
}

lnl.rlogit <- function(param, X, y, Xs, weights = NULL,
                       gradient = TRUE, hessian = FALSE, opposite = TRUE,
                       direction = rep(0, length(param)), initial.value = NULL, stptol = 1e-10,
                       R, seed, id, rpar, correlation, halton){
    panel <- ! is.null(id)
    otime <- proc.time()
    K <- ncol(X[[1]])
    Ktot <- length(param)
    N <- nrow(X[[1]])
    J <- length(X)
    
    if (! is.null(Xs)) Ks <- ncol(Xs) else Ks <- 0
    fpsigma <- Xs
    
    if (panel){
        n <- length(unique(id))
        if (length(weights) == 1) weights <- rep(weights, N)
    }
    colnamesX <- colnames(X[[1]])
    uncorrelated <- setdiff(names(rpar), correlation)
    fixedpar <- setdiff(colnamesX, names(rpar))
    singlepar <- names(rpar)[rpar %in% c("zbu", "zbt")]
    utwopars <- intersect(names(rpar)[! rpar %in% c("zbu", "zbt")], uncorrelated)       
    correlated <-   colnamesX[sort(match(correlation,  colnamesX))]
    uncorrelated <- colnamesX[sort(match(uncorrelated, colnamesX))]
    fixedpar <- colnamesX[sort(match(fixedpar, colnamesX))]
    randompar <- colnamesX[sort(match(names(rpar), colnamesX))]
    singlepar <- colnamesX[sort(match(singlepar, colnamesX))]
    utwopars <- colnamesX[sort(match(utwopars, colnamesX))]

    Kc <- length(correlated)
    Ku <- length(uncorrelated)
    Ktot <- length(param)
    Ko <- Ku - length(utwopars)
    Vf <- match(fixedpar, names(param))
    Va <- match(randompar, names(param))
    Vc <- match(correlated, names(param))
    Vu <- match(utwopars, names(param))
    if (! is.null(Xs)) Vs <- (Ktot - Ks + 1):Ktot else Vs <- numeric(0)
    Vv <- (1:Ktot)[- c(Vf, Va, Vs)]
    Vvu <- Vv[1:(Ku - Ko)]
    Vvc <- Vv[- (0:(Ku - Ko))]

    Kv <- length(Vv)

    Xf <- lapply(X, function(x) x[, fixedpar, drop = FALSE])
    Xc <- lapply(X, function(x) x[, correlated, drop = FALSE])
    Xu <- lapply(X, function(x) x[, utwopars, drop = FALSE])
    Xa <- lapply(X, function(x) x[, randompar, drop = FALSE])

    set.seed(seed)
    random.nb <- make.random.nb(R * ifelse(! is.null(id), n, N), Ku + Kc, halton)
    step <- 2
    repeat{
        step <- step / 2
        if (step < stptol) break
        betaf <- param[Vf] + step * direction[Vf]
        mua <- param[Va] + step * direction[Va]
        siga <- param[Vv] + step * direction[Vv]
        if (! is.null(Xs)){
            lambda <- param[Vs] + step * direction[Vs]
            sigma <- 1 + as.numeric(Xs %*% lambda)
            fpsigma <- Xs
        }
        else sigma <- rep(1, nrow(X[[1]]))
        if (length(betaf)){
            A <- lapply(Xf, function(x) as.vector(crossprod(t(as.matrix(x) / sigma), betaf)))
        }
        else{
            A <- lapply(Xf, function(x) rep(0, N))
        }
        B <- vector(mode = "list", length = J)
        for (j in 1:J) B[[j]] <- matrix(NA, N, R)
        ndraws <- ifelse(panel, n, N)
        NOUVEAU <- FALSE
        if (! NOUVEAU) b <- make.beta(mua, siga, rpar, random.nb, correlation)
        if (panel) theIds <- unique(id)
        for (i in 1:ndraws){
            if (panel){
                anid <- theIds[i]
                theRows <- which(id == anid)
            }
            else theRows <- i
            if (NOUVEAU){
                b <- make.beta(mua, siga, rpar,
                               random.nb[((i - 1) * R + 1):(i * R), ,
                                         drop = FALSE] , correlation)
                betaa <- b$betaa
            }
            else betaa <- b$betaa[((i - 1) * R + 1):(i * R), , drop = FALSE]

            for (j in 1:J){
                B[[j]][theRows,] <-
                               tcrossprod(Xa[[j]][theRows, ,
                                                  drop = FALSE] / sigma[theRows], betaa)
            }
        }
        linpred <- mapply(function(x, y) x + y, A, B, SIMPLIFY = FALSE)
        linpred <- Reduce("cbind", lapply(linpred, apply, 1, mean))
        AB <- mapply(function(x, y) exp(x + y), A, B, SIMPLIFY = FALSE)
        S <- suml(AB)
        P <- lapply(AB, function(x) x / S)
        probabilities <- sapply(P, function(x) apply(x, 1, mean))
        colnames(linpred) <- colnames(probabilities)
        Pch <- suml(mapply("*", P, y, SIMPLIFY = FALSE))
        if (panel) Pch <- apply(Pch, 2, tapply, id, prod)
        pm <- apply(Pch, 1, mean)
        if (panel) lnl <- opposite * sum(weights[! duplicated(id)] * log(pm))
        else lnl <- opposite * sum(weights * log(pm))

        Pch2 <- Pch / rowSums(Pch)
        Pch2 <- as.numeric(t(Pch2))
        if (! NOUVEAU){
            indparam <- Pch2 * b$betaa
            indparam <- apply(indparam, 2, tapply, rep(1:ndraws, each = R), sum)
            if(panel) indparam <- data.frame(id = unique(id), as.data.frame(indparam))
        }
        else indparam <- NULL
        
        if (is.null(initial.value) || lnl <= initial.value) break
    }
    if (gradient){

        Xf.ch <- suml(mapply("*", Xf, y, SIMPLIFY = FALSE))
        Xa.ch <- suml(mapply("*", Xa, y, SIMPLIFY = FALSE))
        if (Kc){
            vecX <- rep(1:Kc, 1:Kc)
            Xc <- lapply(Xc,  function(x) x[, vecX, drop = FALSE])
            Xc.ch <- suml(mapply("*", Xc, y, SIMPLIFY = FALSE))
            colnames(Xc.ch) <- names(param)[Vvc]
            # pb avec Vvc -> integer(0)
        }
        if (Ku - Ko){
            Xu.ch <- suml(mapply("*", Xu, y, SIMPLIFY = FALSE))
            colnames(Xu.ch) <- names(param)[Vvu]
        }
        if ((Ku - Ko) & Kc){
            Xas.ch <- cbind(Xu.ch, Xc.ch)
            Xas <- mapply("cbind", Xu, Xc, SIMPLIFY = FALSE)
        }
        if ((Ku - Ko) & ! Kc){
            Xas.ch <- Xu.ch
            Xas <- Xu
        }
        if (! (Ku - Ko) & Kc){
            Xas.ch <- Xc.ch
            Xas <- Xc
        }
        if (panel){
            Pch <- Pch[as.character(id), ]
            pm <- apply(Pch, 1, mean)
        }
        
        PCP <- lapply(P, function(x) Pch * x)
        PCPs <- lapply(PCP, function(x) apply(x, 1, sum))
        grad.cst <- Xf.ch - suml(mapply("*", Xf, PCPs, SIMPLIFY = FALSE)) / (R * pm)
        PCPm <- PCPs <- vector(mode = "list", length = J)
        Pchm <- matrix(NA, N, Ku + Kc)
        Pchs <- matrix(NA, N, Kv)
        for (j in 1:J){
            PCPm[[j]] <- matrix(NA, N, Ku + Kc)
            PCPs[[j]] <- matrix(NA, N, Kv)
        }
        for (i in 1:ndraws){
            if (panel){
                anid <- theIds[i]
                theRows <- which(id == anid)
            }
            else theRows <- i
            ## b <- make.beta(mua, siga, rpar,
            ##                random.nb[((i - 1) * R + 1):(i * R), , drop = FALSE],
            ##                correlation)
            index.i <- ((i - 1) * R + 1):(i * R)
            Pchm[theRows, ] <- tcrossprod(Pch[theRows, ], t(b$betaa.mu[index.i, , drop = FALSE]))
            if (! is.null(b$betaa.sigma))
                Pchs[theRows, ] <- tcrossprod(Pch[theRows, ], t(b$betaa.sigma[index.i, , drop = FALSE]))
            for (j in 1:J){
                PCPm[[j]][theRows,] <- tcrossprod(PCP[[j]][theRows, ], t(b$betaa.mu[index.i, , drop = FALSE]))
                if (! is.null(b$betaa.sigma))
                    PCPs[[j]][theRows,] <- tcrossprod(PCP[[j]][theRows, ], t(b$betaa.sigma[index.i, , drop = FALSE])) 
            }
        }
        grad.mu <- (Pchm * Xa.ch  - suml(mapply("*", Xa,  PCPm, SIMPLIFY = FALSE))) / (R * pm)
        if (! is.null(b$betaa.sigma)){
            grad.sd <- (Pchs * Xas.ch - suml(mapply("*", Xas, PCPs, SIMPLIFY = FALSE))) / (R * pm)        
            gradi <- matrix(NA, N, K + ncol(grad.sd))
            gradi[, Vf] <- - grad.cst
            gradi[, Va] <- - grad.mu
            gradi[, Vv] <- - grad.sd
        }
        else{
            gradi <- matrix(NA, N, K)
            gradi[, Vf] <- - grad.cst
            gradi[, Va] <- - grad.mu
        }
        if (! is.null(Xs)){
            betas <- numeric(length = length(c(Va, Vf, Vv)))
            betas[Vf] <- betaf
            betas[Va] <- mua
            betas[Vv] <- siga
            gradsi <- sapply(as.data.frame(fpsigma),
                             function(x) apply(x * t(t(gradi) * betas), 1, sum)) / (- sigma ^ 2)
            gradi <- gradi / sigma
            gradi <- cbind(gradi, gradsi)
        }

        if (! is.null(weights)) gradi <- weights * gradi
        colnames(gradi) <- names(param)
        attr(lnl, "gradi") <-  opposite * gradi
        attr(lnl, "gradient") <- apply(gradi, 2, sum)
        attr(lnl, "step") <- step
    }
    if (step < stptol) lnl <- NULL
    else{
        attr(lnl, "probabilities") <- probabilities
        attr(lnl, "fitted") <- pm
        attr(lnl, "step") <- step
        attr(lnl, "indparam") <- indparam
        attr(lnl, "linpred") <- linpred
    }
    lnl
}
