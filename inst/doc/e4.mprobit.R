## -----------------------------------------------------------------------------
#| echo: false
oopts <-options(width = 70)


## -----------------------------------------------------------------------------
library(mlogit)
data("Mode", package="mlogit")
Mo <- dfidx(Mode, choice = "choice", varying = 2:9)


## -----------------------------------------------------------------------------
#| label: probit1
p1 <- mlogit(choice ~ cost + time, Mo, seed = 20, 
             R = 100, probit = TRUE)


## -----------------------------------------------------------------------------
gaze(p1)


## -----------------------------------------------------------------------------
L1 <- matrix(0, 3, 3)
L1[! upper.tri(L1)] <- c(1, coef(p1)[6:10])


## -----------------------------------------------------------------------------
L1 %*% t(L1)


## -----------------------------------------------------------------------------
#| label: probit2
p2 <- mlogit(choice ~ cost + time, Mo, seed = 21, 
             R = 100, probit = TRUE)


## -----------------------------------------------------------------------------
coef(p2)


## -----------------------------------------------------------------------------
actShares <- tapply(Mo$choice, Mo$id2, mean)


## -----------------------------------------------------------------------------
#| collapse: true
predShares <- apply(fitted(p1, outcome = FALSE), 2, mean)
rbind(predShares, actShares)
sum(predShares)


## -----------------------------------------------------------------------------
Mo2 <- Mo
Mo2$cost[Mo2$id2 == 'car'] <- Mo2$cost[Mo2$id2 == 'car'] * 2
newShares <- apply(predict(p1, newdata = Mo2, shape = "wide"), 2, mean)
cbind(original = actShares, new = newShares, 
      change = round((newShares - actShares) / actShares * 100))


## -----------------------------------------------------------------------------
#| echo: false
options(oopts)

