## ----label = setup, include = FALSE----------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, widtht = 65)
options(width = 65)

## --------------------------------------------------------------
library("mlogit")
data("Mode", package="mlogit")
Mo <- dfidx(Mode, choice = "choice", varying = 2:9)

## ----probit1---------------------------------------------------
p1 <- mlogit(choice ~ cost + time, Mo, seed = 20, 
             R = 100, probit = TRUE)

## --------------------------------------------------------------
summary(p1)

## --------------------------------------------------------------
L1 <- matrix(0, 3, 3)
L1[! upper.tri(L1)] <- c(1, coef(p1)[6:10])

## --------------------------------------------------------------
L1 %*% t(L1)

## ----probit2---------------------------------------------------
p2 <- mlogit(choice ~ cost + time, Mo, seed = 21, 
             R = 100, probit = TRUE)

## --------------------------------------------------------------
coef(p2)

## --------------------------------------------------------------
actShares <- tapply(Mo$choice, Mo$id2, mean)

## --------------------------------------------------------------
predShares <- apply(fitted(p1, outcome = FALSE), 2, mean)
rbind(predShares, actShares)
sum(predShares)

## --------------------------------------------------------------
Mo2 <- Mo
Mo2[idx(Mo2, 2) == 'car', 'cost'] <- Mo2[idx(Mo2, 2) == 'car', 'cost'] * 2
newShares <- apply(predict(p1, newdata = Mo2), 2, mean)
cbind(original = actShares, new = newShares, 
      change = round((newShares - actShares) / actShares * 100))

