## ----label = setup, include = FALSE----------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, widtht = 65)
options(width = 65)

## --------------------------------------------------------------
library("mlogit")
data("Mode", package="mlogit")
Mo <- mlogit.data(Mode, choice = 'choice', shape = 'wide', 
                  varying = c(2:9))

## ----probit1, echo = FALSE, results = 'hide'-------------------
strt <- c(1.83086600, -1.28168186,  0.30935104, -0.41344010, -0.04665517,  1, 0.25997237,
          0.73648694,  1.30789474, -0.79818416,  0.43013035)
p1 <- mlogit(choice ~ cost + time, Mo, seed = 20, 
             R = 100, probit = TRUE, start = strt)

## ----eval = FALSE----------------------------------------------
#  p1 <- mlogit(choice ~ cost + time, Mo, seed = 20,
#               R = 100, probit = TRUE)

## --------------------------------------------------------------
summary(p1)

## --------------------------------------------------------------
L1 <- matrix(0, 3, 3)
L1[!upper.tri(L1)] <- c(1, coef(p1)[6:10])

## --------------------------------------------------------------
L1 %*% t(L1)

## ----probit2, echo = FALSE, results = 'hide'-------------------
strt <-  c(1.87149948, -1.28893595,  0.31455318, -0.43068703, -0.04752315,  1, 0.22888163,
           0.69781113,  1.33071717, -0.56802431,  0.71060138)
p2 <- mlogit(choice ~ cost + time, Mo, seed = 21, 
             R = 100, probit = TRUE, start = strt)

## ----eval = FALSE----------------------------------------------
#  p2 <- mlogit(choice ~ cost + time, Mo, seed = 21,
#               R = 100, probit = TRUE)

## --------------------------------------------------------------
coef(p2)

## --------------------------------------------------------------
actShares <- with(Mo, tapply(choice, alt, mean))

## --------------------------------------------------------------
predShares <- apply(fitted(p1, outcome = FALSE), 2, mean)
predShares
sum(predShares)

## --------------------------------------------------------------
Mo2 <- Mo
Mo2[Mo2$alt == 'car', 'cost'] <- Mo2[Mo2$alt == 'car', 'cost'] * 2
newShares <- apply(predict(p1, newdata = Mo2), 2, mean)
cbind(original = actShares, new = newShares, 
      change = round((newShares - actShares) / actShares * 100))

