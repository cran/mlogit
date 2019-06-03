## ----label = setup, include = FALSE----------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, widtht = 65)
options(width = 65)

## --------------------------------------------------------------
library("mlogit")
data("Heating", package = "mlogit")
H <- mlogit.data(Heating, shape = "wide", choice = "depvar", varying = c(3:12))
m <- mlogit(depvar ~ ic + oc | 0, H)
summary(m)

## --------------------------------------------------------------
apply(fitted(m, outcome = FALSE), 2, mean)

## --------------------------------------------------------------
coef(m)["oc"]/coef(m)["ic"]

## --------------------------------------------------------------
H$lcc <- H$ic + H$oc / 0.12
mlcc <- mlogit(depvar ~ lcc | 0, H)
library("lmtest")
lrtest(m, mlcc)
qchisq(0.05, df = 1, lower.tail = FALSE)

## --------------------------------------------------------------
mc <- mlogit(depvar ~ ic + oc, H, reflevel = 'hp')
summary(mc)
apply(fitted(mc, outcome = FALSE), 2, mean)

## --------------------------------------------------------------
wtp <- coef(mc)["oc"] / coef(mc)["ic"]
wtp
r <- 1 / wtp
r

## --------------------------------------------------------------
update(mc, reflevel = "gr")

## --------------------------------------------------------------
mi <- mlogit(depvar ~ oc + I(ic / income), H)
summary(mi)

## --------------------------------------------------------------
mi2 <- mlogit(depvar ~ oc + ic | income, H, reflevel = "hp")

## --------------------------------------------------------------
lrtest(mc, mi2)
waldtest(mc, mi2)
scoretest(mc, mi2)

## --------------------------------------------------------------
X <- model.matrix(mc)
alt <- index(H)$alt
chid <- index(H)$chid
eXb <- as.numeric(exp(X %*% coef(mc)))
SeXb <- tapply(eXb, chid, sum)
P <- eXb / SeXb[chid]
P <- matrix(P, ncol = 5, byrow = TRUE)
head(P)
apply(P, 2, mean)

## --------------------------------------------------------------
apply(fitted(mc, outcome = FALSE), 2, mean)

## --------------------------------------------------------------
Hn <- H
Hn[Hn$alt == "hp", "ic"] <- 0.9 * Hn[Hn$alt == "hp", "ic"]
apply(predict(mc, newdata = Hn), 2, mean)

## --------------------------------------------------------------
X <- model.matrix(mc)
Xn <- X[alt == "ec",]
Xn[, "ic"] <- Xn[, "ic"] + 200
Xn[, "oc"] <- Xn[, "oc"] * 0.75
unchid <- unique(index(H)$chid)
rownames(Xn) <- paste(unchid, 'new', sep = ".")
chidb <- c(chid, unchid)
X <- rbind(X, Xn)
X <- X[order(chidb), ]
eXb <- as.numeric(exp(X %*% coef(mc)))
SeXb <- as.numeric(tapply(eXb, sort(chidb), sum))
P <- eXb / SeXb[sort(chidb)]
P <- matrix(P, ncol = 6, byrow = TRUE)
apply(P, 2, mean)

