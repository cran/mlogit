## -----------------------------------------------------------------------------
#| echo: false
oopts <-options(width = 70)


## -----------------------------------------------------------------------------
#| message: false
library(mlogit)
data("Heating", package = "mlogit")
H <- dfidx(Heating, choice = "depvar", varying = c(3:12, 17:21))
m <- mlogit(depvar ~ ic + oc | 0, H)
gaze(m)


## -----------------------------------------------------------------------------
#| collapse: true
apply(fitted(m, outcome = FALSE), 2, mean)


## -----------------------------------------------------------------------------
#| collapse: true
coef(m)["oc"]/coef(m)["ic"]


## -----------------------------------------------------------------------------
#| message: false
#| collapse: true
H$lcc <- H$ic + H$oc / 0.12
mlcc <- mlogit(depvar ~ lcc | 0, H)
lrtest(m, mlcc) |> gaze()
qchisq(0.05, df = 1, lower.tail = FALSE)


## -----------------------------------------------------------------------------
mc <- mlogit(depvar ~ ic + oc, H, reflevel = 'hp')
gaze(mc)


## -----------------------------------------------------------------------------
#| collapse: true
apply(fitted(mc, outcome = FALSE), 2, mean)


## -----------------------------------------------------------------------------
#| collapse: true
wtp <- unname(coef(mc)["oc"] / coef(mc)["ic"])
wtp
r <- 1 / wtp
r


## -----------------------------------------------------------------------------
update(mc, reflevel = "gr") |> coef()


## -----------------------------------------------------------------------------
mi <- mlogit(depvar ~ oc + I(ic / income), H)
gaze(mi)


## -----------------------------------------------------------------------------
mi2 <- mlogit(depvar ~ oc + ic | income, H, reflevel = "hp")


## -----------------------------------------------------------------------------
#| collapse: true
lrtest(mc, mi2) |> gaze()
waldtest(mc, mi2) |> gaze()
scoretest(mc, mi2) |> gaze()


## -----------------------------------------------------------------------------
X <- model.matrix(mc)
alt <- idx(mc, 2)
chid <- idx(mc, 1)
eXb <- as.numeric(exp(X %*% coef(mc)))
SeXb <- tapply(eXb, chid, sum)
P <- eXb / SeXb[chid]
P <- matrix(P, ncol = 5, byrow = TRUE)
head(P)


## -----------------------------------------------------------------------------
#| collapse: true
apply(P, 2, mean)


## -----------------------------------------------------------------------------
#| collapse: true
apply(fitted(mc, outcome = FALSE), 2, mean)


## -----------------------------------------------------------------------------
#| collapse: true
Hn <- H
Hn$ic[Hn$id2 == "hp"] <- 0.9 * Hn$ic[Hn$id2 == "hp"]
apply(predict(mc, newdata = Hn, shape = "wide"), 2, mean)


## -----------------------------------------------------------------------------
#| collapse: true
X <- model.matrix(mc)
Xn <- X[idx(mc, 2) == "ec",]
Xn[, "ic"] <- Xn[, "ic"] + 200
Xn[, "oc"] <- Xn[, "oc"] * 0.75
unchid <- unique(idx(mc, 1))
rownames(Xn) <- paste(unchid, 'new', sep = ".")
chidb <- c(chid, unchid)
X <- rbind(X, Xn)
X <- X[order(chidb), ]
eXb <- as.numeric(exp(X %*% coef(mc)))
SeXb <- as.numeric(tapply(eXb, sort(chidb), sum))
P <- eXb / SeXb[sort(chidb)]
P <- matrix(P, ncol = 6, byrow = TRUE)
apply(P, 2, mean)


## -----------------------------------------------------------------------------
#| echo: false
options(oopts)

