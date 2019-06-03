## ----label = setup, include = FALSE----------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, widtht = 65)
options(width = 65)

## ----label = 'multinomial logit for the Train data'------------
library("mlogit")
data("Train", package = "mlogit")
Tr <- mlogit.data(Train, shape = "wide", choice = "choice",
                  varying = 4:11, sep = "_", id.var = "id",
                  opposite = c("price", "comfort", "time", "change"))
Tr$price <- Tr$price / 100 * 2.20371
Tr$time <- Tr$time / 60
Train.ml <- mlogit(choice ~ price + time + change + comfort | - 1, Tr)
coef(summary(Train.ml))

## ----label = 'marginal rates of substitution for Train'--------
coef(Train.ml)[- 1] / coef(Train.ml)[1]

## ----label = 'mixed logit estimation for Train (1)'------------
Train.mxlu <- mlogit(choice ~ price + time + change + comfort | - 1, Tr,
panel = TRUE, rpar = c(time = "n", change = "n", comfort = "n"), R = 100,
correlation = FALSE, halton = NA, method = "bhhh")
names(coef(Train.mxlu))

## ----label = 'mixed logit estimation for Train (2)'------------
Train.mxlc <- update(Train.mxlu, correlation = TRUE)
names(coef(Train.mxlc))

## ----label = 'summary of a random parameter in the preference space'----
marg.ut.time <- rpar(Train.mxlc, "time")
summary(marg.ut.time)

## ----label = 'summary of a random parameter in the wtp space'----
wtp.time <- rpar(Train.mxlc, "time", norm = "price")
summary(wtp.time)

## ----label = 'statistics of the random parameter in the wtp space'----
mean(rpar(Train.mxlc, "time", norm = "price"))
med(rpar(Train.mxlc, "time", norm = "price"))
stdev(rpar(Train.mxlc, "time", norm = "price"))

## ----label = 'vcov method for mlogit objects'------------------
vcov(Train.mxlc, what = "rpar")
vcov(Train.mxlc, what = "rpar", type = "cor")
summary(vcov(Train.mxlc, what = "rpar", type = "cor"))
summary(vcov(Train.mxlc, what = "rpar", type = "cov"))

## ----label = 'specific methods for random parameters'----------
cor.mlogit(Train.mxlc)
cov.mlogit(Train.mxlc)
stdev(Train.mxlc)

## ----label = 'mixed logit with a subset of correlated paramaters'----
Train.mxlc2 <- update(Train.mxlc, correlation = c("time", "comfort"))
vcov(Train.mxlc2, what = "rpar", type = "cor")

## ----label = 'convenient statpval function 3', include = FALSE----
statpval <- function(x){
    if (inherits(x, "anova")) 
        result <- as.matrix(x)[2, c("Chisq", "Pr(>Chisq)")]
    if (inherits(x, "htest")) result <- c(x$statistic, x$p.value)
    names(result) <- c("stat", "p-value")
    round(result, 3)
}

## ----label = 'tests of no correlated random effects'-----------
lr.mxc <- lrtest(Train.mxlc, Train.ml)
wd.mxc <- waldtest(Train.mxlc)
library("car")
lh.mxc <- linearHypothesis(Train.mxlc, c("chol.time:time = 0",
"chol.time:change =  0", "chol.time:comfort = 0", "chol.change:change = 0",
"chol.change:comfort = 0", "chol.comfort:comfort = 0"))
sc.mxc <- scoretest(Train.ml, rpar = c(time = "n", change = "n",
comfort = "n"), R = 100, correlation = TRUE, halton = NA, panel = TRUE)
sapply(list(wald = wd.mxc, lh = lh.mxc, score = sc.mxc, lr = lr.mxc),
statpval)

## ----label = 'tests of no correlation'-------------------------
lr.corr <- lrtest(Train.mxlc, Train.mxlu)
lh.corr <- linearHypothesis(Train.mxlc, c("chol.time:change = 0",
"chol.time:comfort = 0", "chol.change:comfort = 0"))
wd.corr <- waldtest(Train.mxlc, correlation = FALSE)
sc.corr <- scoretest(Train.mxlu, correlation = TRUE)
sapply(list(wald = wd.corr, lh = lh.corr, score = sc.corr, lr = lr.corr),
statpval)

## ----label = 'multinomial model for RiskyTransport'------------
data("RiskyTransport", package = "mlogit")
RT <- mlogit.data(RiskyTransport, shape = "long", choice = "choice", 
chid.var = "chid", alt.var = "mode", id.var = "id")
ml.rt <- mlogit(choice ~ cost + risk  + seats + noise + crowdness +
convloc + clientele | 0, data = RT, weights = weight)

## ----label = 'coef of risk and cost'---------------------------
coef(ml.rt)[c("risk", "cost")]

## ----label = 'mixed effects model for RiskyTransport'----------
mx.rt <- mlogit(choice ~ cost + risk  + seats + noise + crowdness +
convloc + clientele | 0, data = RT, weights = weight,
rpar = c(cost = 'zbt', risk = 'zbt'), R = 100, halton = NA, panel = TRUE)

## ----label = 'results for RiskyTransport', results = 'asis'----
library("texreg")
htmlreg(list('Multinomial logit' = ml.rt, 'Mixed logit' = mx.rt),
digits = 3, float.pos = "hbt", label = "tab:risktr", single.row = TRUE,
caption = "Transportation choices.")

## ----label = 'individual parameters'---------------------------
indpar <- fitted(mx.rt, type = "parameters")
head(indpar)

## ----label = 'individal parameters'----------------------------
indpar$VSL <- with(indpar, risk / cost * 100)
quantile(indpar$VSL, c(0.025, 0.975))
mean(indpar$VSL)

## ----label = 'max VSL'-----------------------------------------
max(indpar$cost)
max(indpar$VSL)

