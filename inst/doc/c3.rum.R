## ----label = setup, include = FALSE----------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, widtht = 65)
options(width = 65)

## ----label = 'multinomial logit with a mlogit.data'------------
library("mlogit")
data("ModeCanada", package = "mlogit")
MC <- mlogit.data(ModeCanada, subset = noalt == 4, chid.var = "case",
                  alt.var = "alt", drop.index = TRUE)
ml.MC1 <- mlogit(choice ~ cost + freq + ovt | income | ivt, MC)

## ----label = 'multinomial logit with an ordinary data.frame'----
ml.MC1b <- mlogit(choice ~ cost + freq + ovt | income | ivt, ModeCanada,
subset = noalt == 4, alt.var = "alt", chid.var = "case")

## ----label = 'estimation on a subset of alternatives'----------
MC$time <- with(MC, ivt + ovt)
ml.MC1 <- mlogit(choice ~ cost + freq | income | time, MC, 
alt.subset = c("car", "train", "air"), reflevel = "car")

## ----label = 'summary method for mlogit'-----------------------
summary(ml.MC1)

## ----label = 'fitted method for mlogit'------------------------
head(fitted(ml.MC1, type = "outcome"))
head(fitted(ml.MC1, type = "probabilities"), 4)

## ----label = 'computation of the log likelihood and the market shares'----
sum(log(fitted(ml.MC1, type = "outcome")))
logLik(ml.MC1)
apply(fitted(ml.MC1, type = "probabilities"), 2, mean)

## ----label = 'default behaviour of the predict method'---------
predict(ml.MC1)

## ----label = 'predicting with different data'------------------
NMC <- MC
NMC[index(NMC)$alt == "train", "time"] <- 0.8 *
NMC[index(NMC)$alt == "train", "time"]
Oprob <- fitted(ml.MC1, type = "probabilities")
Nprob <- predict(ml.MC1, newdata = NMC)
rbind(old = apply(Oprob, 2, mean), new = apply(Nprob, 2, mean))

## ----label = 'illustration of the IIA property'----------------
head(Nprob[, "air"] / Nprob[, "car"])
head(Oprob[, "air"] / Oprob[, "car"])

## ----label = 'computation of the initital logsum'--------------
ivbefore <- logsum(ml.MC1)

## ----label = 'computation of the after change logsum'----------
ivafter <- logsum(ml.MC1, data = NMC)

## ----label = 'computation of consumers surplus'----------------
surplus <- - (ivafter - ivbefore) / coef(ml.MC1)["cost"]
summary(surplus)

## ----label = 'marginal effects for an individual specific covariate'----
effects(ml.MC1, covariate = "income", type = "ar")

## ----label = 'marginal effects for an alternative specific covariate'----
effects(ml.MC1, covariate = "cost", type = "rr")

## ----label = 'computation of the marginal rate of substitution'----
coef(ml.MC1)[grep("time", names(coef(ml.MC1)))] /
    coef(ml.MC1)["cost"] * 60 

## ----label = 'estimation of the multinomial logit model for the NOx data'----
data("NOx", package = "mlogit")
NOx$kdereg <- with(NOx, kcost * (env == "deregulated"))
NOxml <- mlogit.data(NOx, chid.var = "chid", alt.var = "alt",
id.var = "id")
ml.pub <- mlogit(choice ~ post + cm + lnb + vcost + kcost + kcost:age |
- 1, subset = available & env == "public", data = NOxml)
ml.reg <- update(ml.pub, subset = available & env == "regulated")
ml.dereg <- update(ml.pub, subset = available & env == "deregulated")
ml.pool <- mlogit(choice ~ post + cm + lnb + vcost + kcost + kcost:age +
kdereg | - 1 | 0 | env, subset = available == 1, data = NOxml,
method = "bhhh")

## ----label = 'results of the NOx estimations', results = 'asis'----
library("texreg")
htmlreg(list(Public = ml.pub, Deregulated = ml.dereg, Regulated = ml.reg,
             Pooled = ml.pool), caption = "Environmental compliance choices.",
        omit.coef = "(post)|(cm)|(lnb)", float.pos = "hbt", label = "tab:nox")

## ----label = 'likelihood ratio test for the NOx data'----------
stat <- 2 * (logLik(ml.dereg) + logLik(ml.reg) +
             logLik(ml.pub) - logLik(ml.pool))
stat
pchisq(stat, df = 9, lower.tail = FALSE)

