## -----------------------------------------------------------------------------
#| echo: false
oopts <-options(width = 70)


## -----------------------------------------------------------------------------
#| label: "multinomial logit with a dfidx"
#| message: false
library(mlogit)
MC <- dfidx(ModeCanada, subset = noalt == 4)
ml.MC1 <- mlogit(choice ~ cost + freq + ovt | income | ivt, MC)


## -----------------------------------------------------------------------------
#| label: "multinomial logit with an ordinary data.frame"
ml.MC1b <- mlogit(choice ~ cost + freq + ovt | income | ivt, ModeCanada,
                  subset = noalt == 4, idx = c("case", "alt"))


## -----------------------------------------------------------------------------
#| label: "estimation on a subset of alternatives"
MC$time <- with(MC, ivt + ovt)
ml.MC1 <- mlogit(choice ~ cost + freq | income | time, MC, 
                 alt.subset = c("car", "train", "air"),
                 reflevel = "car")


## -----------------------------------------------------------------------------
#| label: "summary method for mlogit"
summary(ml.MC1)


## -----------------------------------------------------------------------------
#| label: "fitted method for mlogit"
head(fitted(ml.MC1, type = "outcome"))
head(fitted(ml.MC1, type = "probabilities"), 4)


## -----------------------------------------------------------------------------
#| label: "computation of the log likelihood and the market shares"
#| collapse: true
sum(log(fitted(ml.MC1, type = "outcome")))
logLik(ml.MC1)
apply(fitted(ml.MC1, type = "probabilities"), 2, mean)


## -----------------------------------------------------------------------------
#| label: "default behaviour of the predict method"
#| collapse: true
head(predict(ml.MC1, shape = "wide"))


## -----------------------------------------------------------------------------
#| label: "predicting with different data"
#| message: false
#| warning: false
NMC <- MC
NMC$time[NMC$alt == "train"] <- 0.8 * NMC$time[NMC$alt == "train"]
Oprob <- fitted(ml.MC1, type = "probabilities")
Nprob <- predict(ml.MC1, newdata = NMC, shape = "wide")
rbind(old = apply(Oprob, 2, mean), new = apply(Nprob, 2, mean))


## -----------------------------------------------------------------------------
#| label: "illustration of the IIA property"
#| collapse: true
head(Nprob[, "air"] / Nprob[, "car"])
head(Oprob[, "air"] / Oprob[, "car"])


## -----------------------------------------------------------------------------
#| label: "computation of the initital logsum"
ivbefore <- logsum(ml.MC1)


## -----------------------------------------------------------------------------
#| label: "computation of the after change logsum"
ivafter <- logsum(ml.MC1, data = NMC)


## -----------------------------------------------------------------------------
#| label: "computation of consumers surplus"
#| collapse: true
surplus <- - (ivafter - ivbefore) / coef(ml.MC1)["cost"]
summary(surplus)


## -----------------------------------------------------------------------------
#| label: "marginal effects for an individual specific covariate"
#| collapse: true
effects(ml.MC1, covariate = "income", type = "ar")


## -----------------------------------------------------------------------------
#| label: "marginal effects for an alternative specific covariate"
#| collapse: true
effects(ml.MC1, covariate = "cost", type = "rr")


## -----------------------------------------------------------------------------
#| label: "computation of the marginal rate of substitution"
#| collapse: true
coef(ml.MC1)[grep("time", names(coef(ml.MC1)))] /
    coef(ml.MC1)["cost"] * 60 


## -----------------------------------------------------------------------------
#| label: "estimation of the multinomial logit model for the NOx data"
NOx$kdereg <- with(NOx, kcost * (env == "deregulated"))
NOxml <- dfidx(NOx, idx = c(id = "chid", "alt"))
ml.pub <- mlogit(choice ~ post + cm + lnb + vcost +
                     kcost + kcost:age | - 1,
                 subset = available & env == "public",
                 data = NOxml)
ml.reg <- update(ml.pub, subset = available & env == "regulated")
ml.dereg <- update(ml.pub, subset = available & env == "deregulated")
ml.pool <- ml.dereg
ml.pool <- mlogit(choice ~ post + cm + lnb + vcost + kcost +
                      kcost:age + kdereg | - 1 | 0 | env,
                  subset = available == 1, data = NOxml,
                  method = "bhhh")


## -----------------------------------------------------------------------------
#| label: models
#| echo: false
models <- list(Public = ml.pub,
     Deregulated = ml.dereg,
     Reguated = ml.reg,
     Pooled = ml.pool)
trms <- unique(Reduce("c", lapply(models, function(x) names(coef(x)))))
z <- Reduce("cbind", lapply(models, function(x) coef(x)[trms]))
dimnames(z) <- list(trms, names(models))
print(z, digits = 3, na.print = "-")


## -----------------------------------------------------------------------------
#| label: "likelihood ratio test for the NOx data"
#| collapse: true
stat <- 2 * (logLik(ml.dereg) + logLik(ml.reg) +
             logLik(ml.pub) - logLik(ml.pool))
stat
pchisq(stat, df = 9, lower.tail = FALSE)


## -----------------------------------------------------------------------------
ids <- unique(ModeCanada$case)
set.seed(1L)
ids <- sample(ids, 50)
smallMC <- subset(ModeCanada, case %in% ids)
ml <- mlogit(choice ~ cost | income, smallMC, alt.subset = c("car", "train", "air"))


## -----------------------------------------------------------------------------
.p <- preds(ml)
.s <- slps(ml)
.p
head(.s)


## -----------------------------------------------------------------------------
summary(.p)
summary(.s)


## -----------------------------------------------------------------------------
#| echo: false
options(oopts)

