## -----------------------------------------------------------------------------
#| echo: false
oopts <-options(width = 70)


## -----------------------------------------------------------------------------
#| label: "multinomial logit for the Train data"
#| message: false
library(mlogit)
data("Train", package = "mlogit")
Train$choiceid <- 1:nrow(Train)
Tr <- dfidx(Train, choice = "choice", varying = 4:11, sep = "_",
            opposite = c("price", "comfort", "time", "change"),
            idx = list(c("choiceid", "id")), idnames = c("chid", "alt"))
Tr$price <- Tr$price / 100 / 2.20371
Tr$time <- Tr$time / 60
Train.ml <- mlogit(choice ~ price + time + change + comfort | - 1, Tr)
Train.ml |> gaze()


## -----------------------------------------------------------------------------
#| label: "marginal rates of substitution for Train"
#| collapse: true
coef(Train.ml)[- 1] / coef(Train.ml)[1]


## -----------------------------------------------------------------------------
#| label: "mixed logit estimation for Train (1)"
Train.mxlu <- mlogit(choice ~ price + time + change + comfort | - 1, Tr,
panel = TRUE, rpar = c(time = "n", change = "n", comfort = "n"), R = 100,
correlation = FALSE, halton = NA, method = "bhhh")
names(coef(Train.mxlu))


## -----------------------------------------------------------------------------
#| label: "mixed logit estimation for Train (2)"
Train.mxlc <- update(Train.mxlu, correlation = TRUE)
names(coef(Train.mxlc))


## -----------------------------------------------------------------------------
#| label: "summary of a random parameter in the preference space"
#| collapse: true
marg.ut.time <- rpar(Train.mxlc, "time")
summary(marg.ut.time)


## -----------------------------------------------------------------------------
#| label: "summary of a random parameter in the wtp space"
#| collapse: true
wtp.time <- rpar(Train.mxlc, "time", norm = "price")
summary(wtp.time)


## -----------------------------------------------------------------------------
#| label: "statistics of the random parameter in the wtp space"
#| collapse: true
mean(rpar(Train.mxlc, "time", norm = "price"))
med(rpar(Train.mxlc, "time", norm = "price"))
stdev(rpar(Train.mxlc, "time", norm = "price"))


## -----------------------------------------------------------------------------
#| label: "vcov method for mlogit objects"
vcov(Train.mxlc, what = "rpar")
vcov(Train.mxlc, what = "rpar", type = "cor")
summary(vcov(Train.mxlc, what = "rpar", type = "cor"))
summary(vcov(Train.mxlc, what = "rpar", type = "cov"))


## -----------------------------------------------------------------------------
#| label: "specific methods for random parameters"
cor.mlogit(Train.mxlc)
cov.mlogit(Train.mxlc)
stdev(Train.mxlc)


## -----------------------------------------------------------------------------
#| label: "mixed logit with a subset of correlated paramaters"
Train.mxlc2 <- update(Train.mxlc, correlation = c("time", "comfort"))
vcov(Train.mxlc2, what = "rpar", type = "cor")


## -----------------------------------------------------------------------------
#| label: "tests of no correlated random effects"
#| collapse: true
lrtest(Train.mxlc, Train.ml) |> gaze()
waldtest(Train.mxlc) |> gaze()
car::lht(Train.mxlc,
         c("chol.time:time = 0", "chol.time:change =  0",
           "chol.time:comfort = 0", "chol.change:change = 0",
           "chol.change:comfort = 0", "chol.comfort:comfort = 0")) |>
    gaze()
scoretest(Train.ml,
          rpar = c(time = "n", change = "n", comfort = "n"),
          R = 100, correlation = TRUE, halton = NA, panel = TRUE)


## -----------------------------------------------------------------------------
#| label: "tests of no correlation"
#| collapse: true
lrtest(Train.mxlc, Train.mxlu) |> gaze()
car::lht(Train.mxlc,
         c("chol.time:change = 0", "chol.time:comfort = 0",
           "chol.change:comfort = 0")) |> gaze()
waldtest(Train.mxlc, correlation = FALSE) |> gaze()
scoretest(Train.mxlu, correlation = TRUE) |> gaze()


## -----------------------------------------------------------------------------
#| label: "multinomial model for RiskyTransport"
RT <- dfidx(RiskyTransport, choice = "choice",
            idx = list(c("chid", "id"), "mode"),
            idnames = c("chid", "alt"))
ml.rt <- mlogit(choice ~ cost + risk  + seats + noise + crowdness +
                convloc + clientele | 0, data = RT, weights = weight)


## -----------------------------------------------------------------------------
#| label: "coef of risk and cost"
#| collapse: true
coef(ml.rt)[c("risk", "cost")]


## -----------------------------------------------------------------------------
#| label: "mixed effects model for RiskyTransport"
mx.rt <- mlogit(choice ~ cost + risk  + seats + noise + crowdness +
                convloc + clientele | 0, data = RT, weights = weight,
                rpar = c(cost = 'zbt', risk = 'zbt'), R = 100,
                halton = NA, panel = TRUE)


## -----------------------------------------------------------------------------
#| label: risktr
#| echo: false
models <- list('Multinomial logit' = ml.rt,
               'Mixed logit' = mx.rt)
trms <- unique(Reduce("c", lapply(models, function(x) names(coef(x)))))
z <- Reduce("cbind", lapply(models, function(x) coef(x)[trms]))
dimnames(z) <- list(trms, names(models))
print(z, digits = 3, na.print = "-")


## -----------------------------------------------------------------------------
#| label: "individual parameters"
indpar <- fitted(mx.rt, type = "parameters")
head(indpar)


## -----------------------------------------------------------------------------
#| label: "individal parameters"
#| collapse: true
indpar$VSL <- with(indpar, risk / cost * 100)
quantile(indpar$VSL, c(0.025, 0.975))
mean(indpar$VSL)


## -----------------------------------------------------------------------------
#| label: "max VSL"
#| collapse: true
max(indpar$cost)
max(indpar$VSL)


## -----------------------------------------------------------------------------
#| echo: false
options(oopts)

