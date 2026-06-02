## -----------------------------------------------------------------------------
#| echo: false
oopts <-options(width = 80)


## -----------------------------------------------------------------------------
#| label: "heteroscedastic model for the ModeCanada data"
#| message: false
library(mlogit)
MC <- dfidx(ModeCanada, subset = noalt == 4)
ml.MC <- mlogit(choice ~ freq + cost + ivt + ovt | urban + income, MC, 
                reflevel = 'car', alt.subset = c("car", "train", "air"))
hl.MC <- mlogit(choice ~ freq + cost + ivt + ovt | urban + income, MC, 
                reflevel = 'car', alt.subset = c("car", "train", "air"),
                heterosc = TRUE, hessian = FALSE)
coef(summary(hl.MC))[11:12, ]


## -----------------------------------------------------------------------------
#| label: homoscedasticity tests lr and Wald (1)"
#| eval: false
# lrtest(hl.MC, ml.MC)
# waldtest(hl.MC, heterosc = FALSE)


## -----------------------------------------------------------------------------
#| label: "homoscedasticity tests lr and Wald (2)"
#| collapse: true
lrtest(hl.MC) |> gaze()
waldtest(hl.MC) |> gaze()


## -----------------------------------------------------------------------------
#| label: "homoscedasticity tests: Wald test"
#| collapse: true
car::lht(hl.MC, c('sp.air = 1', 'sp.train = 1')) |> gaze()


## -----------------------------------------------------------------------------
#| label: "homoscedasticity tests: score test"
#| collapse: true
scoretest(ml.MC, heterosc = TRUE) |> gaze()


## -----------------------------------------------------------------------------
#| label: "loading the JapaneseFDI data set"
jfdi <- dfidx(JapaneseFDI, idx = c("firm", country = "region"),
              drop.index = FALSE)


## -----------------------------------------------------------------------------
#| label: "multinomial logit for JapaneseFDI"
ml.fdi <- mlogit(choice ~ log(wage) + unemp + elig + log(area) +
                     scrate + ctaxrate | 0, data = jfdi)


## -----------------------------------------------------------------------------
#| label: "lower model estimation"
lm.fdi <- mlogit(choice ~ log(wage) + unemp + elig + log(area) | 0,
                 data = jfdi,
                 subset = country == choice.c &
                     ! country %in% c("PT", "IE"))


## -----------------------------------------------------------------------------
#| label: "use of the logsum function"
#| collapse: true
lmformula <- formula(lm.fdi)
head(logsum(ml.fdi, data = jfdi, formula = lmformula, type = "group"), 2)
head(logsum(ml.fdi, data = jfdi, formula = lmformula, type = "global"))
head(logsum(ml.fdi, data = jfdi, formula = lmformula, output = "obs"))
head(logsum(ml.fdi, data = jfdi, formula = lmformula, type = "global",
            output = "obs"))


## -----------------------------------------------------------------------------
#| label: "adding the logsum to the data"
JapaneseFDI$iv <- logsum(lm.fdi, data = jfdi, formula = lmformula,
                         output = "obs")


## -----------------------------------------------------------------------------
#| label: "data suitable for the upper model"
JapaneseFDI.c <- subset(JapaneseFDI,
                        select = c("firm", "country", "choice.c",
                                   "scrate", "ctaxrate", "iv"))
JapaneseFDI.c <- unique(JapaneseFDI.c)
JapaneseFDI.c$choice.c <- with(JapaneseFDI.c, choice.c == country)


## -----------------------------------------------------------------------------
#| label: "estimation of the upper model"
jfdi.c <- dfidx(JapaneseFDI.c, choice = "choice.c",
                idnames = c("chid", "alt"))
um.fdi <- mlogit(choice.c ~ scrate + ctaxrate + iv | 0, data = jfdi.c)


## -----------------------------------------------------------------------------
#| label: "upper model with different iv coefficients"
um2.fdi <- mlogit(choice.c ~ scrate + ctaxrate | 0 | iv, data = jfdi.c, 
                  constPar = c("iv:PT" = 1, "iv:IE" = 1))


## -----------------------------------------------------------------------------
#| label: "nested logit models"
nl.fdi <- mlogit(choice ~ log(wage) + unemp + elig + log(area) +
                     scrate + ctaxrate | 0, data = jfdi,
                 nests = TRUE, un.nest.el = TRUE)
nl2.fdi <- update(nl.fdi, un.nest.el = FALSE,
                  constPar = c('iv:PT' = 1, 'iv:IE' = 1))


## -----------------------------------------------------------------------------
#| label: nlogit
#| echo: false
models <- list('Mult. logit' = ml.fdi,
               'Upper model' = um.fdi,
               'Upper model' = um2.fdi,
               'Nested logit' = nl.fdi,
               'Nested logit' = nl2.fdi)
trms <- unique(Reduce("c", lapply(models, function(x) names(coef(x)))))
z <- Reduce("cbind", lapply(models, function(x) coef(x)[trms]))
dimnames(z) <- list(trms, names(models))
print(z, digits = 3, na.print = "-")


## -----------------------------------------------------------------------------
#| label: "test of no nests"
#| collapse: true
lrtest(nl2.fdi) |> gaze()
waldtest(nl2.fdi) |> gaze()
scoretest(ml.fdi, nests = TRUE,
          constPar = c('iv:PT' = 1, 'iv:IE' = 1))  |>
    gaze()


## -----------------------------------------------------------------------------
#| label: "test of no nests with linhyp"
#| collapse: true
car::lht(nl2.fdi, c("iv:BE = 1", "iv:DE = 1", "iv:ES = 1", "iv:FR = 1",
                    "iv:IT = 1", "iv:NL = 1", "iv:UK = 1")) |> gaze()


## -----------------------------------------------------------------------------
#| label: "computing the test for equal iv coefficients"
#| collapse: true
lrtest(nl2.fdi, nl.fdi) |> gaze()
waldtest(nl2.fdi, un.nest.el = TRUE) |> gaze()
scoretest(ml.fdi, nests = TRUE, un.nest.el = FALSE,
          constPar = c('iv:PT' = 1, 'iv:IE' = 1)) |> gaze()
car::lht(nl2.fdi, c("iv:BE = iv:DE", "iv:BE = iv:ES",
                    "iv:BE = iv:FR", "iv:BE = iv:IT",
                    "iv:BE = iv:NL", "iv:BE = iv:UK")) |>
    gaze()


## -----------------------------------------------------------------------------
#| echo: false
options(oopts)

