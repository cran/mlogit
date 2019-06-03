## ----label = setup, include = FALSE----------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, widtht = 65)
options(width = 65)

## ----label = 'heteroscedastic model for the ModeCanada data'----
library("mlogit")
data("ModeCanada", package = "mlogit")
MC <- mlogit.data(ModeCanada, subset = noalt == 4, chid.var = "case",
                  alt.var = "alt", drop.index = TRUE)
ml.MC <- mlogit(choice ~ freq + cost + ivt + ovt | urban + income, MC, 
reflevel = 'car', alt.subset = c("car", "train", "air"))
hl.MC <- mlogit(choice ~ freq + cost + ivt + ovt | urban + income, MC, 
reflevel = 'car', alt.subset = c("car", "train", "air"), heterosc = TRUE)
coef(summary(hl.MC))[11:12, ]

## ----label = 'homoscedasticity tests: lr and Wald (1)'---------
lr.heter <- lrtest(hl.MC, ml.MC)
wd.heter <- waldtest(hl.MC, heterosc = FALSE)

## ----label = 'homoscedasticity tests: lr and Wald (2)', results = 'hide'----
lrtest(hl.MC)
waldtest(hl.MC)

## ----label = 'homoscedasticity tests: Wald test'---------------
library("car")
lh.heter <- linearHypothesis(hl.MC, c('sp.air = 1', 'sp.train = 1'))

## ----label = 'homoscedasticity tests: score test'--------------
sc.heter <- scoretest(ml.MC, heterosc = TRUE)

## ----label = 'convenient statpval function 2', include = FALSE----
statpval <- function(x){
    if (inherits(x, "anova")) 
        result <- as.matrix(x)[2, c("Chisq", "Pr(>Chisq)")]
    if (inherits(x, "htest")) result <- c(x$statistic, x$p.value)
    names(result) <- c("stat", "p-value")
    round(result, 3)
}

## ----label = 'homoscedasticity tests: results'-----------------
sapply(list(wald = wd.heter, lh = lh.heter, score = sc.heter,
lr = lr.heter), statpval)

## ----label = 'loading the JapaneseFDI data set'----------------
data("JapaneseFDI", package = "mlogit")
jfdi <- mlogit.data(JapaneseFDI, chid.var = "firm", alt.var = "region",
group.var = "country")

## ----label = 'multinomial logit for JapaneseFDI'---------------
ml.fdi <- mlogit(choice ~ log(wage) + unemp + elig + log(area) + scrate +
ctaxrate | 0, data = jfdi)

## ----label = 'lower model estimation'--------------------------
lm.fdi <- mlogit(choice ~ log(wage) + unemp + elig + log(area) | 0,
data = jfdi, subset = country == choice.c & ! country %in% c("PT", "IE"))

## ----label = 'use of the logsum function'----------------------
lmformula <- formula(lm.fdi)
head(logsum(ml.fdi, data = jfdi, formula = lmformula, type = "group"), 2)
head(logsum(ml.fdi, data = jfdi, formula = lmformula, type = "global"))
head(logsum(ml.fdi, data = jfdi, formula = lmformula, output = "obs"))
head(logsum(ml.fdi, data = jfdi, formula = lmformula, type = "global",
output = "obs"))

## ----label = 'adding the logsum to the data'-------------------
JapaneseFDI$iv <- logsum(lm.fdi, data = jfdi, formula = lmformula,
output = "obs")

## ----label = 'data suitable for the upper model'---------------
JapaneseFDI.c <- subset(JapaneseFDI, 
select = c("firm", "country", "choice.c", "scrate", "ctaxrate", "iv"))
JapaneseFDI.c <- unique(JapaneseFDI.c)
JapaneseFDI.c$choice.c <- with(JapaneseFDI.c, choice.c == country)

## ----label = 'estimation of the upper model'-------------------
jfdi.c <- mlogit.data(JapaneseFDI.c, choice = "choice.c",
alt.var = "country", chid.var = "firm", shape = "long")
um.fdi <- mlogit(choice.c ~ scrate + ctaxrate + iv | 0, data = jfdi.c)

## ----label = 'upper model with different iv coefficients'------
um2.fdi <- mlogit(choice.c ~ scrate + ctaxrate | 0 | iv, data = jfdi.c, 
constPar = c("iv:PT" = 1, "iv:IE" = 1))

## ----label = 'nested logit models'-----------------------------
nl.fdi <- mlogit(choice ~ log(wage) + unemp + elig + log(area) + scrate +
ctaxrate | 0, data = jfdi, nests = TRUE, un.nest.el = TRUE)
nl2.fdi <- update(nl.fdi, un.nest.el = FALSE, constPar = c('iv:PT' = 1,
'iv:IE' = 1))

## ----label = 'results for the JapaneseFDI data', results = 'asis'----
library("texreg")
htmlreg(list('Mult. logit' = ml.fdi, 'Lower model' = lm.fdi,
             'Upper model' = um.fdi, 'Upper model' = um2.fdi, 'Nested logit' = nl.fdi,
             'Nested logit' = nl2.fdi),
        fontsize = "footnotesize", float.pos = "hbt",  label = "tab:nlogit",
        caption = "Choice by Japanese firms of a european region.")

## ----label = 'test of no nests'--------------------------------
lr.nest <- lrtest(nl2.fdi)
wd.nest <- waldtest(nl2.fdi)
sc.nest <- scoretest(ml.fdi, nests = TRUE, constPar = c('iv:PT' = 1,
'iv:IE' = 1))

## ----label = 'test of no nests with linhyp'--------------------
lh.nest <- linearHypothesis(nl2.fdi, c("iv:BE = 1", "iv:DE = 1",
"iv:ES = 1", "iv:FR = 1", "iv:IT = 1", "iv:NL = 1", "iv:UK = 1"))

## ----label = 'results of the tests of no nests'----------------
sapply(list(wald = wd.nest, lh = lh.nest, score = sc.nest, lr = lr.nest),
statpval)

## ----label = 'computing the test for equal iv coefficients'----
lr.unest <- lrtest(nl2.fdi, nl.fdi)
wd.unest <- waldtest(nl2.fdi, un.nest.el = TRUE)
sc.unest <- scoretest(ml.fdi, nests = TRUE, un.nest.el = FALSE,
constPar = c('iv:IE' = 1, 'iv:PT' = 1))
lh.unest <- linearHypothesis(nl2.fdi, c("iv:BE = iv:DE", "iv:BE = iv:ES", 
"iv:BE = iv:FR", "iv:BE = iv:IT", "iv:BE = iv:NL", "iv:BE = iv:UK"))

## ----label = 'results of the tests of equal iv coefficients'----
sapply(list(wald = wd.unest, lh = lh.unest, score = sc.unest,
lr = lr.unest), statpval)

