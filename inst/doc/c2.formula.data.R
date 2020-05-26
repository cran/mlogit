## ----label = setup, include = FALSE----------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, widtht = 65)
options(width = 65)

## ----label = 'loading mlogit'----------------------------------
library("mlogit")

## ----label = 'Train data'--------------------------------------
data("Train", package = "mlogit")
Train$choiceid <- 1:nrow(Train)
head(Train, 3)

## ----label = 'dfidx for Train'---------------------------------
Tr <- dfidx(Train, shape = "wide", varying = 4:11, sep = "_",
            idx = list(c("choiceid", "id")), idnames = c(NA, "alt"))

## ----label = 'data transformation for Train'-------------------
Tr$price <- Tr$price / 100 * 2.20371
Tr$time <- Tr$time / 60

## ----label = 'head of the transformed Train data set'----------
head(Tr, 3)

## ----label = 'index of the transformed Train data set'---------
head(idx(Tr), 3)

## ----label = 'loading ModeCanada'------------------------------
data("ModeCanada", package = "mlogit")
head(ModeCanada)

## ----label = 'applying dfidx to Modecanada (1)'----------------
MC <- dfidx(ModeCanada, subset = noalt == 4,
            alt.levels = c("train", "air", "bus", "car"))

## ----label = 'applying dfidx to Modecanada (2)'----------------
MC <- dfidx(ModeCanada, subset = noalt == 4, idx = list(NA, "alt"))

## ----label = 'applying dfidx to Modecanada (3)'----------------
MC <- dfidx(ModeCanada, subset = noalt == 4, idx = "case",
            alt.levels = c("train", "air", "bus", "car"))

## ----label = 'applying dfidx to Modecanada (4)'----------------
MC <- dfidx(ModeCanada, subset = noalt == 4, idx = c("case", "alt"))

## ----label = 'ModeCanada without idx'--------------------------
MC <- dfidx(ModeCanada, subset = noalt == 4)

## ----label = 'applying dfidx to Modecanada (5)'----------------
MC <- dfidx(ModeCanada, subset = noalt == 4, idx = c("case", "alt"),
            drop.index = FALSE)
head(MC)

## ----label = 'a three parts formula'---------------------------
library("Formula")
f <- Formula(choice ~ cost | income + urban | ivt)

## ----label = 'ommission of some parts (1)'---------------------
f2 <- Formula(choice ~ cost + ivt | income + urban)
f2 <- Formula(choice ~ cost + ivt | income + urban | 0)

## ----label = 'ommission of some parts (2)'---------------------
f3 <- Formula(choice ~ 0 | income | 0)
f3 <- Formula(choice ~ 0 | income)

## ----label = 'ommission of some parts (3)'---------------------
f4 <- Formula(choice ~ cost + ivt)
f4 <- Formula(choice ~ cost + ivt | 1)
f4 <- Formula(choice ~ cost + ivt | 1 | 0)

## ----label = 'removing the intercept'--------------------------
f5 <- Formula(choice ~ cost | income + 0 | ivt)
f5 <- Formula(choice ~ cost | income - 1 | ivt)

## ----label = 'model.matrix method for Formula objects'---------
f <- Formula(choice ~ cost | income  | ivt)
mf <- model.frame(MC, f)
head(model.matrix(mf), 4)

## ----label = 'convenient statpval function'--------------------
statpval <- function(x){
    if (inherits(x, "anova")) 
        result <- as.matrix(x)[2, c("Chisq", "Pr(>Chisq)")]
    if (inherits(x, "htest")) result <- c(x$statistic, x$p.value)
    names(result) <- c("stat", "p-value")
    round(result, 3)
}

