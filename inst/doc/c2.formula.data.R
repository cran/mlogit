## ----label = setup, include = FALSE----------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, widtht = 65)
options(width = 65)

## ----label = 'loading mlogit'----------------------------------
library("mlogit")

## ----label = 'Train data'--------------------------------------
data("Train", package = "mlogit")
head(Train, 3)

## ----label = 'mlogit.data for Train'---------------------------
Tr <- mlogit.data(Train, shape = "wide", choice = "choice",
                  varying = 4:11, sep = "_", id.var = "id",
                  opposite = c("price", "comfort", "time", "change"))

## ----label = 'data transformation for Train'-------------------
Tr$price <- Tr$price / 100 * 2.20371
Tr$time <- Tr$time / 60

## ----label = 'head of the transformed Train data set'----------
head(Tr, 3)

## ----label = 'index of the transformed Train data set'---------
head(index(Tr), 3)

## ----label = 'loading ModeCanada'------------------------------
data("ModeCanada", package = "mlogit")
head(ModeCanada)

## ----label = 'applying mlogit.data to Modecanada (1)'----------
MC <- mlogit.data(ModeCanada, subset = noalt == 4,
                  alt.levels = c("train", "air", "bus", "car"))

## ----label = 'applying mlogit.data to Modecanada (2)'----------
MC <- mlogit.data(ModeCanada , subset = noalt == 4, alt.var = "alt")

## ----label = 'applying mlogit.data to Modecanada (3)'----------
MC <- mlogit.data(ModeCanada, subset = noalt == 4, chid.var = "case",
                  alt.levels = c("train", "air", "bus", "car"))

## ----label = 'applying mlogit.data to Modecanada (4)'----------
MC <- mlogit.data(ModeCanada, subset = noalt == 4, chid.var = "case",
                  alt.var = "alt")

## ----label = 'applying mlogit.data to Modecanada (5)'----------
MC <- mlogit.data(ModeCanada, subset = noalt == 4, chid.var = "case",
                  alt.var = "alt", drop.index = TRUE)
head(MC)

## ----label = 'a three parts formula'---------------------------
f <- mFormula(choice ~ cost | income + urban | ivt)

## ----label = 'ommission of some parts (1)'---------------------
f2 <- mFormula(choice ~ cost + ivt | income + urban)
f2 <- mFormula(choice ~ cost + ivt | income + urban | 0)

## ----label = 'ommission of some parts (2)'---------------------
f3 <- mFormula(choice ~ 0 | income | 0)
f3 <- mFormula(choice ~ 0 | income)

## ----label = 'ommission of some parts (3)'---------------------
f4 <- mFormula(choice ~ cost + ivt)
f4 <- mFormula(choice ~ cost + ivt | 1)
f4 <- mFormula(choice ~ cost + ivt | 1 | 0)

## ----label = 'removing the intercept'--------------------------
f5 <- mFormula(choice ~ cost | income + 0 | ivt)
f5 <- mFormula(choice ~ cost | income - 1 | ivt)

## ----label = 'model.matrix method for Formula objects'---------
f <- mFormula(choice ~ cost | income  | ivt)
head(model.matrix(f, MC), 4)

## ----label = 'convenient statpval function'--------------------
statpval <- function(x){
    if (inherits(x, "anova")) 
        result <- as.matrix(x)[2, c("Chisq", "Pr(>Chisq)")]
    if (inherits(x, "htest")) result <- c(x$statistic, x$p.value)
    names(result) <- c("stat", "p-value")
    round(result, 3)
}

