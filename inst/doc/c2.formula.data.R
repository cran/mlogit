## -----------------------------------------------------------------------------
#| echo: false
oopts <-options(width = 70)


## -----------------------------------------------------------------------------
#| label: "loading mlogit"
#| message: false
library(mlogit)


## -----------------------------------------------------------------------------
#| label: "Train data"
data("Train", package = "mlogit")
Train$choiceid <- 1:nrow(Train)
head(Train, 3)


## -----------------------------------------------------------------------------
#| label: "dfidx for Train"
Tr <- dfidx(Train, shape = "wide", varying = 4:11, sep = "_",
            choice = "choice",
            idx = c(id = "choiceid"), idnames = c(NA, "alt"),
            opposite = c("price", "comfort", "time", "change"))


## -----------------------------------------------------------------------------
#| label: "data transformation for Train"
Tr$price <- Tr$price / 100 / 2.20371
Tr$time <- Tr$time / 60


## -----------------------------------------------------------------------------
#| label: "head of the transformed Train data set"
print(Tr, n = 3)


## -----------------------------------------------------------------------------
#| label: "index of the transformed Train data set"
head(idx(Tr), 3)


## -----------------------------------------------------------------------------
#| label: "loading ModeCanada"
data("ModeCanada", package = "mlogit")
head(ModeCanada, 3)


## -----------------------------------------------------------------------------
#| label: "applying dfidx to Modecanada (1)"
MC <- dfidx(ModeCanada, subset = noalt == 4,
          alt.levels = c("train", "air", "bus", "car"))


## -----------------------------------------------------------------------------
#| label: "applying dfidx to Modecanada (2)"
MC <- dfidx(ModeCanada, subset = noalt == 4, idx = list(NA, "alt"))


## -----------------------------------------------------------------------------
#| label: "applying dfidx to Modecanada (3)"
MC <- dfidx(ModeCanada, subset = noalt == 4, idx = "case",
            alt.levels = c("train", "air", "bus", "car"))


## -----------------------------------------------------------------------------
#| label: "applying dfidx to Modecanada (4)"
MC <- dfidx(ModeCanada, subset = noalt == 4, idx = c("case", "alt"))


## -----------------------------------------------------------------------------
#| label: "ModeCanada without idx"
MC <- dfidx(ModeCanada, subset = noalt == 4)


## -----------------------------------------------------------------------------
#| label: "applying dfidx to Modecanada (5)"
MC <- dfidx(ModeCanada, subset = noalt == 4, idx = c("case", "alt"),
            drop.index = FALSE)
print(MC, n = 3)


## -----------------------------------------------------------------------------
#| label: "a three parts formula"
library(Formula)
f <- Formula(choice ~ cost | income + urban | ivt)


## -----------------------------------------------------------------------------
#| label: "ommission of some parts (1)"
f2 <- Formula(choice ~ cost + ivt | income + urban)
f2 <- Formula(choice ~ cost + ivt | income + urban | 0)


## -----------------------------------------------------------------------------
#| label: "ommission of some parts (2)"
f3 <- Formula(choice ~ 0 | income | 0)
f3 <- Formula(choice ~ 0 | income)


## -----------------------------------------------------------------------------
#| label: "ommission of some parts (3)"
f4 <- Formula(choice ~ cost + ivt)
f4 <- Formula(choice ~ cost + ivt | 1)
f4 <- Formula(choice ~ cost + ivt | 1 | 0)


## -----------------------------------------------------------------------------
#| label: "removing the intercept"
f5 <- Formula(choice ~ cost | income + 0 | ivt)
f5 <- Formula(choice ~ cost | income - 1 | ivt)


## -----------------------------------------------------------------------------
#| label: "model.matrix method for Formula objects"
MC <- dfidx(ModeCanada, subset = noalt == 4,
            alt.levels = c("train", "air", "bus", "car"),
            pkg = "mlogit")
class(MC)
f <- Formula(choice ~ cost | income  | ivt)
mf <- model.frame(MC, f)
class(mf)


## -----------------------------------------------------------------------------
print(model.matrix(mf), n = 3)


## -----------------------------------------------------------------------------
#| echo: false
options(oopts)

