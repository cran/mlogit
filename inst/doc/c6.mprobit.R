## ----label = setup, include = FALSE----------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, widtht = 65)
options(width = 65)

## --------------------------------------------------------------
library("mlogit")
data("Fishing", package = "mlogit")
Fish <- dfidx(Fishing, varying = 2:9, choice = "mode", idnames = c("chid", "alt"))

## --------------------------------------------------------------
Fish.mprobit <- mlogit(mode~price | income | catch, Fish, probit = TRUE, alt.subset=c('beach', 'boat','pier'))

## --------------------------------------------------------------
summary(Fish.mprobit)

