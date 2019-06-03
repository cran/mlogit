## ----label = setup, include = FALSE----------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, widtht = 65)
options(width = 65)

## --------------------------------------------------------------
library("mlogit")
data("Fishing", package = "mlogit")
Fish <- mlogit.data(Fishing, shape="wide", varying=2:9, choice="mode")

## ----echo = FALSE, results = 'hide'----------------------------
strt <- c(0.725136619390524, 0.623929046825871, -0.012154264563254, 
          2.40052090993936e-06, -6.54190678882061e-05, 1.54787449680641, 
          0.400099854168002, 1.27472767540054, 1, 0.545695205792339, 0.695440476892784)
Fish.mprobit <- mlogit(mode~price | income | catch, Fish, probit = TRUE,
                       alt.subset=c('beach', 'boat','pier'), start = strt)

## ----eval = FALSE----------------------------------------------
#  Fish.mprobit <- mlogit(mode~price | income | catch, Fish, probit = TRUE, alt.subset=c('beach', 'boat','pier'))

## --------------------------------------------------------------
summary(Fish.mprobit)

