## -----------------------------------------------------------------------------
#| echo: false
oopts <-options(width = 70)


## -----------------------------------------------------------------------------
#| label: "fishing"
#| message: false
library(mlogit)
Fish <- dfidx(Fishing, varying = 2:9, choice = "mode",
              idnames = c("chid", "alt"))


## -----------------------------------------------------------------------------
#| label: "mprobit estimation"
Fish.mprobit <- mlogit(mode~price | income | catch, Fish, probit = TRUE,
                       alt.subset=c('beach', 'boat','pier'))


## -----------------------------------------------------------------------------
#| label: "mprobit summary"
summary(Fish.mprobit)


## -----------------------------------------------------------------------------
#| echo: false
options(oopts)

