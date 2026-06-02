## -----------------------------------------------------------------------------
#| echo: false
oopts <-options(width = 70)


## -----------------------------------------------------------------------------
#| label: nest1
library(mlogit)
data("HC", package = "mlogit")
HC <- dfidx(HC, varying = c(2:8, 10:16), choice = "depvar")
cooling.modes <- idx(HC, 2) %in% c('gcc', 'ecc', 'erc', 'hpc')
room.modes <- idx(HC, 2) %in% c('erc', 'er')
# installation / operating costs for cooling are constants, 
# only relevant for mixed systems
HC$icca[! cooling.modes] <- 0
HC$occa[! cooling.modes] <- 0
# create income variables for two sets cooling and rooms
HC$inc.cooling <- HC$inc.room <- 0
HC$inc.cooling[cooling.modes] <- HC$income[cooling.modes]
HC$inc.room[room.modes] <- HC$income[room.modes]
# create an intercet for cooling modes
HC$int.cooling <- as.numeric(cooling.modes)
# estimate the model with only one nest elasticity
nl <- mlogit(depvar ~ ich + och +icca + occa + inc.room +
                 inc.cooling + int.cooling | 0, HC,
             nests = list(cooling = c('gcc','ecc','erc','hpc'), 
             other = c('gc', 'ec', 'er')), un.nest.el = TRUE)
gaze(nl)


## -----------------------------------------------------------------------------
#| collapse: true
unname( (coef(nl)['iv'] - 1) / sqrt(vcov(nl)['iv', 'iv']))


## -----------------------------------------------------------------------------
#| collapse: true
ml <- update(nl, nests = NULL)
lrtest(nl, ml) |> gaze()


## -----------------------------------------------------------------------------
nl2 <- update(nl,
              nests = list(central = c('ec', 'ecc', 'gc', 'gcc', 'hpc'), 
              room = c('er', 'erc')))
gaze(nl2)


## -----------------------------------------------------------------------------
#| label: tstat
#| collapse: true
unname((coef(nl2)['iv'] - 1) / sqrt(vcov(nl2)['iv', 'iv']))
lrtest(nl2, ml) |> gaze()


## -----------------------------------------------------------------------------
#| collapse: true
logLik(nl)
logLik(nl2)


## -----------------------------------------------------------------------------
#| label: nl3
nl3 <- update(nl, un.nest.el = FALSE)


## -----------------------------------------------------------------------------
#| collapse: true
#| label: lrtest1
lrtest(nl, nl3) |> gaze()


## -----------------------------------------------------------------------------
#| label: threentest
#| collapse: true
nl4 <- update(nl, nests=list(n1 = c('gcc', 'ecc', 'erc'), n2 = c('hpc'),
                    n3 = c('gc', 'ec', 'er')))
logLik(nl4)


## -----------------------------------------------------------------------------
#| echo: false
options(oopts)

