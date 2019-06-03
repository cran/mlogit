## ----label = setup, include = FALSE----------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, widtht = 65)
options(width = 65)

## ----nest1-----------------------------------------------------
library("mlogit")
data("HC", package = "mlogit")
HC <- mlogit.data(HC, varying = c(2:8, 10:16), choice = "depvar", shape = "wide")
cooling.modes <- index(HC)$alt %in% c('gcc', 'ecc', 'erc', 'hpc')
room.modes <- index(HC)$alt %in% c('erc', 'er')
# installation / operating costs for cooling are constants, 
# only relevant for mixed systems
HC$icca[!cooling.modes] <- 0
HC$occa[!cooling.modes] <- 0
# create income variables for two sets cooling and rooms
HC$inc.cooling <- HC$inc.room <- 0
HC$inc.cooling[cooling.modes] <- HC$income[cooling.modes]
HC$inc.room[room.modes] <- HC$income[room.modes]
# create an intercet for cooling modes
HC$int.cooling <- as.numeric(cooling.modes)
# estimate the model with only one nest elasticity
nl <- mlogit(depvar ~ ich + och +icca + occa + inc.room + inc.cooling + int.cooling | 0, HC,
             nests = list(cooling = c('gcc','ecc','erc','hpc'), 
             other = c('gc', 'ec', 'er')), un.nest.el = TRUE)
summary(nl)

## --------------------------------------------------------------
 (coef(nl)['iv'] - 1) / sqrt(vcov(nl)['iv', 'iv'])

## --------------------------------------------------------------
# First estimate the multinomial logit model
ml <- update(nl, nests = NULL)
lrtest(nl, ml)

## --------------------------------------------------------------
nl2 <- update(nl, nests = list(central = c('ec', 'ecc', 'gc', 'gcc', 'hpc'), 
                    room = c('er', 'erc')))
summary(nl2)

## ----tstat-----------------------------------------------------
 (coef(nl2)['iv'] - 1) / sqrt(vcov(nl2)['iv', 'iv'])
lrtest(nl2, ml)

## --------------------------------------------------------------
logLik(nl)
logLik(nl2)

## ----nl3-------------------------------------------------------
nl3 <- update(nl, un.nest.el = FALSE)

## ----lrtest1---------------------------------------------------
lrtest(nl, nl3)

## ----threentest------------------------------------------------
nl4 <- update(nl, nests=list(n1 = c('gcc', 'ecc', 'erc'), n2 = c('hpc'),
                    n3 = c('gc', 'ec', 'er')))
summary(nl4)

