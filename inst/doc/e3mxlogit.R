## ----label = setup, include = FALSE----------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, widtht = 65)
options(width = 65)

## --------------------------------------------------------------
library("mlogit")
data("Electricity", package = "mlogit")
Electr <- mlogit.data(Electricity, id = "id", choice = "choice", 
                      varying = 3:26, shape = "wide", sep = "")

## ----Elec.mxl, echo = FALSE, results = 'hide'------------------
strt <- c(-0.9901049, -0.1985715, 2.0572983, 1.5908743, -9.1176130, -9.1946436, 0.2140892, 
          0.4019263, 1.5845577, 1.0130424, 2.5459903, 0.8854062)
strt <- c(-0.9733844, -0.2055565,  2.0757333,  1.4756497, -9.0525423, -9.1037717,
          0.2199450,  0.3783044,  1.4829803,  1.0000609,  2.2894889,  1.1808827)
Elec.mxl <- mlogit(choice ~ pf + cl + loc + wk + tod + seas | 0, Electr, 
              rpar=c(pf = 'n', cl = 'n', loc = 'n', wk = 'n', 
                     tod = 'n', seas = 'n'), 
              R = 100, halton = NA, panel = TRUE, start = strt)

## ----eval = FALSE----------------------------------------------
#  Elec.mxl <- mlogit(choice ~ pf + cl + loc + wk + tod + seas | 0, Electr,
#                rpar=c(pf = 'n', cl = 'n', loc = 'n', wk = 'n',
#                       tod = 'n', seas = 'n'),
#                R = 100, halton = NA, panel = TRUE)

## --------------------------------------------------------------
summary(Elec.mxl)

## --------------------------------------------------------------
coef(Elec.mxl)['cl'] / coef(Elec.mxl)['pf']

## --------------------------------------------------------------
pnorm(- coef(Elec.mxl)['cl'] / coef(Elec.mxl)['sd.cl'])

## --------------------------------------------------------------
pnorm(- coef(Elec.mxl)['pf'] / coef(Elec.mxl)['sd.pf'])

## ----Elec.mxl2, echo = FALSE, results = 'hide'-----------------
strt <- c(-0.8710529, -0.2082428,  2.0460954,  1.4473149,  8.4091315,  8.5432555,  
          0.3687390,  1.5773527,  0.8837160,  2.5638874,  2.0722178)
strt <- c(-0.8799042, -0.2170603,  2.0922916,  1.4908937, -8.5818566, -8.5832956,
          0.3734776,  1.5588576,  1.0508114,  2.6946672,  1.9507270)
Elec.mxl2 <- mlogit(choice ~ pf + cl + loc + wk + tod + seas | 0, Electr, 
                   rpar = c(cl = 'n', loc = 'n', wk = 'n', 
                            tod = 'n', seas = 'n'), 
                   R = 100, halton = NA,  panel = TRUE, start = strt)

## ----eval = FALSE----------------------------------------------
#  Elec.mxl2 <- mlogit(choice ~ pf + cl + loc + wk + tod + seas | 0, Electr,
#                     rpar = c(cl = 'n', loc = 'n', wk = 'n',
#                              tod = 'n', seas = 'n'),
#                     R = 100, halton = NA,  panel = TRUE)

## --------------------------------------------------------------
summary(Elec.mxl2)

## ----Elec.mxl3, echo = FALSE, results = 'hide'-----------------
strt <- c(-0.8685207, -0.2103447,  2.0269971,  1.4773713,  8.3994921,  8.4976319, 
          0.3693250,  1.5862809,  1.5916990,  2.5775540,  2.0405350)
strt <- c(-0.9303806,  -0.2478098,   2.3808084,   1.5921023,  -5.8173333, -10.7742475,
          0.4115700,   1.4761539,   1.3644855,   5.1445848,   2.7185711)
strt <- c(-0.8822317, -0.2171273,  2.0993191,  1.5094101, -8.6070022, -8.6024084,
          0.3810701,  1.5938502,  1.7863766,  2.7190780,  1.9453765)
Elec.mxl3 <- update(Elec.mxl, rpar = c(cl = 'n', loc = 'n', wk = 'u', 
                                       tod = 'n', seas = 'n'), start = strt)

## ----eval = FALSE----------------------------------------------
#  Elec.mxl3 <- update(Elec.mxl, rpar = c(cl = 'n', loc = 'n', wk = 'u',
#                                         tod = 'n', seas = 'n'))

## --------------------------------------------------------------
summary(Elec.mxl3)
rpar(Elec.mxl3, 'wk')
summary(rpar(Elec.mxl3, 'wk'))

## --------------------------------------------------------------
plot(rpar(Elec.mxl3, 'wk'))

## --------------------------------------------------------------
Electr <- mlogit.data(Electricity, id = "id", choice = "choice", 
                      varying = 3:26, shape = "wide", sep = "",
                      opposite = c('tod', 'seas'))

## ----Elec.mxl4, echo = FALSE, results = 'hide'-----------------
strt <- c(-0.8689874, -0.2113327,  2.0238880,  1.4791236,  2.1123811,  2.1242071,
 0.3731202,  1.5485101,  1.5217919,  0.3670763,  0.2753497)
Elec.mxl4 <- mlogit(choice ~ pf + cl + loc + wk + tod + seas | 0, Electr, 
              rpar = c(cl = 'n', loc = 'n', wk = 'u', tod = 'ln', seas = 'ln'), 
              R = 100, halton = NA, panel = TRUE, start = strt)

## ----eval = FALSE----------------------------------------------
#  Elec.mxl4 <- mlogit(choice ~ pf + cl + loc + wk + tod + seas | 0, Electr,
#                rpar = c(cl = 'n', loc = 'n', wk = 'u', tod = 'ln', seas = 'ln'),
#                R = 100, halton = NA, panel = TRUE)

## --------------------------------------------------------------
summary(Elec.mxl4)

## --------------------------------------------------------------
plot(rpar(Elec.mxl4, 'seas'))

## ----Elec.mxl5, echo = FALSE, results = 'hide'-----------------
strt <-  c(-0.917703974, -0.215851727,  2.392570989,  1.747531863,  2.155462393,
           2.169548103,  0.396252325,  0.617497150, -2.071718067,  0.195238185,
           -1.236664544,  0.643190285,  0.001982314,  0.062508396,  0.160672338,
           0.375855648,  0.025996362, -0.001225349,  0.141381623,  0.089990150,
           0.211244575)
#Elec.mxl5 <- update(Elec.mxl4, correlation = TRUE, start = strt)

## ----eval = FALSE----------------------------------------------
#  Elec.mxl5 <- update(Elec.mxl4, correlation = TRUE)

## --------------------------------------------------------------
#summary(Elec.mxl5)
#cor.mlogit(Elec.mxl5)
#lrtest(Elec.mxl5, Elec.mxl4)
#waldtest(Elec.mxl5, correlation = FALSE)
#scoretest(Elec.mxl4, correlation = TRUE)

