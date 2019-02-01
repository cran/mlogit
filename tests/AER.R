# BankWages.Rd

data("BankWages", package = "AER")

## exploratory analysis of job ~ education
## (tables and spine plots, some education levels merged)
xtabs(~ education + job, data = BankWages)
edcat <- factor(BankWages$education)
levels(edcat)[3:10] <- rep(c("14-15", "16-18", "19-21"), c(2, 3, 3))
tab <- xtabs(~ edcat + job, data = BankWages)
prop.table(tab, 1)
spineplot(tab, off = 0)
plot(job ~ edcat, data = BankWages, off = 0)

## fit multinomial model for male employees
library("nnet")
fm_mnl <- multinom(job ~ education + minority, data = BankWages,
                   subset = gender == "male", trace = FALSE)
summary(fm_mnl)
confint(fm_mnl)

## same with mlogit package
library("mlogit")
fm_mlogit <- mlogit(job ~ 1 | education + minority, data = BankWages, subset = gender == "male", shape = "wide", choice = "job", reflevel = "custodial")
summary(fm_mlogit)

# TravelMode.Rd

data("TravelMode", package = "AER")

## overall proportions for chosen mode
with(TravelMode, prop.table(table(mode[choice == "yes"])))

## travel vs. waiting time for different travel modes
library("lattice")
xyplot(travel ~ wait | mode, data = TravelMode)

## Greene (2003), Table 21.11, conditional logit model
if(require("mlogit")) {
TravelMode$incair <- with(TravelMode, income * (mode == "air"))
tm_cl <- mlogit(choice ~ gcost + wait + incair, data = TravelMode,
  shape = "long", alt.var = "mode", reflevel = "car")
summary(tm_cl)
}


# GSOEP9402.Rd

## data
data("GSOEP9402", package = "AER")

## some convenience data transformations
gsoep <- GSOEP9402
gsoep$year2 <- factor(gsoep$year)

## visualization
plot(school ~ meducation, data = gsoep, breaks = c(7, 9, 10.5, 11.5, 12.5, 15, 18))


## Chapter 5, Table 5.1
library("nnet")
gsoep_mnl <- multinom(
  school ~ meducation + memployment + log(income) + log(size) + parity + year2,
  data = gsoep)
#coeftest(gsoep_mnl)[c(1:6, 1:6 + 14),]
 
## alternatively
if(require("mlogit")) {
gsoep_mnl2 <- mlogit(
  school ~ 0 | meducation + memployment + log(income) + log(size) + parity + year2,
  data = gsoep, shape = "wide", reflevel = "Hauptschule")
#coeftest(gsoep_mnl2)[1:12,]
}

# WinkelmannBoes2009.Rd

## data
data("GSOEP9402", package = "AER")

## some convenience data transformations
gsoep <- GSOEP9402
gsoep$meducation2 <- cut(gsoep$meducation, breaks = c(6, 10.25, 12.25, 18),
  labels = c("7-10", "10.5-12", "12.5-18"))
gsoep$year2 <- factor(gsoep$year)

## Chapter 1
## Table 1.4 plus visualizations
gsoep_tab <- xtabs(~ meducation2 + school, data = gsoep)
round(prop.table(gsoep_tab, 1) * 100, digits = 2)
spineplot(gsoep_tab)
plot(school ~ meducation, data = gsoep, breaks = c(7, 10.25, 12.25, 18))
plot(school ~ meducation, data = gsoep, breaks = c(7, 9, 10.5, 11.5, 12.5, 15, 18))


## Chapter 5
## Table 5.1
library("nnet")
gsoep_mnl <- multinom(
  school ~ meducation + memployment + log(income) + log(size) + parity + year2,
  data = gsoep)
#coeftest(gsoep_mnl)[c(1:6, 1:6 + 14),]
 
## alternatively
if(require("mlogit")) {
gsoep_mnl2 <- mlogit(school ~ 0 | meducation + memployment + log(income) +
  log(size) + parity + year2, data = gsoep, shape = "wide", reflevel = "Hauptschule")
coeftest(gsoep_mnl2)[1:12,]
}

##################################
## Choice of Brand for Crackers ##
##################################

## data
if(require("mlogit")) {
data("Cracker", package = "mlogit")
head(Cracker, 3)
crack <- mlogit.data(Cracker, varying = 2:13, shape = "wide", choice = "choice")
head(crack, 12)

## Table 5.6 (model 3 probably not fully converged in W&B)
crack$price <- crack$price/100
crack_mlogit1 <- mlogit(choice ~ price | 0, data = crack, reflevel = "private")
crack_mlogit2 <- mlogit(choice ~ price | 1, data = crack, reflevel = "private")
crack_mlogit3 <- mlogit(choice ~ price + feat + disp | 1, data = crack,
  reflevel = "private")
lrtest(crack_mlogit1, crack_mlogit2, crack_mlogit3)

## IIA test
crack_mlogit_all <- update(crack_mlogit2, reflevel = "nabisco")
crack_mlogit_res <- update(crack_mlogit_all,
  alt.subset = c("keebler", "nabisco", "sunshine"))
hmftest(crack_mlogit_all, crack_mlogit_res)
}
