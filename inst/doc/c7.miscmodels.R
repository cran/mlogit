## ----label = setup, include = FALSE----------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, widtht = 65)
options(width = 65)

## --------------------------------------------------------------
library("mlogit")
data("ModeCanada", package = "mlogit")
busUsers <- with(ModeCanada, case[choice == 1 & alt == 'bus'])
Bhat <- subset(ModeCanada, ! case %in% busUsers & alt != 'bus' & noalt == 4)
Bhat$alt <- Bhat$alt[drop = TRUE]
Bhat <- dfidx(Bhat, idx = c("case", "alt"), choice = "choice", idnames = c("chid", "alt"))
pcl <- mlogit(choice ~ freq + cost + ivt + ovt, Bhat, reflevel = 'car',
              nests = 'pcl', constPar=c('iv:train.air'))
summary(pcl)

## --------------------------------------------------------------
data("Game", package = "mlogit")
data("Game2", package = "mlogit")
head(Game,2)
head(Game2, 7)
nrow(Game)
nrow(Game2)

## --------------------------------------------------------------
G <- dfidx(Game, varying = 1:12, choice = "ch", ranked = TRUE, idnames = c("chid", "alt"))
G <- dfidx(Game2, choice = "ch", ranked = TRUE, idx = c("chid", "platform"),
           idnames = c("chid", "alt"))
head(G)
nrow(G)

## --------------------------------------------------------------
summary(mlogit(ch ~ own | hours + age, G, reflevel = "PC"))

