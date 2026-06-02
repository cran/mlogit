## -----------------------------------------------------------------------------
#| echo: false
oopts <-options(width = 70)


## -----------------------------------------------------------------------------
#| label: pcl
#| message: false
library(mlogit)
busUsers <- with(ModeCanada, case[choice == 1 & alt == 'bus'])
Bhat <- subset(ModeCanada, ! case %in% busUsers &
                           alt != 'bus' & noalt == 4)
Bhat$alt <- Bhat$alt[drop = TRUE]
Bhat <- dfidx(Bhat, idx = c("case", "alt"), choice = "choice",
              idnames = c("chid", "alt"))
pcl <- mlogit(choice ~ freq + cost + ivt + ovt, Bhat, reflevel = 'car',
              nests = 'pcl', constPar=c('iv:train.air'))
gaze(pcl)


## -----------------------------------------------------------------------------
#| label: "game data set"
#| collapse: true
data("Game", package = "mlogit")
data("Game2", package = "mlogit")
head(Game,2)
head(Game2, 7)
nrow(Game)
nrow(Game2)


## -----------------------------------------------------------------------------
#| label: "dfidx game"
G <- dfidx(Game, varying = 1:12, choice = "ch", ranked = TRUE,
           idnames = c("chid", "alt"))
G <- dfidx(Game2, choice = "ch", ranked = TRUE, idx = c("chid", "platform"),
           idnames = c("chid", "alt"))
print(G, n = 3)


## -----------------------------------------------------------------------------
#| label: "rol estimation"
mlogit(ch ~ own | hours + age, G, reflevel = "PC") |> gaze()


## -----------------------------------------------------------------------------
#| echo: false
options(oopts)

