library("mlogit")

data("TravelMode", package = "AER")
hl <- hlogit(choice~wait+travel+vcost,TravelMode,
             shape="long",id.var="individual",alt.var="mode",
             choice="choice",print.level=3,method="bfgs")

data("Fishing",package="mlogit")
Fish <- mlogit.data(Fishing,varying=c(4:11),shape="wide",
        choice="mode",opposite=c("pr"))
rlf<- rlogit(mode~pr+ca,data=Fish,rpar=c(ca="n"),R=100,
       halton=NA,print.level=3,norm="pr",method="bhhh")

data("Train",package="Ecdat")
Train <- mlogit.data(Train,choice="choice",varying=4:11,sep="",
 alt.levels=c("ch1","ch2"),shape="wide",
 opposite=c("price","change","comfort","time"))
rlt <- rlogit(choice~price+time+change+comfort-1,data=Train,
            rpar=c(change="n",comfort="n",time="n"),R=20,
            halton=NA,print.level=3,id="id",
            correlation=TRUE,norm="price",method="bfgs")


data("TravelMode", package = "AER")
TravelMode$avincome <- with(TravelMode, income * (mode == "air"))
TravelMode$time <- with(TravelMode, travel+wait)/60
TravelMode$timeair <- with(TravelMode,time*I(mode=="air"))
TravelMode$income <- with(TravelMode,income/10)
#Heiss p.231
nl <- nlogit(choice~time+timeair|income,TravelMode,choice="choice",shape="long",alt.var="mode",print.level=3,method="bfgs",nest=list(public=c("train","bus"), other=c("air","car")))

save(file="longExamples.rda",rlt,rlf,hl,nl)
