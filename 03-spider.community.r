# add environmental data to the dataset
# need the original environmental dataset
# need a dataset called dataset with the lat and long variables (caleld lat and long)
############################################################################


windows(height=11,width=7,rescale="fixed")
setwd("C:/Travail/Projets/ANTARCTIQUE/Brest2011/Sophie/R")
load("env.data.RData")
source("base scripts//load env data.r")



############################################################################
# get environmental data from the dataset
############################################################################

dataset <- read.table("C:/Travail/Projets/ANTARCTIQUE/Brest2011/Sophie/Data/data_butch_update.txt", dec=".", sep="\t", header=T)
names(dataset)
summary(dataset)
names(dataset)[1:2]<-c("lat","long")

res1 <- load.env (dataset=dataset,env.dat=NA,path="C:\\Sophie\\Antarctic\\")
dataset<-res1[["dataset"]]


############################################################################
# checks - some on land
############################################################################

summary(dataset)
length(dataset$bathymetry[dataset$bathymetry>0]) #34 points on land according to bathymetry...
nrow(dataset)
dataset <- dataset[dataset$bathymetry<0,]
nrow(dataset)
summary(dataset)
rm(res1)
save.image("spider.comm.RData")




############################################################################
# do gdm
############################################################################

library(MASS)
require(cluster)
source("base scripts//gdmfuncs.1.1.R")
source("base scripts//lookup.names.variables.R")
source("base scripts//gdm.function.r")

predvar<-names(env.data)[3:ncol(env.data)]
resp.vars <-names(dataset)[! names(dataset) %in% c(names(env.data),"lat","long")]

# chose one in 5 for now
index <- c(1:nrow(env.data))
index <- index %% 5
index <- index == 0
env.sm <- env.data[index,]
env.sm <-env.sm[!is.na(env.sm$bathymetry),]

result <- do.gdm(dat=dataset,resp.vars=resp.vars,predvar=predvar,samp=5000,pred.data=env.sm, plotname = "spider.gdm.effects",image.name="spider.gdm.file")



############################################################################
# do brt on most important ones
############################################################################

library(gbm)
source("base scripts//gbm.functions.2010SM.r")
source("base scripts//brt.function.r")

result <- do.brt(dat=dataset,resp.vars=c("Austropallene.brachyura","Austropallene.cornigera"),predvar=predvar,distrib ="bernoulli", wghts=NULL,monotone=NULL,n.boot = 1, plotname = "spider.brt.effects",image.name="spider.brt.file",n.pred=1,pred.data=env.sm)


#save.image("spider.comm.RData")





