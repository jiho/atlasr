# add environmental data to the dataset
# need the original environmental dataset
# need a dataset called dataset with the lat and long variables (caleld lat and long)
############################################################################


dev.new(height=11,width=7,rescale="fixed")
# setwd("C:/Travail/Projets/ANTARCTIQUE/Brest2011/Sophie/R")
load("env.data.RData")
source("base scripts/load env data.r")

dataset <- read.csv("../Data/squid community.csv")
for (i in (3:28)){
  dataset[,i]<-ifelse( dataset[,i]>1,1,dataset[,i])
}

res1 <- load.env(dataset=dataset,env.dat=NA,path="../Antarctic/")
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

env.data<- env.data[!is.na(env.data$bathymetry),]

# save.image("squid.comm.RData")




############################################################################
# do gdm - works on a subset...
############################################################################

library(MASS)
require(cluster)
source("base scripts/gdmfuncs.1.1.R")
source("base scripts/lookup.names.variables.R")
source("base scripts/gdm.function.r")

predvar <- names(env.data)[3:52] #use all of them
resp.vars <-names(dataset)[! names(dataset) %in% c(names(env.data),"lat","long")]

# chose one in 5 for now
index <- c(1:nrow(env.data))
index <- index %% 5
index <- index == 0
env.sm <- env.data[index,]


dat <- dataset[sample(nrow(dataset),2000),]

result <- do.gdm(dat=dat,resp.vars=resp.vars,predvar=predvar,samp=500,pred.data=env.sm, plotname = "squid.gdm.effects",image.name="squid.gdm.file")



