# add environmental data to the dataset
# this one uses the new functionwarnings
############################################################################


# dev.new(height=11, width=7, rescale="fixed")
# setwd("C:\\Projects\\Specific projects\\2011\\Bioregionalisation\\R")
source("base scripts/load env data.r")



############################################################################
# load the dataset and make the
############################################################################

dataset <- read.csv("../Data/Psychroteuthis glacialis.csv")
names(dataset) <- c("long","lat","presence")


#res1 <- load.env (dataset=dataset,env.dat=NA,path="C:\\Projects\\Specific projects\\ANT and Arctic env layers\\Antarctic\\")
res2 <- load.env(dataset=NA, env.dat=T, path="../Antarctic/")

summary(res2[["env.dat"]])
env.data <- res2[["env.dat"]]
env.data$bathymetry[env.data$bathymetry > 0] <- NA

save(env.data, file="env.data.RData")
