###########################################################################
# Wrap around BRT function
# does model, effects, bootstraps, predictions in one function
# Sophie Mormede May 2011 s.mormede@niwa.co.nz
###########################################################################


#setwd("C:\\Projects\\Specific projects\\2011\\Bioregionalisation\\R\\base scripts")
#load("..\\fish.data.RData")

# library(bgm)
# source("gbm.functions.2010SM.r")
# example of use
# result <- do.brt(dat=fam,resp.vars=resp.vars[3],predvar=predvar,n.pred=1,pred.data=env.Ant)


do.brt <- function (dat, resp.vars, predvar, int = 2, distrib ="bernoulli", wghts=NULL,monotone=NULL,n.boot = NA, plotname = "brt.effects",image.name="brt.file",n.pred=NA,pred.data ) {

## Base case BRT running - it will output a list of objects with for each response variable: the brt object, the responses, the interactions, predictions etc.

# dat has your dataset of responses
# resp.vars are the names of different variables you want to predict
# predvar are your predictive variable names
# int is the number of interactions, or tree complexity
# distrib is the distribution family: bernoulli (=binomial), poisson, laplace or gaussian
# if adding weights, give the name of the variable in dat dataset which corresponds to weights
# if adding monotone, a vector of 0 and 1 of the lengths of predvar, 1 is monotone increasing and -1 monotone decreasing, 0 is asbitrary
# n.boot is the number of bootstraps for the predictions, NA means single point estimate. Beware, bootstraps are slow to run, and if less than 50-100 bootstraps, will fail.
# plotname = what name will the plot be saved under, can include a relative path or full path, ignore extension, it will be .eps
# image.name = what name the R workspace will be saved under, can include a relative path or full path, ignore extension, it will be .RData
# n.pred is number of bootstraps for predictions, 1 if no bootstraps wanted, and NA if no predictions wanted
# pred.data: predictive dataframe, needs lat and long, and environmental variables as in predvar



result <- list()

  for (resp in resp.vars) {

  result[[resp]] <- list()


  Median<-function(x) {median(x,na.rm=T)}
  Mean<-function(x) {mean(x,na.rm=T)}
  cv<-function(x) {return(sqrt(var(x,na.rm=T))/mean(x,na.rm=T))}


###########################################################################
# setup the dataset
###########################################################################

 dat<-dat[!is.na(dat[,resp]),]
 if (distrib=="bernoulli") {dat[,resp]<-dat[,resp]>0}

 if (!is.null(wghts)) {wgths <- dat$wghts} else {wghts <- rep(1,nrow(dat))}
 if (is.null(monotone)) {monotone <- rep(0,length(predvar))}



###########################################################################
# run the brt without weights
###########################################################################

  cat(paste(resp,": optimising BRT model \n",sep=""))


    lr <- 0.05
    no.trees <- 0

    while ( no.trees < 1000 & lr > 0.0005 ) {

      try( obj <-  gbm.step( data = dat,
                          gbm.x = match(predvar, names(dat)),
                          gbm.y = match(resp, names(dat)),
                          learning.rate = lr,
                          tree.complexity = int,
                          site.weights=wghts,
                          family = distrib,
                          max.trees = 10000,
                          var.monotone=monotone,
                          n.trees = 50,
                          silent=T
                          )
      )
      if ( !is.null (obj) ) {
         if ( object.size( obj ) > 0 ) {
            no.trees <- obj$gbm.call$best.trees
      } }
      else {
        no.trees <- 0
      }
      lr <- lr / 2
    }

  result[[resp]]$obj <- obj


###########################################################################
# write results
###########################################################################
cat("writing results\n")

  temp<-c(obj$gbm.call$family,  obj$gbm.call$tree.complexity,round(as.numeric(as.character(obj$cv.statistics$deviance.mean))/as.numeric(as.character(obj$self.statistics$mean.null)),2),round(min(obj$cv.roc.matrix),2))
  names(temp)<-c("family","number of trees","perc deviance explained","AUC")
  result[[resp]]$deviance <- temp

  # write contributions
  temp<-t(obj$contributions)
  result[[resp]]$contributions <- temp
  rm(temp)




###########################################################################
# plot effects, bootstrapped or not
###########################################################################
cat("plot effects\n")
graphics.off() # this just closes old windows

  if (!is.na(n.boot)) {
    if (n.boot > 1) {
       boot <- NULL
       cat(paste(resp,": running bootstrap on BRT model \n",sep=""))
       try(boot <- gbm.bootstrap(obj,n.reps=n.boot,verbose=F))
       if (!is.null(boot)) {
        result[[resp]]$boot <- boot
        gbm.plot.boot(obj,boot)
      } else {
       gbm.plot(obj)
       cat("Need more bootstraps to provide a confidence interval\n")
       }
     } else {
       gbm.plot(obj)
     }
   }


#save all the plots

#for (i in 1:length(dev.list())) {
#  savePlot(paste(plotname,i,".wmf",sep=""))
#  dev.off(which = dev.cur())
#}




###########################################################################
# run the predictions
###########################################################################

cat("make predictions\n")

  # make predictions

  if (!is.na(n.pred)) {
    if (n.pred>1) {
      temp<-gbm.bootstrap(obj, pred.data = pred.data, return.pred.matrix = T,n.reps=n.pred,verbose=F)
      temp<- temp$pred.matrix
      out <- data.frame(long = pred.data$long, lat = pred.data$lat, Mpred = apply(temp,1,mean), CVpred = apply(temp,1,cv))
      } else {
      out <- data.frame(cbind( pred.data[,c("long","lat")], pred = predict( obj, newdata = pred.data, n.trees = obj$n.trees, type="response")))
      }
   result[[resp]]$pred <- out
  }
save.image(paste(image.name,".RData",sep=""))

}  # end of resp.vars loop

return(result)

} #end of function





############################################################################
# some general functions
############################################################################


Sum <- function(x) {sum(x,na.rm=T)}
Max <-function(x) {max(x,na.rm=T)}
Min<-function(x) {min(x,na.rm=T)}
Median<-function(x) {median(x,na.rm=T)}
Mean<-function(x) {mean(x,na.rm=T)}
Length<-function(x)  length(x[!is.na(x)])


Table <- function(...) {
  a <-table(...,useNA="ifany")
  b <- length(dim(a))
  if (b==2) {
    a<-cbind(a,Sum = margin.table(a,1))
  } else {
    a<-c(a,Sum = margin.table(a))
  }
  return(a)
  }

