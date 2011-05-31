###########################################################################
# Wrap around GDM function
# does model, optimisation of clusters and predictions in one function
# Sophie Mormede May 2011
###########################################################################

## modified May 2011 - BR
## returns a list including the prediction results, the gdm object, the interpolated cluster
##  labels of the training data, and the indicator species information
## note that the training data cluster labels and the indicator species info uses simple spatial
## interpolation, and will not work with time-varying models


# library(MASS)
require(cluster)

# source("gdmfuncs.1.1.R")


# example of use
# result <- do.gdm(dat=fam,resp.vars=resp.vars,predvar=predvar,samp=15000,pred.data=env.Ant)

do.gdm <- function (dat, resp.vars, predvar, samp = 10000, plotname = "gdm.effects",image.name="gdm.file",pred.data, plot.type="wmf", do.indicator.species=F,n.clust=NA) {

    ## Base case GDM running
    ## dat has your dataset of responses
    ## resp.vars are the names of different variables you want to predict
    ## predvar are your predictive variable names
    ## samp is your sample size, if too big the function will fall over. The function will randomly pick that number of samples from the dataset
    ## lookup names variables is where I put general names to your variables names, for plotting.
    ## plotname = what name will the plot be saved under, can include a relative path or full path, ignore extension, it will be .eps
    ## image.name = what name the R workspace will be saved under, can include a relative path or full path, ignore extension, it will be .RData
    ## pred.data: predictive dataframe, needs lat and long, and environmental variables as in predvar
    ## plot.type: can be "wmf", "png" or any of the other types supported by savePlot (default="wmf")
    ## n.clust = NA if automatically chose cluster number, or give a number
    ## do.indicator.species: calculate indicator species for clusters?


    pred.var.col <- which(names(dat) %in% predvar)
    resp.var.col <- which(names(dat) %in% resp.vars)

    datJ <- dat
    for (i in c(pred.var.col,resp.var.col)) {
        datJ<-datJ[!is.na(datJ[,i]),]
    }


    ## #################################
    ## run gdm
    no.gdm <- gdm.fit(datJ[,pred.var.col],datJ[,resp.var.col],sample=samp)


    ## and use gdm.transform to create the curves

    ## make environmental ranges to plot response curves

    env.ranges<-dat[1:200,]
    for (i in pred.var.col) {
        env.ranges[,i]<-seq(from=(Min(dat[,i])),to=Max((dat[,i])),length=200)
    }
    env.ranges<-env.ranges[,pred.var.col]
    no.curves <- gdm.transform(no.gdm,env.ranges)

    ## sort them by max to min
    tp<-names(rev(sort(apply(no.curves,2,max))))
    temp<-match(tp,names(no.curves))


    ## then plot them out
    windows(height=7.4,width=11.6,rescale="fixed")
    par(mfrow=c(3,4))
    par(cex=0.9)
    j<-1
    for (i in temp) {
        if(names(no.curves)[i] %in% names(lookup.names.variables)) {
            tp<-as.character(lookup.names.variables[names(no.curves)[i]])
        } else {
            tp<-names(no.curves)[i]
        }
        plot(env.ranges[,i],no.curves[,i],type='l',cex=1.5,
             xlab=tp, ylab = "transform")
        rug(quantile(dat[,names(no.curves)[i]], probs = seq(0, 1, 0.1), na.rm = TRUE))
        if (j %%12 ==0) {
            savePlot(paste(plotname,j%/%12,".",plot.type,sep=""))
            windows(height=7.4,width=11.6,rescale="fixed")
            par(mfrow=c(3,4))
            par(cex=0.9)
        }
        j<- j+1
    }
    savePlot(paste(plotname,j%/%12+1,".",plot.type,sep=""))




    ## #######################################
    ## now cluster the data and write the results

    pred<-dat[1,pred.var.col]

    ## save some things before we start deleting rows - BR
    pred.data.originaldim=dim(pred.data)
    loncol=which(names(pred.data) %in% c('long','lon','longitude'))[1] ## try to allow for variants of longitude name
    latcol=which(names(pred.data) %in% c('lat','latitude'))[1] ## try to allow for variants of latitude name
    nlon=length(unique(pred.data[,loncol]))
    nlat=length(unique(pred.data[,latcol]))
    ## reconstruct lon and lat grids
    temp=rep(NA,pred.data.originaldim[1])
    temp[as.numeric(rownames(pred.data))]=pred.data[,loncol]
    longrid=matrix(temp,nrow=nlon)
    temp=rep(NA,pred.data.originaldim[1])
    temp[as.numeric(rownames(pred.data))]=pred.data[,latcol]
    latgrid=matrix(temp,nrow=nlon)

    for (i in c(predvar)) {
        pred.data<-pred.data[!is.na(pred.data[,i]),]
    }

    for (i in seq(from=1,to=nrow(pred.data)%/%1000*1000+1, by=1000)) {
        pred<- rbind(pred,gdm.transform(no.gdm,pred.data[i:min(i+999,nrow(pred.data)),colnames(pred.data)%in%colnames(dat)[pred.var.col]]))
    }
    pred<-pred[-1,]


    ## optimise the cluster
   if (is.na(n.clust)) {
     res<-array(NA,c(10,8,4))

      for (cl in 4:10) {
          for (samples in 3:8) {
              for (sampsz in 2:4) {
                  temp <- clara(pred,cl,metric="manhattan",keep.data=FALSE,samples=samples,sampsize=min(nrow(pred), 40 +sampsz*cl))
                  res[cl,samples,sampsz]<-temp$silinfo$avg.width
              }
          }
      }

      temp<-which(res==max(res,na.rm=T),arr.ind=T)
      if (nrow(temp)>1) temp<-temp[1,]
      kclust <- clara(pred,temp[1],metric="manhattan",keep.data=FALSE,samples=temp[2],sampsize=min(nrow(pred), 40 +temp[3]*temp[1]))
    } else {
      kclust <- clara(pred,n.clust,metric="manhattan",keep.data=FALSE)
    }

    pred.data[,"cluster"]<-kclust$clustering

    ## calculate indicator species using dufrene-legendre method
    if (do.indicator.species) {
        require(labdsv)
        ## use predicted clusters on gridded data, and then
        ## interpolate station cluster labels from those gridded cluster labels
        ## assume that grid is regular and rectangular
        ## this code is probably quite fragile! beware! - ben
        ## note that pred.data has been modified to delete rows with NAs, so need to reconstruct full matrix including NAs (in order to get correct shape of grid)
        temp=rep(NA,pred.data.originaldim[1])
        temp[as.numeric(rownames(pred.data))]=pred.data$cluster
        cgrid=matrix(temp,nrow=nlon)
        ## which grid cells are closest to each of our training data points?
        loncol=which(names(dat) %in% c('long','lon','longitude'))[1] ## try to allow for variants of longitude name
        latcol=which(names(dat) %in% c('lat','latitude'))[1] ## try to allow for variants of latitude name
        datlonidx=round(approx(longrid[,1],1:dim(longrid)[1],dat[,loncol])$y)
        ## datlonidx holds the row indices into the lon/lat grids of the training data
        datlatidx=round(approx(latgrid[1,],1:dim(longrid)[2],dat[,latcol])$y)
        ## datlonidx holds the col indices into the lon/lat grids of the training data
        dat.cluster=rep(NA,dim(dat)[1])
        for (k in 1:dim(dat)[1]) {
            dat$cluster[k]=cgrid[datlonidx[k],datlatidx[k]]
        }
        nnanidx=which(!is.na(dat$cluster)) ## only include non-NA clusters
        ## calculate indicator species stuff
        frodo.baggins=indval(dat[nnanidx,resp.var.col],clustering=dat$cluster[nnanidx])
    }

    # save.image(paste(image.name,".RData",sep=""))

    if (do.indicator.species) {
      return(list('predicted'=pred.data,'model'=no.gdm,'dat.cluster'=dat$cluster,'indval'=frodo.baggins))
    } else {
      return(list('predicted'=pred.data,'model'=no.gdm))
      }

}





############################################################################
## some general functions
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





