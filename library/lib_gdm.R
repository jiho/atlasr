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

  # setting up a lookup table for names of variables


  lookup.names.variables<-c("day number","tow path","tow configuration","number of turns in the tow","latitude (deg)","longitude (deg)","net depth (m)","depth of the bottom (m)","sea surface temperature (degC)","anomaly in sea surface temperature (degC)","sex unknown","day number of the year", "distance to the nearest rookery (km)","day / night","tow category","tow on the bottom","line of the tow","turn in the tow","fishing year","day number of the fishing year","distance from the nearest rookery (km)","type of processing plant onboard","distance from nearest shore","gear","fishing day number","fishing day number / days in the year","target group","12hr fishing pressure (tows)","fishing year","fishing day","nearest hill","duration of tow","longitude at tow start","month","name of vessel","10 day fishing pressure (tows)","2 day fishing pressure (tows)","depth bin","depth bin","catch weight","longitude","hill number","latitude","target group","hook code","hook size","bait species","soak time","cumulative catch in 5deg * 5deg squares","cumulative catch in 1 deg * 1 deg squares","toothfish CPUE (kg/1000 baited hooks)","temperature","wind direction","wind speed","depth of gear at start","maximum depth of gear","start time","latitude","longitude","bottom depth","minimum water temperature","maximum water temperature","salinity","water temperature","depth","slope","ice15 mean","chla mean","sst sd","SiO mean", "ice15 sd","chla sd","sst mean","PhO mean","depth","slope","mean seaice","mean summer chl","mean summer sst","gradient summer sst","distance from polar front","mean chl of the month","mean sst of the month","gradient sst of the month","normalised time of day","seaice of the month","days since seaice melt","depth","summer chlorophyll","bottom speed","rugosity","bottom temperature","bottom salinity","depth",  "slope","summer chlorophyll", "distance the the ice edge", "distance to the shelf", "floor temperature", "mixed layer depth in summer", "oxygen 500m in summer", "salinity 500m in summer", "si 500m in summer",  "seaice gt85" , "distance to antarctic coast"    )
names(lookup.names.variables)<-c("day.no","towpath","towconfig","no.turns","lat","long","netdepth","botdepth","sst","sst.anom","unk.sex","day.no","dist.col","DN","towcat","towbot","towline","towturn","fyear","dayn","d.col","proctype","d.shore","gear.checked","fdayn","fract.fdayn","targ.g","pres.12h","fyr","fday","name.hill","tow.dur","longs","mth","vessel.name","pres.10d","pres.2d","depth.bin","depth.bin2","catch_weight","long_s","RegNO.hill","lat_s","target.gp","hook.code","hook.size","bait.species","soak.time","vlscc","lscc","lCPUE.Toothfish","temp","wind.dir","wind.sp","gear_s","max_gdepth","time_s","lat_s_dec","long_s_dec","bot_gs","ctd_mint","ctd_maxt", "TSG_Sal", "W_Temp_Hi", "bathymetry_gebco_depth", "bathymetry_gebco_slope", "sea_ice_fraction_time_lt_15_percent_ice_long_term_mean",  "chla_seawifs_annual_maximum_long_term_mean", "sst_oiv2_summer_long_term_sd", "nutrient_Si0_mean", "sea_ice_fraction_time_lt_15_percent_ice_long_term_sd","chla_seawifs_annual_maximum_long_term_sd","sst_oiv2_summer_long_term_mean" ,"nutrient_Ph0_mean","bathymetry_depth_smith_sandwell","bathymetry_slope_smith_sandwell","seaice_gt85_long_term_mean","chl_modis_summer_climatology_mean","sst_modis_summer_climatology_mean","sst_modis_summer_climatology_gradient","distance_polar_front","chl_seawifs_mean_bymonth","sst_oidaily_mean_bymonth","sst_oidaily_gradient_bymonth" ,"normalised_time_of_day","seaice_smmrssmi_mean_bymonth","seaice_smmrssmi_days_since_melt_bymonth","gebco8_modified" ,"chla_summer",     "speed_bottom",    "rugosity2" ,"T_HIGEM_bottom" , "S_HIGEM_bottom","gebco8_mid", "bathymetry_slope", "chl_summer_climatology", "distance_max_ice_edge" ,"distance_shelf", "floor_temperature", "mixed_layer_depth_summer_climatology", "oxygen_500_summer_climatology", "salinity_500_summer_climatology", "si_500_summer_climatology", "seaice_gt85" ,"distance_antarctica"    )



gdm <- function(file, variables, lat.min=-80, lat.max=-30, lat.step=0.1, lon.min=-180, lon.max=180, lon.step=0.5, path="env_data", ...) {

  env <- read.env.data(variables=variables, path=path)

  CPR <- read.data(file)

  # only 10 most abundant
  # resp.vars <- names(sort(colSums(CPR[,!names(CPR) %in% c("lat", "lon")]), decreasing=T)[1:10])
  resp.vars <- names(CPR[,!names(CPR) %in% c("lat", "lon")])

  dat <- associate.env.data(CPR,env)

  pred.data <- build.grid(lat.min=lat.min, lat.max=lat.max, lat.step=lat.step, lon.min=lon.min, lon.max=lon.max, lon.step=lon.step)

  pred.data <- associate.env.data(pred.data, env)

  predvar <- names(env)

  CPR.gdm.all <- compute.gdm(dat=dat, resp.vars=resp.vars, predvar=predvar, samp=1000, pred.data=pred.data, ...)

  return(invisible(CPR.gdm.all))
}



# example of use
# result <- do.gdm(dat=fam,resp.vars=resp.vars,predvar=predvar,samp=15000,pred.data=env.Ant)

compute.gdm <- function(dat, resp.vars, predvar, samp = 10000, plot.name = "GDM", pred.data, do.indicator.species=F, n.clust=NA,...) {

    suppressPackageStartupMessages(require("cluster", quietly=TRUE))

    ## Base case GDM running
    ## dat has your dataset of responses
    ## resp.vars are the names of different variables you want to predict
    ## predvar are your predictive variable names
    ## samp is your sample size, if too big the function will fall over. The function will randomly pick that number of samples from the dataset
    ## plot.name = what name will the plot be saved under, can include a relative path or full path, ignore extension, it will be .eps
    ## pred.data: predictive dataframe, needs lat and long, and environmental variables as in predvar
    ## n.clust = NA if automatically chose cluster number, or give a number
    ## do.indicator.species: calculate indicator species for clusters?
    ## image.name superseded, not saved anymore, left in for code running elsewhere


  # create a new PDF file for the plots
  if (!is.null(plot.name)) {
    pdf(paste(plot.name, ".pdf" ,sep=""), width=10, height=8, useDingbats=FALSE, onefile=TRUE)
  }


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
    }




    #########################################################
    ## now cluster the data and write the results

    pred<-dat[1,pred.var.col]
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

    ## plot the results
    p = polar.ggplot(pred.data, aes(colour=cluster), geom="point")
    # add a title
    p = p + opts(title="GDM")
    # display the plot
    print(p)
    # NB: suppress warnings about missing values: they are necessary to split the coastline in several bits

    ## calculate indicator species using dufrene-legendre method
    if (do.indicator.species) {
        require(labdsv)
        library(reshape)
        temp <- as.matrix(cast(pred.data, lat ~  long, value = "cluster",add.missing=T))
        datlonidx <- round(approx(as.numeric(colnames(temp)),1:ncol(temp),dat$long)$y)
        datlatidx <- round(approx(as.numeric(rownames(temp)),1:nrow(temp),dat$lat)$y)

        dat$cluster <- NA
        for (i in 1:length(datlonidx)){
          if (!is.na(datlonidx) && !is.na(datlatidx)) {
            dat[i,"cluster"] <- temp[datlatidx[i],datlonidx[i]]
          }
        }

        nnanidx=which(!is.na(dat$cluster)) ## only include non-NA clusters
        ## calculate indicator species stuff
        frodo.baggins=indval(dat[nnanidx,resp.var.col],clustering=dat$cluster[nnanidx])
    }

    pred.data$cluster<-factor(pred.data$cluster)
    # dat$cluster <- factor(dat$cluster)


  # close the PDF file
  if (!is.null(plot.name)) {
    dev.off()
  }

    if (do.indicator.species) {
      res=list('predicted'=pred.data,'model'=no.gdm,'plot.obj'=p,'dat.cluster'=dat$cluster,'indval'=frodo.baggins)
    } else {
      res=list('predicted'=pred.data,'model'=no.gdm,'plot.obj'=p)
      }

      return(invisible(res))
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





