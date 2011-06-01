# add environmental data to the dataset
# need the original environmental dataset
# need a dataset called dataset with the lat and long variables (caleld lat and long)
############################################################################


load.env <- function(dataset=NA,env.dat=NA,path="C:\\Projects\\Specific projects\\ANT and Arctic env layers\\Antarctic\\",lat.lim=c(-80,-30),lat.step=0.1,long.lim=c(-180,180),long.step=0.5) {


  # if dataset = NA, no dataset to match
  # if env.dat = NA, no environmental dataset to extract, anything else will bring it out
  # path is where the environmental dataset is
  #lat.lim = limits of lat for the predictions
  #lat.step = distance between predicitve points in lat
  #long.lim = limits of lat for the predictions
  #long.step = distance between predicitve points in long




############################################################################
# make the environmental prediction dataset
# every 0.5 deg long and 0.1 deg lat (Ben's data is 0.1 by 0.1)
############################################################################
 if (!is.na(env.dat)) {
    lat <- seq(from=lat.lim[1],to=lat.lim[2],by=lat.step)
    long <- seq(from=long.lim[1],to=long.lim[2],by=long.step)

    env.dat <- data.frame(lat=rep(lat,each=length(long)),long=rep(long,times=length(lat)))
    nrow(env.dat) # 217,021 for now

    rm(lat,long)
  }

############################################################################
# get the environmental variables names
############################################################################

env.data <- c("bathymetry","bathymetry_slope","chl_summer_climatology","distance_colony","distance_max_ice_edge","distance_shelf","floor_temperature","mixed_layer_depth_summer_climatology",
"nox_0_summer_climatology","nox_0_winter_climatology","nox_50_summer_climatology","nox_50_winter_climatology","nox_200_summer_climatology","nox_200_winter_climatology","nox_500_summer_climatology","nox_500_winter_climatology",
"oxygen_0_summer_climatology","oxygen_0_winter_climatology","oxygen_50_summer_climatology","oxygen_50_winter_climatology","oxygen_200_summer_climatology","oxygen_200_winter_climatology","oxygen_500_summer_climatology","oxygen_500_winter_climatology",
"salinity_0_summer_climatology","salinity_0_winter_climatology","salinity_50_summer_climatology","salinity_50_winter_climatology","salinity_200_summer_climatology","salinity_200_winter_climatology","salinity_500_summer_climatology","salinity_500_winter_climatology",
"si_0_summer_climatology","si_0_winter_climatology","si_50_summer_climatology","si_50_winter_climatology","si_200_summer_climatology","si_200_winter_climatology","si_500_summer_climatology","si_500_winter_climatology",
"t_0_summer_climatology","t_0_winter_climatology","t_50_summer_climatology","t_50_winter_climatology","t_200_summer_climatology","t_200_winter_climatology","t_500_summer_climatology","t_500_winter_climatology",
"nox_0_interpolated_summer_climatology","nox_0_interpolated_winter_climatology","nox_50_interpolated_summer_climatology","nox_50_interpolated_winter_climatology","nox_200_interpolated_summer_climatology","nox_200_interpolated_winter_climatology","nox_500_interpolated_summer_climatology","nox_500_interpolated_winter_climatology",
"oxygen_0_interpolated_summer_climatology","oxygen_0_interpolated_winter_climatology","oxygen_50_interpolated_summer_climatology","oxygen_50_interpolated_winter_climatology","oxygen_200_interpolated_summer_climatology","oxygen_200_interpolated_winter_climatology","oxygen_500_interpolated_summer_climatology","oxygen_500_interpolated_winter_climatology",
"salinity_0_interpolated_summer_climatology","salinity_0_interpolated_winter_climatology","salinity_50_interpolated_summer_climatology","salinity_50_interpolated_winter_climatology","salinity_200_interpolated_summer_climatology","salinity_200_interpolated_winter_climatology","salinity_500_interpolated_summer_climatology","salinity_500_interpolated_winter_climatology",
"si_0_interpolated_summer_climatology","si_0_interpolated_winter_climatology","si_50_interpolated_summer_climatology","si_50_interpolated_winter_climatology","si_200_interpolated_summer_climatology","si_200_interpolated_winter_climatology","si_500_interpolated_summer_climatology","si_500_interpolated_winter_climatology",
"t_0_interpolated_summer_climatology","t_0_interpolated_winter_climatology","t_50_interpolated_summer_climatology","t_50_interpolated_winter_climatology","t_200_interpolated_summer_climatology","t_200_interpolated_winter_climatology","t_500_interpolated_summer_climatology","t_500_interpolated_winter_climatology",
"seaice_gt85","sst_summer_climatology","distance_antarctica")
env.nam <- c("depth", "slope","mass_concentration_of_chlorophyll_in_sea_water", "distance", "distance", "distance", "temperature", "ocean_mixed_layer_thickness", "nitrate", "nitrate", "nitrate", "nitrate", "nitrate", "nitrate", "nitrate", "nitrate", "oxygen", "oxygen", "oxygen", "oxygen", "oxygen", "oxygen", "oxygen", "oxygen", "salinity", "salinity", "salinity", "salinity", "salinity", "salinity", "salinity", "salinity", "silicate", "silicate", "silicate", "silicate", "silicate", "silicate", "silicate", "silicate", "temperature", "temperature", "temperature", "temperature", "temperature", "temperature", "temperature", "temperature", "nitrate", "nitrate", "nitrate", "nitrate", "nitrate", "nitrate", "nitrate", "nitrate", "oxygen", "oxygen", "oxygen", "oxygen", "oxygen", "oxygen", "oxygen", "oxygen", "salinity", "salinity", "salinity", "salinity", "salinity", "salinity", "salinity", "salinity", "silicate", "silicate", "silicate", "silicate", "silicate", "silicate", "silicate", "silicate", "temperature", "temperature", "temperature", "temperature", "temperature", "temperature", "temperature", "temperature",              "sea_ice", "sea_surface_temperature","distance")




############################################################################
# get the environmental values for the Antarctic
############################################################################


library(ncdf)
library(fields)

for (i in 1:length(env.data)) {
  filename = paste(path,env.data[i],"_0.1_0.1.nc",sep="")
  ncid=open.ncdf(filename)

  # get the lon and lats
  lon=get.var.ncdf(ncid,"lon")
  lat=get.var.ncdf(ncid,"lat")

  # now retrieve the actual data
  dat=get.var.ncdf(ncid,env.nam[i]) #note, variable name has to match the one in the file
  # get the missing value code
  missing_value=att.get.ncdf(ncid,env.nam[i],"_FillValue")
  # replace any such data with NA
  dat[which(dat==missing_value$value)]=NA

  close(ncid) # close the file now that we're finished with it

  # display data as a reality check
  #image(lon,lat,dat,main=env.data[i])

  # obtain interpolated data at desired locations
  dat2=list(x=lon,y=lat,z=dat)
  if (length(dataset)>1 && !is.na(dataset)) {
    zhat=interp.surface(dat2,dataset[,c("long","lat")])
    dataset[,env.data[i]]<-zhat
  }

if (length(env.dat)> 1 && !is.na(env.dat)) {
  zhat=interp.surface(dat2,env.dat[,c("long","lat")])
  env.dat[,env.data[i]]<-zhat
}

}

temp<-list()
temp[["dataset"]] <- dataset
temp[["env.dat"]]<-env.dat
return(temp)

}



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

