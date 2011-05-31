################################################################
# plotting routines to be adapted, use polar
# Sophie Mormede, May 2011 s.mormede@niwa.co.nz
################################################################

library(polar)
source("C:/Projects/Library/computing/R functions/polar.plots.R")
source("C:/Projects/Library/computing/R functions/image.scale.R")
source("C:/Projects/Library/computing/R functions/Hist3d.R")


################################################################
# plot for BRT
# item called result
# use this routine after the main function
################################################################

  # Plot the results - here for presence / absence

  res <- result[[1]]$pred

  breaks<-seq(from=0,to=1,length=6)
  temp<-cut(res$pred,breaks,labels=F)
  temp<-factor(temp,levels=c(1:5),labels=c("dark blue", "green",  "yellow","orange","red"))

  dev.new(width=7, height=7,rescale="fixed")
  polar.plot(xlim = c(0, 360), ylim = c(-50, -90),rotate=F)
  polar.points(res$long,res$lat,col=as.character(temp),pch=15,cex=0.5)
  polar.coastline()
  polar.clip()
  polar.axis(cex=0.7,lwd=1.5,line=0.4)
  invisible()
  leg<-paste(breaks,"-",breaks[-1],sep="")[-6]
  polar.legend(9,-75,pch=15,col=rev(levels(temp)),legend=rev(leg),cex=0.8,bty="n")






################################################################
# plot for GDM
# item called result
# use this routine after the main function
################################################################


  library(polar)
  source("C:/Projects/Library/computing/R functions/polar.plots.R")
  source("C:/Projects/Library/computing/R functions/image.scale.R")
  source("C:/Projects/Library/computing/R functions/Hist3d.R")

  Col<-c(9,13,16)
  Median<-function(x) {median(x,na.rm=T)}
  Mean<-function(x) {mean(x,na.rm=T)}

  dev.new(height=7,width=7,rescale="fixed")
  polar.plot(xlim = c(0,360), ylim = c(-45, -90),rotate=F)
  polar.points(result$long,result$lat,col=as.numeric(as.character(result$cluster)),pch=".",cex=5)
  polar.depth(2000,col=gray(0.8))
  polar.depth(1000,col=gray(0.8))
  polar.coastline()
  polar.clip()
  polar.axis(cex=0.7,lwd=1.5)
  invisible()
  savePlot("gdm.pred.png",type="png")

