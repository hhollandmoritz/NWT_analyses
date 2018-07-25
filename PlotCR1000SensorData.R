# Plots sample niwot data
# here data from new CR1000 loggers
# uploaded by Jen July 2018
# manuall added header to each column for ease of plotting
remove(list=ls())

help(package = CampbellLogger) 
#install.packages("devtools")
library(devtools)
#devtools::install_github("MarkusLoew/CampbellLogger")
library(CampbellLogger)
library(zoo)
library(polynom)

estimate_mode <- function(x) {
  d <- density(x, na.rm=T)
  d$x[which.max(d$y)]
}


dir  <- "/Users/wwieder/Desktop/Working_files/Niwot/Sensors/Sensor_data/sdl_sn_dataJuly18/"
setwd(dir)

for (i in 6:21) {
  if (i <  10) {Node <- paste('0',i ,sep='')}
  if (i >= 10) {Node <- paste(    i ,sep='')}
  fin  <- paste(dir,"sn_0",Node,"_cr1000x_tenminute.dat", sep='')
  dat  <- CampbellFileImport(fin, time.zone = "MST", skip.rows = 4)
# names(dat)  
# datetime            <- as.POSIXct(paste(dat$date, dat$time), format="%Y-%m-%d %H:%M:%S")
  datetime            <- dat$TIMESTAMP
  T5                  <- dat$soiltemp_5cm_Avg
  T30                 <- dat$soiltemp_30cm_Avg
  M5a                 <- dat$vwca_5cm_Avg
  M30a                <- dat$vwca_30cm_Avg
  height              <- dat$dts
  height[height>=4900]<- NA
  height[height<=500] <- NA
  height_ma           <- rollapply(height,24*6, mean, na.rm=TRUE)
  good <- range(height, na.rm=T)[1]
  if (good != "Inf") {      #no data for #10, 13
    plot(datetime, height, main=Node, type='l')
    plot(height, main=Node, type='l')
    
    # manually remove values > max heightt
    # could also use median (or mean) height, but 
    # looks like mode works best

    if (i==6 || i==9 || i==11) {
      maxH <- median(height_ma[1:5000], na.rm=T) + 50
      modH <- estimate_mode(height[1:5000]) 
    } else if (i ==14) {
      maxH <- median(height_ma[1e4:1.5e4], na.rm=T) + 50
      modH <- estimate_mode(height[1e4:1.5e4]) 
    } else if (i ==16) {
      maxH <- median(height_ma[3e4:5.5e4], na.rm=T) + 50
      modH <- estimate_mode(height[3e4:4.5e4]) 
    } else if (i ==21) {
      maxH <- median(height_ma[5e3:1e4], na.rm=T) + 50
      modH <- estimate_mode(height[5e3:1e4]) 
    } else  {
      maxH <- median(height_ma[1:100], na.rm=T) + 50
      modH <- estimate_mode(height[1:1000]) 
    }

    plot(height, type='l')
    lines(height_ma, col=2)
    abline(h=maxH)
    abline(h=modH, col=4)
    
    plot(height, type='l', xlim=c(0,5000))
    lines(height_ma, col=2)
    abline(h=maxH)
    abline(h=modH, col=4)
    
    height_ma <- rollapply(height,24*6, mean, na.rm=TRUE, align="center")
    
    plot(height, type='l')
    lines(height_ma, col=4)
    
    snowDP    <- (modH-height)/10 #cm
    ntime     <- length(datetime)
    snowDP_ma <- rep(NA, ntime)
    snowDP_ma[144:ntime] <- (modH-height_ma)/10 #cm
    plot(datetime,snowDP, type='l')
    lines(datetime,snowDP_ma, col=4)
    
    # caclulate monthly averages for snow depth and surface temps
    mo      <- format.Date(datetime, '%Y-%m')
    dat     <- cbind.data.frame(datetime,mo,snowDP,snowDP_ma,T5,T30)
    names(dat)
    meanT <- tapply(dat$T5, dat$mo,mean, na.rm=T)
    meanS <- tapply(dat$snowDP_ma, dat$mo,mean, na.rm=T)
    }  

  maxT <- length(height)
  ntime     <- length(datetime)

  #-----------------------------------------
  #--- now plot everything on one plot ----
  #-----------------------------------------
  fout <- paste(dir,"sample_node_",i,".pdf", sep='')
  pdf(fout, width=7, height=4)

  par(mfrow=c(3,1),mar=c(0,5,2,1), las=2)
  if (good != "Inf") {      #no data for #10, 13, 19
    plot(datetime,snowDP, type='l',xaxt="n",xlab="", ylab="Snow Depth (cm)",
         main=paste("node ",i,sep=""))
    lines(datetime,snowDP_ma, col=4)
    daterange=c(as.POSIXlt(min(datetime,na.rm=T)), as.POSIXlt(max(datetime,na.rm=T)))
    abline(h=0, lty=2)
  } else {
    plot(datetime,T5, type='n',xaxt="n",xlab="", ylab="Snow Depth (cm)",
         main=paste("node ",i,sep=""))
    daterange=c(as.POSIXlt(min(datetime,na.rm=T)), as.POSIXlt(max(datetime,na.rm=T)))
    abline(h=0, lty=2)
  }
  trange=c(-5,5) 
  par(mar=c(0,5,0,1))
  plot(datetime,T5, type='l',xaxt="n",xlab="", ylim=trange,col=2,
       ylab='Soil Temp.')
  lines(datetime,T30, lty=3, lwd=3, col=2)
  abline(h=0, lty=2)

  par(mar=c(5,5,0,1))
  mrange <- range(c(M5a,M30a),na.rm=T)
  plot(datetime,M5a, type='l',xaxt="n",xlab="", ylim=mrange, col=4,
       ylab='Soil Miosture')
  lines(datetime,M30a, lty=3, lwd=3, col=4)
  axis.POSIXct(1, at=seq(daterange[1], daterange[2], by="week"), format="%b %d")
  dev.next()

  par(mar = c(5,5,2,5), mfrow=c(1,1))
  if (good != "Inf") {      #no data for #10, 13, 19
    plot(datetime,snowDP, lty=2, lwd=0.5,
         type='l',xaxt="n",xlab="", 
         ylab="Snow Depth (cm)", 
        main=paste("node ",i,sep=""))
    lines(datetime,snowDP_ma, col=1, lwd=4)
    axis.POSIXct(1, at=seq(daterange[1], daterange[2], by="week"), format="%b %d")
    abline(h=0, lty=2)
  } else {
    plot(datetime,T5, lty=2, lwd=0.5,
         type='n',xaxt="n",xlab="", 
         ylab="Snow Depth (cm)", 
         main=paste("node ",i,sep=""))
    axis.POSIXct(1, at=seq(daterange[1], daterange[2], by="week"), format="%b %d")
    abline(h=0, lty=2)
  }
  par(new = T)
  plot(T5, type="l", lwd=1,col=2, 
       ylim=trange, axes=F, 
       xlab=NA, ylab=NA, cex=1.2)
  lines(T30, lty=3, col=2,lwd=3)
  abline(h=0, lty=2)
  axis(side = 4)
  mtext(side = 4, line = 2,las=0, 'Soil Temperature')
  dev.off()

  print(paste('finished ',i))

} # close big for loop

#-----------------------------------------
