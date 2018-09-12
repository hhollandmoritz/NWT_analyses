# Plots sample niwot data
# here data from new CR1000 loggers
# uploaded by Jen July 2018
# manuall added header to each column for ease of plotting
remove(list=ls())

help(package = CampbellLogger) 
#install.packages("devtools")
#install.packages("ForeCA")

library(devtools)
library(pracma)
library(ForeCA)
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
plots    <- seq(6,21,1)
nplots   <- length(plots)
sensor   <- c("VWC_5a","VWC_5b","VWC_5c","VWC_30a","VWC_30b","VWC_30c")
nsensors <- length(sensor)
omegaX   <- array(NA, c(nplots,nsensors), dimnames = list(plots, sensor))
omegaY   <- array(NA, c(nplots,nsensors), dimnames = list(plots, sensor))

#for (i in 7:7) {
for (i in 6:21) {
  if (i <  10) {Node <- paste('0',i ,sep='')}
  if (i >= 10) {Node <- paste(    i ,sep='')}
  fin  <- paste(dir,"sn_0",Node,"_cr1000x_tenminute.dat", sep='')
  dat  <- CampbellFileImport(fin, time.zone = "MST", skip.rows = 4)
# names(dat)  
# datetime  <- as.POSIXct(paste(dat$date, dat$time), format="%Y-%m-%d %H:%M:%S")
  datetime  <- dat$TIMESTAMP
  T5        <- dat$soiltemp_5cm_Avg
  T30       <- dat$soiltemp_30cm_Avg
  M5a       <- dat$vwca_5cm_Avg
  M30a      <- dat$vwca_30cm_Avg
  M5b       <- dat$vwcb_5cm_Avg
  M30b      <- dat$vwcb_30cm_Avg
  M5c       <- dat$vwcc_5cm_Avg
  M30c      <- dat$vwcc_30cm_Avg
  ntime     <- length(datetime)

# mask moisture when tsoi < 0 C
  M5A       <- M5a
  M5B       <- M5b
  M5C       <- M5c
  M30A      <- M30a
  M30B      <- M30b
  M30C      <- M30c
  M5A[T5 < 0] <- NA
  M5B[T5 < 0] <- NA
  M5C[T5 < 0] <- NA
  M30A[T30<0] <- NA
  M30B[T30<0] <- NA
  M30C[T30<0] <- NA
# Mask out 2017 data
  y       <- format.Date(datetime, '%Y')
  M5A[y == '2017'] <- NA
  M5B[y == '2017'] <- NA
  M5C[y == '2017'] <- NA
  M30A[y =='2017'] <- NA
  M30B[y =='2017'] <- NA
  M30C[y =='2017'] <- NA
  DATETIME <- datetime
  DATETIME [ is.na(M5A)] <- NA

  par(mfrow=c(1,1), mar=c(5,5,2,5))
  # no growing season data for #13, 17 19
  # values for 8 5a,5b make no sense
  if (i !=8 & i != 13 & i !=17 & i !=19) { 
  plot(M5A~DATETIME,  lwd=2, type='l',
       main=paste("node ",i,sep=""))
  lines(M5B~DATETIME, lwd=1.5)
  lines(M5C~DATETIME, lwd=1)
  
  max5A <- which.max(M5A)
  max5B <- which.max(M5B)
  max5C <- which.max(M5C)
  abline(v=datetime[max5A], lwd=2)
  abline(v=datetime[max5B], lwd=1.5)
  abline(v=datetime[max5C], lwd=1)
# mask out other points where soils likely frozen
  M5A [ index(datetime) < max5A & M5A < (max(M5A, na.rm=T) * 0.95)] <- NA
  M5B [ index(datetime) < max5B & M5B < (max(M5B, na.rm=T) * 0.95)] <- NA
  M5C [ index(datetime) < max5C & M5C < (max(M5C, na.rm=T) * 0.95)] <- NA
  lines(M5A ~ DATETIME, col=4 )
  lines(M5B ~ DATETIME, col=4 )
  lines(M5C ~ DATETIME, col=4 )

# Repeat for 30 cm Data  
  plot(M30A~DATETIME,  lwd=2, type='l', ylim=c(0.1,0.6),
       main=paste("node ",i,sep=""))
  lines(M30B~DATETIME, lwd=1.5)
  lines(M30C~DATETIME, lwd=1)
  
  max30A <- which.max(M30A)
  max30B <- which.max(M30B)
  max30C <- which.max(M30C)
  abline(v=datetime[max30A], lwd=2)
  abline(v=datetime[max30B], lwd=1.5)
  abline(v=datetime[max30C], lwd=1)
  M30A [ index(datetime) < max30A & M30A < (max(M30A, na.rm=T) * 0.95)] <- NA
  M30B [ index(datetime) < max30B & M30B < (max(M30B, na.rm=T) * 0.95)] <- NA
  M30C [ index(datetime) < max30C & M30C < (max(M30C, na.rm=T) * 0.95)] <- NA
  lines(M30A ~ DATETIME, col=4 )
  lines(M30B ~ DATETIME, col=4 )
  lines(M30C ~ DATETIME, col=4 )

  DATETIME2 <- datetime
  DATETIME2 [ is.na(M5A)] <- NA
  plot(na.omit(M30A))
  # approximate entropy calculated by pracma package
  # see https://stats.stackexchange.com/questions/126829/how-to-determine-forecastability-of-time-series
  
  all.series <- list(series1 = na.omit(M5A),
                     series2 = na.omit(M5B),
                     series3 = na.omit(M5C),
                     series4 = na.omit(M30A),
                     series5 = na.omit(M30B),
                     series6 = na.omit(M30C) )
# this is really slow and may take several minutes.
#  temp <- approx_entropy(na.omit(M30A), r=0.2*sd(na.omit(M30A)))
#  tempC <- approx_entropy(na.omit(M30C), r=0.2*sd(na.omit(M30C)))
#  print(paste(temp,tempC))

#  sample_A <- sample_entropy(na.omit(M30A),r=0.2*sd(na.omit(M30A)), tau = 48)
#  sample_C <- sample_entropy(na.omit(M30C),r=0.2*sd(na.omit(M30C)), tau = 48)
  # results not as expected, with lower entropy for M30C

  library(ForeCA)  
  # omega function [0,1], with 0 = white noise
  foreA <- Omega(na.omit(M30C), entropy.control = list(threshold = 1/40))
  foreB <- Omega(na.omit(M30B), spectrum.control = list(method = "wosa"))
  foreC <- Omega(na.omit(M30C), spectrum.control = list(method = "wosa"))
  X <- sapply(all.series,
              Omega, spectrum.control = list(method = "wosa"))
  omegaX[(i-5),] <- X  
  #gives lower score for 7_30c than 7_30a, which makes sense... but why?
  Y <- sapply(all.series,
              Omega, entropy.control = list(threshold = 1/40)) 
  omegaY[(i-5),] <- Y  
  
  entA <- spectral_entropy(na.omit(M30A), entropy.control = list(threshold = 1/40))
  
#  library(forecast)
  plot(M5A~DATETIME,  lwd=2, type='l', ylim=c(0.1,0.6),
       main=paste("node ",i,sep=""), ylab='VWC, 5 & 30')
  lines(M5B~DATETIME, lwd=1.5)
  lines(M5C~DATETIME, lwd=1)
  
  lines(M30A ~ DATETIME, col=4 )
  lines(M30B ~ DATETIME, col=4 )
  lines(M30C ~ DATETIME, col=4 )
  
  length(M5A[!is.na(M5A)])
  length(M5A)
  # caclulate monthly averages 
  mo      <- format.Date(datetime, '%Y-%m')
  dat     <- cbind.data.frame(datetime,mo,T5,T30,M5a,M5b,M5c,M30a,M30b,M30c)
  names(dat)
  meanT <- tapply(dat$T5, dat$mo,mean, na.rm=T)
  

  ntime     <- length(datetime)

  #-----------------------------------------
  #--- now plot everything on one plot ----
  #-----------------------------------------
  trange=c(-4,10) 
  
  fout <- paste(dir,"moisture_",i,".pdf", sep='')
  pdf(fout, width=7, height=4)

  par(mfrow=c(2,1),mar=c(0,5,3,5), las=2)
  plot(datetime,M5a, type='l',xaxt="n",xlab="", ylab="VWC, 5cm",
         main=paste("node ",i,sep=""), lwd=2)
  lines(datetime,M5b, lty=2, lwd=1.5)
  lines(datetime,M5c, lty=3, lwd=1)
  
  par(new = T)
  plot(T5, type="l", lwd=1,col=2, 
       ylim=trange, axes=F, 
       xlab=NA, ylab=NA, cex=1.2)
  abline(h=0, lty=2)
  axis(side = 4)
  mtext(side = 4, line = 2,las=0, 'Soil T, 5cm')
  
  par(mar=c(3,5,0,5), las=2)
  plot(datetime,M30a,  col=4, lty=1, lwd=2, type='l', ylab="VWC, 30cm")
  lines(datetime,M30b, col=4, lty=2, lwd=1.5)
  lines(datetime,M30c, col=4, lty=3, lwd=1)
  
  par(new = T)
  plot(T30, type="l", lwd=1,col=2, 
       ylim=trange, axes=F, 
       xlab=NA, ylab=NA, cex=1.2)
  abline(h=0, lty=2)
  axis(side = 4)
  mtext(side = 4, line = 2,las=0, 'Soil T, 30cm')
  dev.next()
  
  
# growing season only
  ylim=c(0.05,0.55)
  par(mfrow=c(2,1),mar=c(0,5,3,5), las=2)
  plot(DATETIME2,M5A, type='l',xaxt="n",xlab="", ylab="VWC, 5cm",
       main=paste("node ",i,sep=""), lwd=2, ylim=ylim)
  lines(DATETIME2,M5B, lty=2, lwd=1.5)
  lines(DATETIME2,M5C, lty=3, lwd=1)
  

  par(mar=c(3,5,0,5), las=2)
  plot(DATETIME2,M30A,  col=4, lty=1, lwd=2, type='l', 
       ylim=ylim,ylab="VWC, 30cm")
  lines(DATETIME2,M30B, col=4, lty=2, lwd=1.5)
  lines(DATETIME2,M30C, col=4, lty=3, lwd=1)
  dev.off()

  } # close if statement for #13
  print(paste('finished ',i))
} # close big for loop

# sensor at 15a not working
  omegaX[10,1] <- NA
  omegaY[10,1] <- NA
  write.table(omegaX,"omega_wosa.csv", sep=',')
  write.table(omegaY,"omega_ent1/40.csv", sep=',')
  boxplot(omegaX)
  boxplot(omegaY)
  #-----------------------------------------

  rotate <- function(x) t(apply(x, 2, rev))
  plot_omega <- rotate(omegaX)
boxplot(plot_omega)  
  