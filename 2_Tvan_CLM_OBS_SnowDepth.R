# Compares observed and simulated fluxes, soil temperature 

# Uses ReddyProc to gap fill Eddy-Flux data
# https://www.bgc-jena.mpg.de/bgi/index.php/Services/REddyProcWeb
# install.packages("raster", repos="http://R-Forge.R-project.org", type="source")

remove(list=ls())
library(REddyProc)
library(ncdf)
library(lattice)
library(lubridate)

#-----------------------------------------------------------------------
#--------------Snow Depth from saddle grid------------------------------
#-----------------------------------------------------------------------
# Plots community veg. data from Jane via Emily Farrer 
# uses long term community analysis 'class_3'
dir            <- ('/Users/wwieder/Desktop/Working_files/Niwot/NR_fluxes/TVan/Tvan_PTCLM45/')
data.veg       <- read.csv( paste(dir,'NWT_SnowXprod.csv',sep='') )
DATA.VEG       <- subset(data.veg,class_3!='rock' & class_3!='SF' & class_3!='ST',select=names(data.veg))
DATA.VEG       <- subset(DATA.VEG,DATA.VEG$year>=2007 & DATA.VEG$year<=2013,select=names(DATA.VEG))
good           <- list('FF','DM','MM','WM','SB')
DATA.VEG$class <- factor(DATA.VEG$class_3, levels=good)
VEG.PLOTS      <- as.numeric(levels(as.factor(DATA.VEG$plot)))
VEG.TYPE       <- DATA.VEG$class[1:length(VEG.PLOTS)]
  
# Read in full snow survey dataset
Data.sno <- read.csv( paste(dir,'saddsnow.dw.data.csv',sep='') )
DATA.SNO <- subset(Data.sno,Data.sno$YEAR>=2008 & Data.sno$YEAR<=2013,select=names(Data.sno))
names(DATA.SNO)
SAD_year  <- DATA.SNO$YEAR
SAD_mo    <- DATA.SNO$MO
SAD_date  <- DATA.SNO$date
SAD_depth <- DATA.SNO$mean_depth
SAD_plot  <- DATA.SNO$point_ID
SAD_TYPE  <- as.numeric(levels(SAD_plot))
SAD_veg   <- rep(NA, length(SAD_plot))
SAD_DATE  <- as.Date(as.character(SAD_date), format="%m/%d/%y")
obstime   <- levels(as.ordered(SAD_DATE))
ntime     <- length(obstime)
  
# assign vegetation classification
for (i in 1:length(VEG.PLOTS)) {
  SAD_veg[SAD_plot == VEG.PLOTS[i]] <- as.character(VEG.TYPE[i])
}
SAD_VEG <- factor(SAD_veg, levels=good)
  
names(DATA.VEG)
DATA.VEG$plot[1:78]
DATA.VEG$class
date_mean    <- tapply(SAD_DATE, INDEX=list(SAD_VEG, SAD_DATE), FUN=mean, na.rm=T)
depth_mean   <- tapply(SAD_depth, INDEX=list(SAD_VEG, SAD_DATE), FUN=mean, na.rm=T)
depth_sd     <- tapply(SAD_depth, INDEX=list(SAD_VEG, SAD_DATE), FUN=sd,   na.rm=T)

#-----------------------------------------------------------------------
#----Snow Depth from Tvan ----------------------------------------------
#-----------------------------------------------------------------------
# Plots community veg. data from Jane via Emily Farrer 
# uses long term community analysis 'class_3'
DATA.TVAN  <- read.csv( paste(dir,'SnowDepth_daily.csv', sep='') )
names(DATA.TVAN)
plot(DATA.TVAN$SnoDep_med)
TVAN_year <- as.character(DATA.TVAN$year_mean)
TVAN_mo   <- as.character(DATA.TVAN$mo_mean)
TVAN_day  <- as.character(DATA.TVAN$day_mean)
TVAN_depth<- DATA.TVAN$SnoDep_med
TVAN_date <- paste(TVAN_year,TVAN_mo,TVAN_day,sep="-")

TVAN_DATE  <- as.Date(TVAN_date)#, format="%m/%d/%y")
TVAN_depth[TVAN_year==2008] <- NA #obs seem poor

#-----------------------------------------------------------------------
#----Soil Moisture from Tvan -------------------------------------------
#-----------------------------------------------------------------------
Data.flx <- fLoadTXTIntoDataframe( paste(dir,'Tvan_flux_OBS_2008-2013b.txt',sep='') )
names(Data.flx)

SoiMoi   <- Data.flx$SoilMoisture/100
GPP_obs  <- Data.flx$GPP_f * 12 / (44 * 1e3) #convert from mgCO2/m2/s to gC/m2/s
NEE_obs  <- Data.flx$NEE   * 12 / (44 * 1e3) #convert from mgCO2/m2/s to gC/m2/s
y_obs    <- as.character(Data.flx$Year)
m_obs    <- as.character(Data.flx$MO)
d_obs    <- as.character(Data.flx$DD)
FLX_date <- paste(y_obs,m_obs,d_obs,sep="-")
FLX_DATE <- as.Date(FLX_date)

SOIMOI   <- tapply(SoiMoi , FLX_DATE, mean, na.rm=T)
GPP_OBS  <- tapply(GPP_obs, FLX_DATE, mean, na.rm=T) * 3600 * 24    #gC/m2/d
NEE_OBS  <- tapply(NEE_obs, FLX_DATE, mean, na.rm=T) * 3600 * 24
FLX_DATE <- as.Date(levels(as.factor(FLX_DATE)))

Data.DM <-  read.csv('/Users/wwieder/Desktop/Working_files/Niwot/NR_fluxes/TVan/Soil_moisture/sdlcr23x-cr1000.hourly.ml.data.csv')
Data.MM <-  read.csv('/Users/wwieder/Desktop/Working_files/Niwot/NR_fluxes/TVan/Soil_moisture/SoilMoisture_MM_14.csv')
Data.WM <-  read.csv('/Users/wwieder/Desktop/Working_files/Niwot/NR_fluxes/TVan/Soil_moisture/SoilMoisture_WM_18.csv')

names(Data.DM)
Data.DM$soiltemp_avg[Data.DM$soiltemp_avg =='NaN'] <- NA
Data.DM$vwc_avg[     Data.DM$vwc_avg      =='NaN'] <- NA
Data.DM$vwc_10cm_avg[Data.DM$vwc_10cm_avg =='NaN'] <- NA
Data.DM$vwc_20cm_avg[Data.DM$vwc_20cm_avg =='NaN'] <- NA
DATE_DM      <- as.Date(Data.DM$date.time, "%Y-%m-%d %H:%M")
SOIMOI_DM    <- tapply(Data.DM$vwc_avg, DATE_DM, mean, na.rm=T)
TSOI_DM      <- tapply(Data.DM$soiltemp_avg, DATE_DM, mean, na.rm=T)
DATE_DM      <- as.Date(levels(as.factor(DATE_DM)))
plot(SOIMOI_DM~DATE_DM)

names(Data.MM)
Data.MM$SoilMoist[Data.MM$SoilMoist=='NaN'] <- NA
DATE_MM   <- as.Date(Data.MM$date.time, "%m/%d/%y %H:%M")
SOIMOI_MM <- tapply(Data.MM$SoilMoist, DATE_MM, mean, na.rm=T)
DATE_MM   <- as.Date(levels(as.factor(DATE_MM)))
plot(SOIMOI_MM~DATE_MM)

names(Data.WM)
Data.WM$SoilMoist[Data.WM$SoilMoist=='NaN'] <- NA
DATE_WM   <- as.Date(Data.WM$date.time, "%m/%d/%y %H:%M")
SOIMOI_WM <- tapply(Data.WM$SoilMoist, DATE_WM, mean, na.rm=T)
DATE_WM   <- as.Date(levels(as.factor(DATE_WM)))
plot(SOIMOI_WM~DATE_WM)

#-----------------------------------------------------------------------
#----Soil Temperature from Tvan -------------------------------------------
#-----------------------------------------------------------------------
Temp.Path <- '/Users/wwieder/Desktop/Working_files/Niwot/NR_fluxes/TVan/Soil_temp/'
Temp.FF   <- read.csv(paste(Temp.Path, 'SoilTemp_FF_2.csv' , sep=''))
Temp.DM   <- read.csv(paste(Temp.Path, 'SoilTemp_DM_13.csv', sep=''))
Temp.MM   <- read.csv(paste(Temp.Path, 'SoilTemp_MM_14.csv', sep=''))
Temp.WM   <- read.csv(paste(Temp.Path, 'SoilTemp_WM_18.csv', sep=''))

Temp.FF$Temp_avg[Temp.FF$Temp_avg =='NaN']  <- NA
Temp.DM$Temp_avg[Temp.DM$Temp_avg =='NaN']  <- NA
Temp.MM$Temp_avg[Temp.MM$Temp_avg =='NaN']  <- NA
Temp.WM$Temp_avg1[Temp.WM$Temp_avg1 =='NaN']<- NA
Temp.WM$Temp_avg2[Temp.WM$Temp_avg2 =='NaN']<- NA
Temp.WM$Temp_avg <- rowMeans(cbind(Temp.WM$Temp_avg1,Temp.WM$Temp_avg2), na.rm=T)

TEMP_DATE_FF <- as.Date(Temp.FF$Time, "%m/%d/%y %H:%M")
TEMP_DATE_DM <- as.Date(Temp.DM$Time, "%m/%d/%y %H:%M")
TEMP_DATE_MM <- as.Date(Temp.MM$Time, "%m/%d/%y %H:%M")
TEMP_DATE_WM <- as.Date(Temp.WM$Time1,"%m/%d/%y %H:%M")

TSOI_FF      <- tapply(Temp.FF$Temp_avg, TEMP_DATE_FF, mean, na.rm=T)
TSOI_DM2     <- tapply(Temp.DM$Temp_avg, TEMP_DATE_DM, mean, na.rm=T)
TSOI_MM      <- tapply(Temp.MM$Temp_avg, TEMP_DATE_MM, mean, na.rm=T)
TSOI_WM      <- tapply(Temp.WM$Temp_avg, TEMP_DATE_WM, mean, na.rm=T)

TEMP_DATE_FF <- as.Date(levels(as.factor(TEMP_DATE_FF)))
TEMP_DATE_DM <- as.Date(levels(as.factor(TEMP_DATE_DM)))
TEMP_DATE_MM <- as.Date(levels(as.factor(TEMP_DATE_MM)))
TEMP_DATE_WM <- as.Date(levels(as.factor(TEMP_DATE_WM)))

plot( TEMP_DATE_FF, TSOI_FF,  type="l")
lines(TEMP_DATE_DM, TSOI_DM2, col=2)
lines(TEMP_DATE_MM, TSOI_MM,  col=3)
lines(TEMP_DATE_WM, TSOI_WM,  col=4)

#-----------------------------------------------------------------------
#---------------read in CLM variables----------------------------------
#-----------------------------------------------------------------------
#path   <- paste(dir,'CLM_nc_files/bgc_clm5phys/',sep='')
#case   <- ('Tvan_PTCLM50bgc_allPHYS_noLW_lowVCmax')
#suf <- ('_run.clm2.h1.2008-01-01-00000.nc')

# --- to use PTrespmods --- 
path     <- paste(dir,'CLM_nc_files/PTrespmods/lowWATSAT/',sep='')
pre      <- ('Tvan_PTrespmods_allPHYS_noLW_lowVCmax_')
suf      <- ('_lowRESIST_run.clm2.h1.2008-01-01-00000.nc')
case     <- c('010dry_ff_noRHpsn_lowWATSAT','010dry_dm_noRHpsn_lowWATSAT',
              '100_mm_noRHpsn','075_wm_noRHpsn',
              '200_sb_noRHpsn_lowWATSAT')
#case     <- c('010dry','010dry_dm','100','075_wm','200')
#case     <- c('015dry','0150dry_dm','100','075_wm','200')
nrows  <- length(case)
nyears <- 6
ndays  <- 365 * nyears + 2
nlevsoi     <- 25
dims            <- c(nrows,ndays)
dims2           <- c(nrows,nlevsoi,ndays)
DATE_mean       <- array(NA,dim=dims)
GPP_mean        <- array(NA,dim=dims)
NPP_mean        <- array(NA,dim=dims)
HR_mean         <- array(NA,dim=dims)
NEE_mean        <- array(NA,dim=dims)
SNOW_mean       <- array(NA,dim=dims)
RAIN_mean       <- array(NA,dim=dims)
BTRAN_mean      <- array(NA,dim=dims)
FPG_mean        <- array(NA,dim=dims)
PPT_sum         <- array(NA,dim=dims)
SNOWICE_mean    <- array(NA,dim=dims)
SNOW_DEPTH_mean <- array(NA,dim=dims)
H2OSOI_mean     <- array(NA,dim=dims2)       
TSOI_mean       <- array(NA,dim=dims2)
SOILLIQ_mean    <- array(NA,dim=dims2)

for (e in 1:length(case)) {  
  infile   <- paste(path,pre,case[e],suf,sep='')
  Data.clm <- open.ncdf(infile)   
  print(paste('read',infile))  
  #print(Data.clm)
  levgrnd     <- get.var.ncdf(Data.clm,'levgrnd')
  MCDATE      <- get.var.ncdf(Data.clm, "mcdate") # getting/naming the variable
  MCSEC       <- get.var.ncdf(Data.clm, "mcsec") 
  GPP         <- get.var.ncdf(Data.clm, "GPP")          #GPP (gC/m2/s)
  NPP         <- get.var.ncdf(Data.clm, "NPP")          #GPP (gC/m2/s)
  HR          <- get.var.ncdf(Data.clm, "HR")           #HR (gC/m2/s)
  NEE         <- get.var.ncdf(Data.clm, "NEE")          #NEE (gC/m2/s)
  TSOI        <- get.var.ncdf(Data.clm, "TSOI")         #soil temperature (K)
  H2OSOI      <- get.var.ncdf(Data.clm, "H2OSOI")	     #Volumetric soil moisture
  SNOW   		  <- get.var.ncdf(Data.clm, "SNOW")       
  RAIN     	  <- get.var.ncdf(Data.clm, "RAIN")       
  BTRAN       <- get.var.ncdf(Data.clm, "BTRAN")       
  FPG         <- get.var.ncdf(Data.clm, "FPG")       
  SOILLIQ	    <- get.var.ncdf(Data.clm, "SOILLIQ")      #soil liquid water 
  SNOWICE	    <- get.var.ncdf(Data.clm, "SNOWICE")      #snow ice water 
  SNOW_DEPTH  <- get.var.ncdf(Data.clm, "SNOW_DEPTH")	 #snow height of snow covered area
  H2OSOIunits <- att.get.ncdf(Data.clm, "H2OSOI", "units")
  H2OSOIunits
  mean(H2OSOI)
  TSOI   <- TSOI - 273.15
  BTRAN[GPP<=0] <- NA      #so only + values count in the calculation
  FPG[  GPP<=0] <- NA
  nsteps <- length(MCDATE)
  MCDATE[nsteps]
  MCSEC[2]
  MCDATE[nsteps-48]
  n1 <- 365*48  #for no-leap   data
  n2 <- n1 + 48 #for leap year data (leap)
  
  shift      <- 15  #adjust CLM to MST  
  mcdate     <- MCDATE[1:(nsteps-shift+1)]
  mcsec      <- MCSEC[1:(nsteps-shift+1)]
  mchour     <- mcsec/3600
  nee        <- NEE[1+shift:nsteps-1]   
  gpp        <- GPP[1+shift:nsteps-1]   
  npp        <- NPP[1+shift:nsteps-1]   
  hr         <- HR[ 1+shift:nsteps-1]   
  snow       <- SNOW[1+shift:nsteps-1]
  rain       <- RAIN[1+shift:nsteps-1]
  btran      <- BTRAN[1+shift:nsteps-1]
  fpg        <- FPG[1+shift:nsteps-1]
  snowice    <- SNOWICE[1+shift:nsteps-1]
  snow_depth <- SNOW_DEPTH[1+shift:nsteps-1]
  h2osoi     <- H2OSOI[,1+shift:nsteps-1]	
  tsoi       <- TSOI[,1+shift:nsteps-1]   
  soilliq    <- SOILLIQ[, 1+shift:nsteps-1]
  nstep2     <- length(mcdate)
  precip     <- rain + snow
  dim(h2osoi)
  mcdate[1]
  mcdate[nstep2]
  day <- as.factor(mcdate)
  
  DATE_mean[e,]    <- tapply(mcdate, day, mean)
  GPP_mean[e,]     <- tapply(gpp,    day, mean) * 3600 * 24 
  NPP_mean[e,]     <- tapply(npp,    day, mean) * 3600 * 24 
  HR_mean[e,]      <- tapply(hr ,    day, mean) * 3600 * 24 
  NEE_mean[e,]     <- tapply(nee,    day, mean) * 3600 * 24
  SNOW_mean[e,]    <- tapply(snow,   day, mean) * 3600 * 24
  RAIN_mean[e,]    <- tapply(rain,   day, mean) * 3600 * 24
  BTRAN_mean[e,]   <- tapply(btran,  day, mean)
  FPG_mean[e,]     <- tapply(fpg,  day, mean)
  SNOWICE_mean[e,] <- tapply(snowice,day, mean)
  SNOW_DEPTH_mean[e,] <- tapply(snow_depth,day, mean)
  PPT_sum[e,]      <- SNOW_mean[e,] + RAIN_mean[e,]
#  SWE_mean[e,]     <- SNOWLIQ_mean[e,] + SNOWICE_mean[e,]
  for (i in 1:nlevsoi) {
    SOILLIQ_mean[e,i,] <- tapply(soilliq[i,],day, mean)       
    H2OSOI_mean[e,i,]  <- tapply(h2osoi[i,], day, mean)       
    TSOI_mean[e,i,]    <- tapply(tsoi[i,],   day, mean)       
  }

  close.ncdf(Data.clm)  
}
dim(GPP_mean)
rowMeans(GPP_mean)  * 365
rowMeans(NPP_mean)  * 365
rowMeans(SNOW_mean) * 365 
rowMeans(RAIN_mean) * 365  
length(HR_mean[1,])
plot(HR_mean[1,(200:500)], type='l')
  for (e in 1:5) { lines(HR_mean[e,(200:500)], col=e)}
abline(h=0)
dim(TSOI_mean)

# look at ranges in soil temperature
for (e in 1:5) {
  rangeTSOI <- range(TSOI_mean[e,3,], na.rm=T)
  print(paste('CLM ',levels(VEG.TYPE)[e],'=',sum(abs(rangeTSOI))[[1]],',',rangeTSOI[[1]], ',',rangeTSOI[[2]],sep=''))
}

print(paste(sum(abs(range(TSOI_FF, na.rm=T)))))
print(paste(sum(abs(range(TSOI_DM, na.rm=T)))))
print(paste(sum(abs(range(TSOI_DM2, na.rm=T)))))
print(paste(sum(abs(range(TSOI_MM, na.rm=T)))))
print(paste(sum(abs(range(TSOI_WM, na.rm=T)))))
#---------------------------------------------------------
#---------------------------------------------------------
#   Plot model vs. obs.
#---------------------------------------------------------
#---------------------------------------------------------
s2y       <- 3600 * 24 * 365
ncase     <- length(case)
year      <- seq(2007,2013,1) 
nyear     <- length(year)
OBS_TIME  <- as.Date(obstime)
OBS_YEAR  <- as.numeric(format.Date(OBS_TIME, "%Y"))
OBS_MONTH <- as.numeric(format.Date(OBS_TIME, "%m"))
OBS_DATE  <- yday(OBS_TIME)
at        <- as.Date(c('2008-01-01','2009-01-01','2010-01-01','2011-01-01','2012-01-01','2013-01-01'))

CLM_date <- as.Date(as.character(DATE_mean[1,]),format="%Y%m%d")
CLM_year <- as.numeric(format.Date(CLM_date, "%Y"))
CLM_month<- as.numeric(format.Date(CLM_date, "%m"))
CLM_dm   <- format.Date(CLM_date, "%m%d")
CLM_yday <- yday(CLM_date)
Ylim     <- c(min(SOIMOI,na.rm=T), max(H2OSOI_mean,na.rm=T))
plot(H2OSOI_mean[1,3,] ~ CLM_date, type="l", ylim=Ylim, lwd=3, col=2)
OBS_yday <- yday(OBS_TIME)   # julian day of observations
FLX_yday <- yday(FLX_DATE)
FLX_year <- as.numeric(format.Date(FLX_DATE, "%Y"))

# annual totals for each case
RAIN_SUM   <- array(NA, c(ncase,nyear))
SNOW_SUM   <- array(NA, c(ncase,nyear))
GPP_SUM    <- array(NA, c(ncase,nyear))
BTRAN_MEAN <- array(NA, c(ncase,nyear))
FPG_MEAN   <- array(NA, c(ncase,nyear))
for ( e in 1:ncase) {
  for ( y in 2:nyear) {
    GPP_temp      <- GPP_mean[e,] 
    SNOW_temp     <- SNOW_mean[e,]
    RAIN_temp     <- RAIN_mean[e,]
    BTRAN_temp    <- BTRAN_mean[e,]
    FPG_temp      <- FPG_mean[e,]
    GPP_SUM[e,y]  <- sum(GPP_temp[CLM_year == year[y]]) 
    SNOW_SUM[e,y] <- sum(SNOW_temp[CLM_year == year[y]])
    RAIN_SUM[e,y] <- sum(RAIN_temp[CLM_year == year[y]])
    BTRAN_MEAN[e,y] <- mean(BTRAN_temp[CLM_year == year[y] & GPP_temp > 0])
    FPG_MEAN[e,y] <- mean(FPG_temp[CLM_year == year[y] & GPP_temp > 0])
  }
}

temp <- transform(RAIN_SUM, Mean=apply(RAIN_SUM,1, mean, na.rm = TRUE),
                  SD=apply(RAIN_SUM,1, sd, na.rm = TRUE))
cbind(temp$Mean, temp$SD)

temp <- transform(SNOW_SUM, Mean=apply(SNOW_SUM,1, mean, na.rm = TRUE),
                  SD=apply(SNOW_SUM,1, sd, na.rm = TRUE))
cbind(temp$Mean, temp$SD)

temp <- transform(BTRAN_MEAN, Mean=apply(BTRAN_MEAN,1, mean, na.rm = TRUE),
                  SD=apply(BTRAN_MEAN,1, sd, na.rm = TRUE))
cbind(temp$Mean, temp$SD)

temp <- transform(FPG_MEAN, Mean=apply(FPG_MEAN,1, mean, na.rm = TRUE),
                  SD=apply(FPG_MEAN,1, sd, na.rm = TRUE))
cbind(temp$Mean, temp$SD)

ann_SNOW_DEPTH_mean <- array(NA,c(ncase,366)) 
ann_SNOW_DEPTH_sd   <- array(NA,c(ncase,366)) 
ann_GPP_mean        <- array(NA,c(ncase,366)) 
ann_GPP_sd          <- array(NA,c(ncase,366)) 
ann_BTRAN_mean      <- array(NA,c(ncase,366)) 
ann_BTRAN_sd        <- array(NA,c(ncase,366)) 
ann_FPG_mean      <- array(NA,c(ncase,366)) 
ann_FPG_sd        <- array(NA,c(ncase,366)) 
ann_dday            <- CLM_yday[1:366]

for (e in 1:ncase) {
  ann_SNOW_DEPTH_mean[e,] <- tapply(SNOW_DEPTH_mean[e,]*100,CLM_dm,mean)
  ann_SNOW_DEPTH_sd[e,]   <- tapply(SNOW_DEPTH_mean[e,]*100,CLM_dm,sd)
  ann_GPP_mean[e,]        <- tapply(GPP_mean[e,],CLM_dm,mean)
  ann_GPP_sd[e,]          <- tapply(GPP_mean[e,],CLM_dm,sd)
  ann_FPG_mean[e,]        <- tapply(FPG_mean[e,],CLM_dm,mean)
  ann_FPG_sd[e,]          <- tapply(FPG_mean[e,],CLM_dm,sd)
  ann_BTRAN_mean[e,]      <- tapply(BTRAN_mean[e,],CLM_dm,mean)
  ann_BTRAN_sd[e,]        <- tapply(BTRAN_mean[e,],CLM_dm,sd)
  ann_FPG_mean[e,]      <- tapply(FPG_mean[e,],CLM_dm,mean)
  ann_FPG_sd[e,]        <- tapply(FPG_mean[e,],CLM_dm,sd)
}

ann_BTRAN_mean[ann_GPP_mean == 0] <- NA
ann_BTRAN_sd[ann_GPP_mean == 0] <- NA
ann_FPG_mean[ann_GPP_mean == 0] <- NA
ann_FPG_sd[ann_GPP_mean == 0] <- NA

plot(ann_BTRAN_mean[5,], type='l', ylim=c(0,1))
lines(ann_FPG_mean[5,], col=1, lty=3)

lines(ann_BTRAN_mean[3,], col=3)
lines(ann_BTRAN_mean[4,], col=4)
lines(ann_BTRAN_mean[5,], col=5)
dim(ann_SNOW_DEPTH_mean)

fout <- paste(path,'GPP_',case[1],'SnowDepth.pdf', sep='')
pdf(fout, width=6, height=6, compress = FALSE)

par(mfrow=c(5,1), mar=c(0,0,0,0), oma=c(5,5,3,2))
for (e in 1:length(case)) {
  if ( e <= 2) { ylim=c(0,160) }
  if ( e >= 3) { ylim=c(0,275) }
  if ( e == 5) { ylim=c(0,450)}
  plot(ann_SNOW_DEPTH_mean[e,] ~ ann_dday, 
       xlab=NA, ylab=NA, xaxt='n',
       col=2, lwd=3, type='l', ylim=ylim) 
  xx1 <- c(ann_dday, rev(ann_dday))
  yy1 <- c((   ann_SNOW_DEPTH_mean[e,] + ann_SNOW_DEPTH_sd[e,]), 
           rev(ann_SNOW_DEPTH_mean[e,] - ann_SNOW_DEPTH_sd[e,]))
  polygon(xx1, yy1, col = rgb(1, 0, 0,0.35), border = NA)
  at  <- seq(0,350,50)
  axis(1, tck = 0.05, labels = FALSE, at=at)
  text(x=20,y=ylim[2]*0.9,labels=good[e],cex=2)
  
  for (y in 2:nyear) {
    lines(SNOW_DEPTH_mean[e,CLM_year==year[y]]*100 ~ CLM_yday[CLM_year==year[y]], col=2, lwd=1)
    if (e == 3 ) {mtext("Snow Depth (cm)", side = 2, line = 3)}
    if (e == 5 ) {
      mtext("Julian Day", side = 1, line = 3)
      axis(1, tck = 0.05, at=at, outer=TRUE)
    }

  }  # close year loop
  points(depth_mean[e,] ~ OBS_yday , cex=1.3,pch=4, lwd=2)
  }    # close case loop
dev.next()

par(mfrow=c(5,1), mar=c(0,0,0,0), oma=c(5,2,3,5))
for (e in 1:length(case)) {
  plot(ann_GPP_mean[e,] ~ ann_dday, 
       xlab=NA, ylab=NA, xaxt='n', yaxt='n',
       col=2, lwd=3, type='l', ylim=range(GPP_mean)) 
  xx1 <- c(ann_dday, rev(ann_dday))
  yy1 <- c((   ann_GPP_mean[e,] + ann_GPP_sd[e,]), 
           rev(ann_GPP_mean[e,] - ann_GPP_sd[e,]))
  polygon(xx1, yy1, col = rgb(1, 0, 0,0.35), border = NA)
  at  <- seq(0,350,50)
  axis(1, tck = 0.05, labels = FALSE, at=at)
  axis(4, at=seq(0,9,3))
  for (y in 2:nyear) {
    lines(GPP_mean[e,CLM_year==year[y]] ~ CLM_yday[CLM_year==year[y]], col=2, lwd=1)
    if (e == 3 ) {
      mtext(expression(paste("GPP (gC ",m^-2," ",d^-1,")",sep="")), 
            adj = 0, side = 4, line = 3)
    }
    if (e == 5 ) {
      mtext("Julian Day", side = 1, line = 3)
      axis(1, tck = 0.05, at=at, outer=TRUE)
    }
    if ( e == 1 ) {
      lines(GPP_OBS[FLX_year==year[y]] ~ FLX_yday[FLX_year==year[y]], col=1, lwd=1, type='l') 
    }
  }  # close year loop
}    # close case loop
dev.off()

print(paste('wrote',fout,sep=" "))

#------------------------------------
#------------------------------------


# create arrays to store max snow depth & snow off date
max_OBS_depth <- array(NA, c(ncase,nyears)) 
max_CLM_depth <- array(NA, c(ncase,nyears)) 
max_OBS_day   <- array(NA, c(ncase,nyears)) 
max_CLM_day   <- array(NA, c(ncase,nyears)) 
min_OBS_depth <- array(NA, c(ncase,nyears)) 
min_CLM_depth <- array(NA, c(ncase,nyears)) 
min_OBS_day   <- array(NA, c(ncase,nyears)) 
min_CLM_day   <- array(NA, c(ncase,nyears)) 
ma <- function(x,n){filter(x,rep(1/n,n), sides=2, circular = TRUE)}

for (e in 1:ncase) {
  for (y in 1:nyears) {
    temp_obs_depth <- depth_mean[e,]
    temp_obs_depth <- temp_obs_depth[OBS_YEAR == year[(y+1)]] 
    temp_obs_date  <- OBS_DATE[ OBS_YEAR == year[(y+1)]]
    temp_obs_month <- OBS_MONTH[OBS_YEAR == year[(y+1)]]
    
    temp_clm_depth <- SNOW_DEPTH_mean[e,]*100
    temp_clm_depth <- temp_clm_depth[CLM_year == year[(y+1)]] 
    ma_clm_depth   <- ma(temp_clm_depth,11)
    temp_clm_date  <- CLM_yday[CLM_year == year[(y+1)]]
    temp_clm_month <- CLM_month[CLM_year == year[(y+1)]]

    max_OBS_depth[e,y] <- max(temp_obs_depth)
    max_CLM_depth[e,y] <- max(ma_clm_depth[temp_clm_month < 8], na.rm=T)
    max_OBS_day[e,y] <- temp_obs_date[temp_obs_depth == max_OBS_depth[e,y]]
    max_CLM_day[e,y] <- max(temp_clm_date[ma_clm_depth == max_CLM_depth[e,y]], na.rm=T)
    
    min_OBS_depth[e,y] <- min(temp_obs_depth[temp_obs_date > max_OBS_day[e,y] &
                                               temp_obs_month < 9])
    min_CLM_depth[e,y] <- min(ma_clm_depth[temp_clm_date > max_CLM_day[e,y] & 
                                             temp_clm_month < 9], na.rm=T)
    min_OBS_day[e,y] <- min(temp_obs_date[temp_obs_date > max_OBS_day[e,y] &
                                            temp_obs_depth == min_OBS_depth[e,y]], na.rm=T)
    min_CLM_day[e,y] <- min(temp_clm_date[ma_clm_depth == min_CLM_depth[e,y]], na.rm=T)
    
    plot(temp_clm_depth~temp_clm_date, type='l', main=paste(case[e],year[(y+1)]))
    lines(ma_clm_depth~temp_clm_date, col=2)
    points(temp_obs_depth~temp_obs_date)
    points(max_OBS_depth[e,y]~max_OBS_day[e,y], pch='X')
    points(min_OBS_depth[e,y]~min_OBS_day[e,y], pch='X')
    
  }  
}
ylim = range(c(max_OBS_depth,max_CLM_depth)) * c(0.9,1.1)
xlim = range(c(min_OBS_day,min_CLM_day)) * c(0.9,1.1)

fout <- paste(path,'FigX_',case[1],'meanSnowDepth.pdf', sep='')
pdf(fout, width=9, height=7, compress = FALSE)
par(mfrow=c(2,2), mar=c(2,2,2,2), oma=c(3,3,1,1))
plot(max_CLM_depth[1,]~min_CLM_day[1,], pch=0, 
     ylim=ylim, xlim=xlim, 
     ylab='Max Snow Depth (mm)',
     xlab='Melt out (Julian Day)')
mtext('max_depth (cm)', 2, line=2)
mtext('melt out (d)', 1, line=2)

for (e in 1:ncase){
  points(max_CLM_depth[e,]~min_CLM_day[e,], pch=e, col=2 ) 
  points(max_OBS_depth[e,]~min_OBS_day[e,], pch=e, col=1 ) 
}
rowMeans(max_CLM_depth)
rowMeans(max_OBS_depth)
rowMeans(min_CLM_day)
rowMeans(min_OBS_day) 

obs_depth   <- apply(max_OBS_depth,1, mean, na.rm = TRUE)
clm_depth   <- apply(max_CLM_depth,1, mean, na.rm = TRUE)
obs_depthSD <- apply(max_OBS_depth,1, sd  , na.rm = TRUE)
clm_depthSD <- apply(max_CLM_depth,1, sd  , na.rm = TRUE)

obs_day   <- apply(min_OBS_day,1, mean, na.rm = TRUE)
clm_day   <- apply(min_CLM_day,1, mean, na.rm = TRUE)
obs_daySD <- apply(min_OBS_day,1, sd,   na.rm = TRUE)
clm_daySD <- apply(min_CLM_day,1, sd,   na.rm = TRUE)

lim=c(0,375)
pch=levels(VEG.TYPE)
plot(obs_depth~clm_depth, ylim=lim, xlim=lim, pch=16, cex=1.8)
arrows(clm_depth, obs_depth-obs_depthSD, clm_depth,obs_depth+obs_depthSD, angle=90, length=0.1, lwd=2, code=3)
arrows(clm_depth-clm_depthSD, obs_depth, clm_depth+clm_depthSD,obs_depth, angle=90, length=0.1, lwd=2, code=3)
abline(0,1, lty=2)
mtext('obs max depth (cm)', 2, line=2)
mtext('clm max depth (cm)', 1, line=2)

lim=c(100,225)
pch=levels(VEG.TYPE)
plot(obs_day~clm_day, ylim=lim, xlim=lim, pch=16, cex=1.8)
arrows(clm_day, obs_day-obs_daySD, clm_day,obs_day+obs_daySD, angle=90, length=0.1, lwd=2, code=3)
arrows(clm_day-clm_daySD, obs_day, clm_day+clm_daySD,obs_day, angle=90, length=0.1, lwd=2, code=3)
abline(0,1, lty=2)
mtext('obs melt out (d)', 2, line=2)
mtext('clm melt out (d)', 1, line=2)

dev.off()

temp <- transform(min_CLM_day, Mean=apply(min_CLM_day,1, mean, na.rm = TRUE),
          SD=apply(min_CLM_day,1, sd, na.rm = TRUE))
cbind(temp$Mean, temp$SD)

temp <- transform(min_OBS_day, Mean=apply(min_OBS_day,1, mean, na.rm = TRUE),
          SD=apply(min_OBS_day,1, sd, na.rm = TRUE))
cbind(temp$Mean, temp$SD)

temp <- transform(max_OBS_depth, Mean=apply(max_OBS_depth,1, mean, na.rm = TRUE),
                  SD=apply(max_OBS_depth,1, sd, na.rm = TRUE))
cbind(temp$Mean, temp$SD)

temp <- transform(max_CLM_depth, Mean=apply(max_CLM_depth,1, mean, na.rm = TRUE),
                  SD=apply(max_CLM_depth,1, sd, na.rm = TRUE))
cbind(temp$Mean, temp$SD)


# --- Plots of full inter annual variability ---
fout <- paste(path,'Fig3_Tvan_full_SnowDepth.pdf', sep='')
pdf(fout, width=9, height=7, compress = FALSE)
par(mfrow=c(5,1), mar=c(0,0,0,0), oma=c(5,5,3,2))

at       <- as.Date(c('2008-01-01','2009-01-01','2010-01-01',
                      '2011-01-01','2012-01-01','2013-01-01'))
for (e in 1:ncase) {
  CLM_date <- as.Date(as.character(DATE_mean[e,]),format="%Y%m%d")
  ylim     <- c(0, max(depth_mean+depth_sd))
  xx1 <- c(OBS_TIME, rev(OBS_TIME))
  yy1 <- c((depth_mean[e,] + depth_sd[e,]), rev(depth_mean[e,] - depth_sd[e,]))
  Ylim=c(0,400)

  plot(depth_mean[e,] ~ OBS_TIME, type="l", ylim=Ylim, lwd=3, 
       xlab=NA, ylab=NA, xaxt='n')
  polygon(xx1, yy1, col = rgb(0, 0, 0,0.35), border = NA)
#  plot moving average of simulations
  lines(ma((SNOW_DEPTH_mean[e,]*100),11) ~ CLM_date, col=2, lwd=2, type='l') 
  points(depth_mean[e,] ~ OBS_TIME, cex=1.2,pch=4)
  text(x=CLM_date[30],y=Ylim[2]*0.9,labels=good[e],cex=2)
  if (e == 3 ) {mtext("Snow Depth (cm)", side = 2, line = 3)}
  if (e == 5 ) {
    mtext("Year", side = 1, line = 3)
    axis(1, tck = 0.05, labels=year[2:7], at=at)
  }
  if (e == 1) {
    lines(TVAN_depth ~ TVAN_DATE, col=4, lwd=1, type='l') 
    legend('topright',c("CLM4.5","Saddle Obs.","TVan Obs."), 
           col=c(2,1,4), lty=1, lwd=3, bty="n")
  }  
}           # close case loop

dev.off()

# ---- soil moisture ---- 
fout <- paste(path,'FigSI_1_Tvan_SoilMoisture.pdf', sep='')
pdf(fout, width=9, height=7, compress = FALSE)
par(mfrow=c(5,1), mar=c(0,0,0,0), oma=c(5,5,3,2))

for (e in 1:length(case)) {
  CLM_date <- as.Date(as.character(DATE_mean[e,]),format="%Y%m%d")
  Ylim     <- c(min(SOIMOI,na.rm=T), max(0.6,na.rm=T))
  tempY    <- H2OSOI_mean[e,,]
  tempY[TSOI_mean[e, ,] < 0]  <- NA 
  tempX    <- CLM_date
  tempX[TSOI_mean[e,3,] < 0 ] <- NA
  plot(tempY[3,] ~ CLM_date, type="l", ylim=Ylim, lwd=3, col=2, 
       xlab=NA, ylab=NA, xaxt='n') 
  text(x=CLM_date[30],y=Ylim[2]*0.9,labels=good[e],cex=2)
  if (e == 3) {
    mtext(expression(paste("Soil Moisture (",cm^3," ",cm^-3,")",sep="")), 
          side = 2, line = 3)
  }
  if (e == 5) {
    mtext("Year", side = 1, line = 3)
    axis(1, tck = 0.05, labels=year[2:7], at=at)
  }
  if (e == 1) {
    lines(SOIMOI ~ FLX_DATE, col=4, lwd=1) 
    legend('topright',c("CLM4.5","Saddle","Tvan"), 
           col=c(2,1,4), lty=1, lwd=3, bty="n")
  }  
  if (e == 2 ) {
    lines(SOIMOI_DM ~ DATE_DM, col=1, lwd=1) 
  }  
  if (e == 3) {
    lines(SOIMOI_MM ~ DATE_MM, col=4, lwd=1) 
  }  
  if (e == 4) {
    lines(SOIMOI_WM ~ DATE_WM, col=4, lwd=1) 
  }  
}           # close case loop
dev.off()


#  ---- Relative soil moisture ---- 
fout <- paste(path,'FigSI_X_More_xSite_comps.pdf', sep='')
pdf(fout, width=9, height=7, compress = FALSE)
par(mfrow=c(5,1), mar=c(0,0,0,0), oma=c(5,5,3,2))

for (e in 1:length(case)) {
  CLM_date <- as.Date(as.character(DATE_mean[e,]),format="%Y%m%d")
  Ylim     <- c(0, 1)
  tempY1   <- BTRAN_mean[e,] 
  tempY2   <- FPG_mean[e,] 
  tempY1[GPP_mean[e,] <= 0] <- NA
  tempY2[GPP_mean[e,] <= 0] <- NA
  plot(tempY1 ~ CLM_date, type="l", 
       ylim=Ylim, lwd=3, col=4, 
       xlab=NA, ylab=NA, xaxt='n')
  lines(tempY2 ~ CLM_date, lwd=3, col=6)
  text(x=CLM_date[30],y=Ylim[2]*0.9,labels=good[e],cex=2)
  if (e == 3) {
    mtext("Limiting factor", side = 2, line = 3)
  }
  if (e == 5) {
    mtext("Year", side = 1, line = 3)
    axis(1, tck = 0.05, labels=year[2:7], at=at)
  }
  if (e == 1) {
    #legend('topright',c("water","nitrogen"), col=c(4,6), lty=1, lwd=3, bty="n")
  } 
}           # close case loop
dev.next()


for (e in 1:length(case)) {
  CLM_date <- as.Date(as.character(DATE_mean[e,]),format="%Y%m%d")
  Ylim     <- c(0, 2)
  plot(log(H2OSOI_mean[e,3,]/min(H2OSOI_mean[e,3,])) ~ CLM_date, type="l", 
       ylim=Ylim, lwd=3, col=2, 
       xlab=NA, ylab=NA, xaxt='n')
  text(x=CLM_date[30],y=Ylim[2]*0.9,labels=good[e],cex=2)
  if (e == 3) {
    mtext("Relative Soil Moisture", side = 2, line = 3)
  }
  if (e == 5) {
    mtext("Year", side = 1, line = 3)
    axis(1, tck = 0.05, labels=year[2:7], at=at)
  }
  if (e == 1) {
    yplot <- log(SOIMOI/min(SOIMOI, na.rm=T))
    lines(yplot ~ FLX_DATE, col=4, lwd=1) 
    legend('topright',c("CLM4.5","Obs."), col=c(2,4), lty=1, lwd=3, bty="n")
  } 
  if (e == 2) {
    yplot <- log(SOIMOI_DM/min(SOIMOI_DM, na.rm=T))
    lines(yplot ~ DATE_DM, col=4, lwd=1) 
  }  
  if (e == 3) {
    yplot <- log(SOIMOI_MM/min(SOIMOI_MM, na.rm=T))
    lines(yplot ~ DATE_MM, col=4, lwd=1) 
  }  
  if (e == 4) {
    yplot <- log(SOIMOI_WM/min(SOIMOI_WM, na.rm=T))
    lines(yplot ~ DATE_WM, col=4, lwd=1) 
  }  
}           # close case loop
dev.next()

#  ---- Soil Temp  ---- 
for (e in 1:length(case)) {
  CLM_date <- as.Date(as.character(DATE_mean[e,]),format="%Y%m%d")
  Ylim     <- c(min(TSOI_mean,na.rm=T), max(TSOI_mean,na.rm=T))
  
  plot(TSOI_mean[e,3,] ~ CLM_date, type="l", ylim=Ylim, lwd=3, col=2, 
       xlab=NA, ylab=NA, xaxt='n')
  abline(h=0,lty=2)
  text(x=CLM_date[30],y=Ylim[2]*0.9,labels=good[e],cex=2)
  if (e == 1) {
    legend('topright',c("CLM4.5","Saddle Obs.","TVan"), 
           col=c(2,1,4), lty=1, lwd=3, bty="n")
    lines(TSOI_FF ~ TEMP_DATE_FF, col=4, lwd=1) 
  } 
  if (e == 2) {
    lines(TSOI_DM ~ DATE_DM, col=1, lwd=1) 
    lines(TSOI_DM2 ~ TEMP_DATE_DM, col=4, lwd=1) 
  }  
  if (e == 3) {
    mtext("Soil Temperature (˚C)", side = 2, line = 3)
    lines(TSOI_MM ~ TEMP_DATE_MM, col=4, lwd=1) 
  }
  if (e == 4) {
    lines(TSOI_WM ~ TEMP_DATE_WM, col=4, lwd=1) 
  }
  if (e == 5) {
    mtext("Year", side = 1, line = 3)
    axis(1, tck = 0.05, labels=year[2:7], at=at)
  }
}           # close case loop
dev.next()


#  ---- GPP ---- 
for (e in 1:length(case)) {
  CLM_date <- as.Date(as.character(DATE_mean[e,]),format="%Y%m%d")
  Ylim     <- c(-2, max(GPP_mean))
  plot(GPP_mean[e,] ~ CLM_date, type="l", ylim=Ylim, lwd=3, col=2, 
       xlab=NA, ylab=NA, xaxt='n')
  text(x=CLM_date[30],y=Ylim[2]*0.9,labels=good[e],cex=2)
  if (e == 3) {
    mtext(expression(paste("GPP (gC",m^-2," ",d^-1,")",sep="")), 
          side = 2, line = 3)
  }
  if (e == 5) {
    mtext("Year", side = 1, line = 3)
    axis(1, tck = 0.05, labels=year[2:7], at=at)
  }
  if (e == 1) {
    lines(GPP_OBS ~ FLX_DATE, col=4, lwd=1, type='l') 
    legend('topright',c("CLM4.5","TVan Obs."), col=c(2,4), 
           lty=1, lwd=3, bty="n")
  }  
  abline(h=0, lty=2)
}           # close case loop
dev.next()

#  ---- NEE ---- 
for (e in 1:length(case)) {
  CLM_date <- as.Date(as.character(DATE_mean[e,]),format="%Y%m%d")
  Ylim     <- c(min(NEE_mean), max(NEE_mean)*1.2)
  plot(NEE_mean[e,] ~ CLM_date, type="l", ylim=Ylim, lwd=3, col=2, 
       xlab=NA, ylab=NA, xaxt='n')
  text(x=CLM_date[30],y=Ylim[2]*0.9,labels=good[e],cex=2)
  if (e == 3) {
    mtext(expression(paste("NEE (gC",m^-2," ",d^-1,")",sep="")), 
          side = 2, line = 3)
  }
  if (e == 5) {
    mtext("Year", side = 1, line = 3)
    axis(1, tck = 0.05, labels=year[2:7], at=at)
  }
  if (e == 1) {
    lines(NEE_OBS ~ FLX_DATE, col=4, lwd=1, type='l') 
    legend('bottomright',c("CLM4.5","TVan Obs."), col=c(2,4), 
           lty=1, lwd=3, bty="n")
  }  
  abline(h=0, lty=2)
}           # close case loop
dev.off()
print(paste('wrote',fout))


# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# mean annual cycle of limitation
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
month   <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","")
days    <- c( 31  ,  28 ,  31 ,  30 ,   31,   30,   31,   31,   30,   31,   30,  31 ) 
totdays <- rep(0, 13)
for (m in 2:13) { totdays[m] <- sum(days[1:m-1]) }

clmdate   <- as.Date(as.character(mcdate),format='%Y%m%d')
clmday    <- (as.character(clmdate, format = "%m-%d")) 
clmyear   <- as.numeric(as.character(clmdate, format = "%Y")) 
clmmonth  <- as.numeric(format.Date(clmdate, "%m"))
clmDATE   <- clmdate[mchour==1.0]
clmDAY    <- as.character(clmDATE, format = "%m-%d") 

# arrays to store results
dims       <- c(ncase, 366)
GPP_ann    <- array(NA, dims)
GPP_sd     <- array(NA, dims)
FPG_ann    <- array(NA, dims)
FPG_sd     <- array(NA, dims)
BTRAN_ann  <- array(NA, dims)
BTRAN_sd   <- array(NA, dims)

# calculate daily statistics over all years
for (e in 1:ncase) {
  GPP_ann[e,]   <- tapply(GPP_mean[e,], clmDAY, mean)
  GPP_sd[e,]    <- tapply(GPP_mean[e,], clmDAY, sd)
  FPG_ann[e,]   <- tapply(1-FPG_mean[e,], clmDAY, mean)
  FPG_sd[e,]    <- tapply(1-FPG_mean[e,], clmDAY, sd)
  BTRAN_ann[e,] <- tapply(1-BTRAN_mean[e,], clmDAY, mean)
  BTRAN_sd[e,]  <- tapply(1-BTRAN_mean[e,], clmDAY, sd)
}
  
# mask out where GPP = 0
GPP_rel       <- GPP_ann / max(GPP_ann) 
BTRAN_ann[GPP_ann == 0] <- NA
FPG_ann[GPP_ann == 0]   <- NA
BTRAN_sd[GPP_ann == 0]  <- NA
FPG_sd[GPP_ann == 0]    <- NA

fout <- paste(path,'FigSI_2_Seasonal_limitation.pdf', sep='')
pdf(fout, width=9, height=7, compress = FALSE)
par(mfrow=c(5,1), mar=c(0,0,0,0), oma=c(5,5,3,5))

for (e in 1:ncase) {
  # plot mean annual cycle
  Ylim=c(0,0.9)
  Xlim=c(110,310)
  plot(BTRAN_ann[e,], type='l', col=4, lwd=3,
       xlim=Xlim,ylim =Ylim, 
       xlab=NA, ylab=NA, 
       xaxt='n',yaxt='n')
  text(x=Xlim[1]+2,y=Ylim[2]*0.9,labels=good[e],cex=2)
  
  x        <- seq(1,length(GPP_ann[e,]),1)
  xx       <- c(x, rev(x))
  yy0      <- c((GPP_ann[e,]  +GPP_sd[e,])  ,rev(GPP_ann[e,]  -GPP_sd[e,])  )
  yy1      <- c((BTRAN_ann[e,]+BTRAN_sd[e,]),rev(BTRAN_ann[e,]-BTRAN_sd[e,]))
  yy2      <- c((FPG_ann[e,]  +FPG_sd[e,]),  rev(FPG_ann[e,]  -FPG_sd[e,])  )
  polygon(na.omit(cbind(x = xx, y = yy1)),col = rgb(0,  0, 1 ,0.3), border = NA)
  polygon(na.omit(cbind(x = xx, y = yy2)),col = rgb(0,158/255,115/255, 0.5),border = NA)
  lines(  FPG_ann[e,], col = "#009E73", lwd=3)
  lines(BTRAN_ann[e,], col=4, lwd=3)
  #lines(  GPP_rel[e,], col=2, lwd=3)
  at <- seq(0.0,1,0.5)
  axis(1, tck = 0.05, at=c(totdays[4:10]), labels=NA)  
  axis(2, at=at)
  if (e == 3) { 
    mtext('Limitation (unitless)', side = 2, line = 3) 
    mtext(expression(paste("GPP (gC ",m^-2," ",d^-1,")")), 
          side = 4, line = 3) 
  }
  if (e == 5) {
#    mtext("Day of Year", side = 1, line = 3)
    axis(1, tck = 0.05, at=c(totdays[4:10]), labels=month[4:10], cex.lab=1.2)
  }
  
  # add GPP to second y axis
  par(new = T)
  at=seq(0,9,3)
  plot(GPP_ann[e,], col=2, lwd=2, type='l', 
       axes=F, xlab=NA, ylab=NA, 
       xlim=c(110,310),ylim =c(0,8))
  axis(side = 4, at=at)
  if (e == 3) { 

  }
  polygon(na.omit(cbind(x = xx, y = yy0)),col = rgb(1, 0, 0,0.2), border = NA)
  lines(GPP_ann[e,], lwd=2, col=2)  
  abline(h=0, lty=2)
  
}

dev.off()



# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# plot FF data only for soil moisture, temperature, GPP & NEE
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
fout <- paste(path,'Tvan_',case[1],'annual_obs.pdf', sep='')
pdf(fout, width=9, height=7, compress = FALSE)


CLM_date <- as.Date(as.character(DATE_mean[e,]),format="%Y%m%d")
at       <- as.Date(c('2008-01-01','2009-01-01','2010-01-01',
                      '2011-01-01','2012-01-01','2013-01-01','2014-01-01'))
par(mfrow=c(4,1), mar=c(0,0,0.5,0), oma=c(3,5,1,2))

plot(TSOI_mean[1,3,] ~ CLM_date, type="l", lwd=3, col=2, 
     xlab=NA, ylab=NA, xaxt='n')
lines(TSOI_FF ~ TEMP_DATE_FF, col=1, lwd=3) 
lines(TSOI_mean[1,3,] ~ CLM_date, col=2, lwd=1.5)
abline(h=0,lty=2)
mtext("Soil Temp.", side = 2, line = 3.5)
mtext(parse(text = paste("(.",'^o', "*C)", sep="")), side = 2, line=2.2, cex=0.75)
axis(1, tck = 0.05, labels = FALSE, at=at)
#legend('topright',c("CLM4.5","TVan Obs"), 
#       col=c(2,1), lty=1, lwd=3, bty="n")

Ylim=c(0,1.2)
#Ylim=c(0,0.4)
yplot0 <- log(H2OSOI_mean[1,3,]/min(H2OSOI_mean[1,3,]))
#yplot0 <- H2OSOI_mean[1,2,]
#yplot0[TSOI_mean[1,3,] < -4] <- NA
plot(yplot0 ~ CLM_date, type="l", 
     ylim=Ylim, lwd=3, col=2, 
     xlab=NA, ylab=NA, xaxt='n')
yplot <- log(SOIMOI/min(SOIMOI, na.rm=T))
#yplot <- SOIMOI
#yplot[TSOI_FF< -4] <- NA
lines(yplot ~ FLX_DATE, col=1, lwd=3) 
lines(yplot0 ~ CLM_date,lwd=3, col=2)
mtext("Soil Moist.", side = 2, line = 2.2)
mtext("Norm.", side=2, line=3.5)
axis(1, tck = 0.05, labels = FALSE, at=at)

Ylim     <- c(-1.5, max(GPP_mean[1,]*1.05))
plot(GPP_mean[1,] ~ CLM_date, type="l", ylim=Ylim, lwd=3, col=2, 
     xlab=NA, ylab=NA, xaxt='n')
abline(h=0,lty=2)    
lines(GPP_OBS ~ FLX_DATE, col=1, lwd=3, type='l') 
lines(GPP_mean[1,]~CLM_date, col=2, lwd=2)
mtext("GPP", side = 2, line = 3.5)
mtext(expression(paste("(gC ",m^-2," ",d^-1,")")), side=2, line=1.8, cex=0.75)
axis(1, tck = 0.05, labels = FALSE, at=at)

Ylim     <- c(-4,4)
plot(NEE_mean[1,] ~ CLM_date, ylim=Ylim, type="l", lwd=3, col=2, 
     xlab=NA, ylab=NA, xaxt='n')
abline(h=0,lty=2)    
lines(NEE_OBS ~ FLX_DATE, col=1, lwd=3, type='l') 
lines(NEE_mean[1,]~CLM_date, col=2, lwd=3)
mtext("NEE", side = 2, line = 3.5)
mtext(expression(paste("(gC ",m^-2," ",d^-1,")")), side=2, line=1.8, cex=0.75)
axis(1, tck = 0.05, labels=seq(2008,2014,1), at=at, outer=F)

dev.off()
print(paste('wrote',fout))
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# plot annual cycle of GPP values
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
remove(year)
year   <- seq(2008,2013,1) 
nyear  <- length(year)
hour <- seq(0,23.5,0.5)
month  <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
days   <- c( 31  ,  28 ,  31 ,  30 ,   31,   30,   31,   31,   30,   31,   30,  31 ) 
DOY    <- c(15,40,60,75,90,115,150,215,315) # random days to plot below
MONTH   <- c('01','02','03','04','05','06','07','08','09','10','11','12')
nMONTH  <- length(MONTH)
ydays   <- c(366,365,365,365,366,365)
sum(ydays)
paste(year[1],'-',MONTH[1],'-01',sep='')

clmdate  <- as.Date(as.character(DATE_mean[1,]),format='%Y%m%d')
clmday      <- (as.character(clmdate, format = "%m-%d")) 
length(clmday)
clmyear     <- as.numeric(as.character(clmdate, format = "%Y")) 

ndays       <- 366
dims        <- c(nrows,ndays)
dims2       <- c(nrows,nyear,ndays)
GPP_ANN_mean<- array(NA,dim=dims)
GPP_i_mean  <- array(NA,dim=dims2)
GPP_ANN_sd  <- array(NA,dim=dims)
for(e in 1:nrows) {
  GPP_ANN_mean[e,] <- tapply(GPP_mean[e,], clmday, mean) 
  GPP_ANN_sd[e,]   <- tapply(GPP_mean[e,], clmday, sd)   
  di               <- 1
  for (y in 1:nyear) {
    df  <- di + ydays[y] - 1
    GPP_i_mean[e,y,1:ydays[y]] <- GPP_mean[e,di:df]
    di  <- df + 1
  }
}

s2y          <- 365 * 24 * 3600


Data.flx$day <- as.character(paste(Data.flx$MO,"-",Data.flx$DD, sep=""), format ="%m-%d")
Data.flx$day <- as.Date(Data.flx$day, format ="%m-%d")
OBS_DAY      <- tapply(as.numeric(Data.flx$day), Data.flx$day, mean, na.rm=T)
OBS_GPP_mean <- tapply(Data.flx$GPP_f, Data.flx$day, mean, na.rm=T)
OBS_GPP_sd   <- tapply(Data.flx$GPP_f, Data.flx$day, sd,   na.rm=T)
OBS_GPP_mean <- OBS_GPP_mean * 1e-3 * 3600 * 24 * 12/44    #gC/m2/d from mgCO2/m2/s
OBS_GPP_sd   <- OBS_GPP_sd   * 1e-3 * 3600 * 24 * 12/44   #gC/m2/d
max(OBS_GPP_mean)

obsJuly <- Data.flx$GPP_f[Data.flx$MO==7] 
mean(obsJuly, na.rm=T) * 1e-3 * 3600 * 24 * 12/44

fout <- paste(path,'FigX_Tvan_Annual_GPP.pdf', sep='')
pdf(fout, width=7, height=6, compress = FALSE)
par(mar=c(4,4,1,1), oma=c(0,0,0,0),mgp=c(2.5,1,0), cex=1.3)
plot(OBS_GPP_mean, type="l", ylim=c(-1,7),lwd=4, xlab="Day",
     #       xlim=c(100,300),
     ylab=expression(paste("GPP (gC ",m^-2," ",d^-1,")")))
#x  <- seq(1,length(OBS_DAY),1)
#xx <- c(x, rev(x))
#yy <- c((OBS_GPP_mean+OBS_GPP_sd),rev(OBS_GPP_mean-OBS_GPP_sd))
#polygon(xx,yy,col = rgb(0, 0, 0,0.5), border = NA)
abline(h=0, lty=2) 

#x  <- seq(1,length(GPP_mean),1)
#xx <- c(x, rev(x))
#yy <- c((GPP_mean+GPP_sd),rev(GPP_mean-GPP_sd))
#polygon(xx,yy,col = rgb(1, 0, 0,0.3), border = NA)
for(e in 1:nrows) {
  lines(GPP_ANN_mean[e,], col=e+1, lwd=4)   
#  for (y in 1:nyear) {
#    lines(GPP_i_mean[e,y,], col=e+1, lwd=1)
#  }
}

dim(GPP_ANN_mean)
lines(GPP_ANN_mean, col=2,lwd=3)

legend('topleft',legend=c("CLM", "TVan obs"), lwd=4, col=c(2,1),bty='n',cex=1.1)        

dev.off()

sum(GPP_mean)
sum(OBS_GPP_mean)

GPP_mean2    <- array(NA, dim=c(nyear,366))
GPP_sd2      <- array(NA, dim=c(nyear,366))
length(clmyear)

quit()
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# isolate annual plots 
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
for (e in 1:length(case)) {
  infile   <- paste(path,pre,case[e],suf,sep='')
  Data.clm <- open.ncdf(infile)   

  #print(Data.clm)
  MCDATE         = get.var.ncdf(Data.clm, "mcdate") # getting/naming the variable
	MCSEC          = get.var.ncdf(Data.clm, "mcsec") 
  TSOI           = get.var.ncdf(Data.clm, "TSOI")         #soil temperature (K)
  H2OSOI         = get.var.ncdf(Data.clm, "H2OSOI")	     	#Volumetric soil moisture
  SNOW   		     = get.var.ncdf(Data.clm, "SNOW")       
  RAIN     	     = get.var.ncdf(Data.clm, "RAIN")       
	SOILLIQ	       = get.var.ncdf(Data.clm, "SOILLIQ")      #snow liquid water 
  SNOWICE	       = get.var.ncdf(Data.clm, "SNOWICE")      #snow ice water 
  SNOW_DEPTH     = get.var.ncdf(Data.clm, "SNOW_DEPTH")	#snow height of snow covered area

  H2OSOIunits <- att.get.ncdf(Data.clm, "H2OSOI", "units")
  H2OSOIunits
  mean(H2OSOI)
  TSOI   <- TSOI - 273.15
	nsteps <- length(MCDATE)
	MCDATE[nsteps]
	MCSEC[2]
	MCDATE[nsteps-48]
  n1 <- 365*48  #for no-leap   data
  n2 <- n1 + 48 #for leap year data (leap)

  shift      <- 15  #adjust CLM to MST  
  mcdate     <- MCDATE[1:(nsteps-shift+1)]
  mcsec      <- MCSEC[1:(nsteps-shift+1)]
  mchour     <- mcsec/3600
  tsoi       <- TSOI[4,1+shift:nsteps-1]   #select 4th soil layer
  snow       <- SNOW[1+shift:nsteps-1]
  rain       <- RAIN[1+shift:nsteps-1]
  soilliq    <- SOILLIQ[1+shift:nsteps-1]
  snowice    <- SNOWICE[1+shift:nsteps-1]
  snow_depth <- SNOW_DEPTH[1+shift:nsteps-1]
  h2osoi     <- H2OSOI[,1+shift:nsteps-1]	
  nstep2     <- length(mcdate)
  precip     <- rain + snow
  dim(h2osoi)
  mcdate[1]
  mcdate[nstep2]
  day <- as.factor(mcdate)
    
  DATE_mean    <- tapply(mcdate, day, mean)
  TSOI_mean    <- tapply(tsoi,   day, mean)
  SNOW_sum     <- tapply(snow,   day, sum)
  RAIN_sum     <- tapply(rain,   day, sum)
  SOILLIQ_mean <- tapply(soilliq,day, mean)
  SNOWICE_mean <- tapply(snowice,day, mean)
  SNOW_DEPTH_mean <- tapply(snow_depth,day, mean)
  nlevsoi     <- 10
  ndaily      <- length(DATE_mean)
  H2OSOI_mean <- array(NA, c(nlevsoi,ndaily))	
  for (i in 1:nlevsoi) {
  	H2OSOI_mean[i,]  <- tapply(h2osoi[i,],day, mean)     	
  }
  PPT_sum      <- SNOW_sum + RAIN_sum
#  SWE_mean     <- SNOWLIQ_mean + SNOWICE_mean
  
  close.ncdf(Data.clm)  
  
#---------------------------------------------------------
#   Plot model vs. obs.
#---------------------------------------------------------
  year     <- seq(2007,2013,1) 
  nyear    <- length(year)
  CLM_date <- as.Date(as.character(DATE_mean),format="%Y%m%d")
  OBS_TIME <- as.Date(obstime)
  ylim     <- c(0, max(depth_mean[e,]+depth_sd[e,]))

  fout <- paste(path,'Figs/Tvan',case[e],'_SnowDepth.pdf', sep='')
  pdf(fout, width=9, height=7, compress = FALSE)

  par()
  xx1 <- c(OBS_TIME, rev(OBS_TIME))
  yy1 <- c((depth_mean[e,] + depth_sd[e,]), rev(depth_mean[e,] - depth_sd[e,]))
  Ylim=c(0,max(yy1))

  plot(depth_mean[e,] ~ OBS_TIME, type="l", ylim=Ylim, lwd=3, 
       xlab=NA, ylab='depth (cm)', main=good[e])
  polygon(xx1, yy1, col = rgb(0, 0, 0,0.35), border = NA)
  lines(SNOW_DEPTH_mean*100 ~ CLM_date, col=2, lwd=3, type='l') 
  points(depth_mean[e,] ~ OBS_TIME, cex=1.5,pch=16)
  if (e<=2) {
    lines(TVAN_depth ~ TVAN_DATE, col=4, lwd=3, type='l') 
    legend('topright',c("CLM4.5","Saddle Obs.","TVan Obs."), col=c(2,1,4), lty=1, lwd=3, bty="n")
  } else {  
    legend('topright',c("CLM4.5","Saddle Obs."), col=c(2,1), lty=1, lwd=3, bty="n")
  }
  dev.next()

  par(mfrow=c(3,2), mar=c(0,0,0,0), oma=c(5,5,3,2))

  for (y in 1:(nyear-1)) {
    #select first and last time points for each snow year
    f_clm <- year[y] * 1e4 + 0930
    l_clm <- (year[y]+1) * 1e4 + 1001
    f_obs <- as.Date(paste(year[y],'-09-30',sep=""))
    l_obs <- as.Date(paste(year[y]+1,'-10-01',sep=""))
    #subset data top plot  
    CLM_x  <- CLM_date[DATE_mean >= f_clm & DATE_mean <= l_clm]
    CLM_y  <- SNOW_DEPTH_mean[DATE_mean >= f_clm & DATE_mean <= l_clm]
    OBS_x  <- OBS_TIME[OBS_TIME >= f_obs & OBS_TIME <= l_obs] 
    OBS_y1 <- depth_mean[e,][OBS_TIME >= f_obs & OBS_TIME <= l_obs] #select veg community
    OBS_y2 <- depth_sd[e,][OBS_TIME >= f_obs & OBS_TIME <= l_obs] 
    TVAN_x <- TVAN_DATE[TVAN_DATE >= f_obs & TVAN_DATE <= l_obs]
    TVAN_y <- TVAN_depth[TVAN_DATE >= f_obs & TVAN_DATE <= l_obs]

    xx <- c(OBS_x, rev(OBS_x))
    yy <- c((OBS_y1 + OBS_y2), rev(OBS_y1 - OBS_y2))
    yy[is.na(yy)] <- 0
    
    plot(OBS_y1 ~ OBS_x, type="l", ylim=ylim, lwd=3, yaxt = "n", 
         xlab=NA, ylab=NA, xlim=c(f_obs,l_obs))
    polygon(xx, yy, col = rgb(0, 0, 0,0.35), border = NA)
    lines(CLM_y * 100 ~ CLM_x, col=2, lwd=3)  
    points(OBS_y1 ~ OBS_x, cex=1.5, pch=16)
    text(x=CLM_x[330],y=ylim[2]*0.1,labels=paste(year[y],'-',year[y]+1,sep=""))
    if (e <= 2) {
      lines(TVAN_y ~ TVAN_x, col=4, lwd=3, type='l') 
    }

    if (y == 1 || y == 3 || y == 5 ) { axis(2, tck = 0.05) } 
    else { axis(2, labels=NA, tck = 0.05)  }  
  
    #if (y >= 4) { axis(1, ,tck = 0.05) } 
    #else { axis(1, labels=mo, tck = 0.05) }
  
    if (y == 3 ) {mtext("Snow Depth (cm)", side = 2, line = 3)}
  
    if (y == 1) {
      if (e<=2) {
        legend('topright',c("CLM4.5","Saddle Obs.","TVan Obs."), col=c(2,1,4), lty=1, lwd=3, bty="n")
      } else {  
        legend('topright',c("CLM4.5","Saddle Grid"), col=c(2,1), lty=1, lwd=2, bty="n")
      }
    }
    
  }         # close annual loop
  print(paste('fished',case[e]))
  dev.off()

}           # close case loop



xx <- c(OBS_TIME, rev(OBS_TIME))
yy1 <- c((depth_mean[1,] + depth_sd[1,]), rev(depth_mean[1,] - depth_sd[1,]))
yy2 <- c((depth_mean[2,] + depth_sd[2,]), rev(depth_mean[2,] - depth_sd[2,]))

ylim=c(0,max(yy1))
plot(depth_mean[1,] ~ OBS_TIME, type="l", ylim=ylim, lwd=3, 
     xlab=NA, ylab=NA, main="'FF & DM obs")
polygon(xx, yy1, col = rgb(0, 0, 0,0.35), border = NA)
polygon(xx, yy2, col = rgb(0, 0, 1,0.35), border = NA)
lines(depth_mean[2,] ~ OBS_TIME, col=4, lwd=3)  


xx <- c(OBS_TIME, rev(OBS_TIME))
yy1 <- c((depth_mean[3,] + depth_sd[3,]), rev(depth_mean[3,] - depth_sd[3,]))
yy2 <- c((depth_mean[4,] + depth_sd[3,]), rev(depth_mean[4,] - depth_sd[4,]))
ylim=c(0,max(yy1))
plot(depth_mean[3,] ~ OBS_TIME, type="l", ylim=ylim, lwd=3, 
     xlab=NA, ylab=NA, main="'MM & WM obs")
polygon(xx, yy1, col = rgb(0, 0, 0,0.35), border = NA)
polygon(xx, yy2, col = rgb(0, 0, 1,0.35), border = NA)
lines(depth_mean[4,] ~ OBS_TIME, col=4, lwd=3)  


fday  <- 365 * 1
lday  <- 365 * 2
nbot  <- 11
e <- 5
x <- CLM_date[fday:lday]
y <- -rev(levgrnd[1:nbot])
z <- H2OSOI_mean[e,1:nbot,fday:lday]
z <- t(z[nrow(z):1,]) 
filled.contour(x,y, z, color = function(z)rev(rainbow(z)), nlevels=20,
               main = paste(case[e], 'H2OSOI (C)'), zlim=range(H2OSOI_mean, finite = TRUE)*0.8)

dim(TSOI_mean[3,,])
dim(TSOI_mean)
plot()
#---------------------------------------------------------
#----------------  END  ----------------------------------
#---------------------------------------------------------
  
  
