# Compares vertically resolved temperature and moisture profiles

remove(list=ls())
library(REddyProc)
library(ncdf)
library(lattice)

#-----------------------------------------------------------------------
#---------------read in CLM variables----------------------------------
#-----------------------------------------------------------------------
dir      <- ('/Users/wwieder/Desktop/Working_files/Niwot/NR_fluxes/TVan/Tvan_PTCLM45/')
path     <- paste(dir,'CLM_nc_files/PTrespmods/lowWATSAT/',sep='')
pre      <- ('Tvan_PTrespmods_allPHYS_noLW_lowVCmax_')
suf      <- ('.clm2.h1.2008-01-01-00000.nc')
exp      <- c('075_wm_noRHpsn_lowRESIST_run',
              '075_wm_noRHpsn_lowRESIST_M-Awarm',
              '075_wm_noRHpsn_lowRESIST_BlackSand')


nrows    <- length(exp)
nlevgrnd <- 25
nyears   <- 6
ndays    <- 366 # to calculate DOY average
#ndays    <- 365 * nyears + 2 # for look at interannual variability
dims1           <- c(nrows,ndays)
dims            <- c(nrows,nlevgrnd,ndays)
DATE_mean       <- array(NA,dim=dims1)
LH_mean         <- array(NA,dim=dims1)
FSH_mean        <- array(NA,dim=dims1)
GPP_mean        <- array(NA,dim=dims1)
SNOW_mean       <- array(NA,dim=dims1)
SNOW_DEPTH_mean <- array(NA,dim=dims1)
QVEGT_mean      <- array(NA,dim=dims1)
QVEGE_mean      <- array(NA,dim=dims1)
QSOIL_mean      <- array(NA,dim=dims1)
QRUNOFF_mean    <- array(NA,dim=dims1)
TV_mean         <- array(NA,dim=dims1)
TV_min          <- array(NA,dim=dims1)
TV_max          <- array(NA,dim=dims1)
TG_min          <- array(NA,dim=dims1)
TSOI_mean       <- array(NA,dim=dims)
SOILPSI_mean    <- array(NA,dim=dims)
SOILLIQ_mean    <- array(NA,dim=dims)
H2OSOI_mean     <- array(NA,dim=dims)       

TSOI_sd       <- array(NA,dim=dims)
SOILPSI_sd    <- array(NA,dim=dims)
SOILLIQ_sd    <- array(NA,dim=dims)
H2OSOI_sd     <- array(NA,dim=dims)       

for (e in 1:length(exp)) {  
  infile   <- paste(path,pre,exp[e],suf,sep='')
  Data.clm <- open.ncdf(infile)   
  print(paste('read',infile))  
  #print(Data.clm)
  nbedrock    <- get.var.ncdf(Data.clm,'nbedrock')
  levgrnd     <- get.var.ncdf(Data.clm,'levgrnd')
  MCDATE      <- get.var.ncdf(Data.clm, "mcdate") # getting/naming the variable
  MCSEC       <- get.var.ncdf(Data.clm, "mcsec") 
  QVEGT       <- get.var.ncdf(Data.clm, "QVEGT") 
  QVEGE       <- get.var.ncdf(Data.clm, "QVEGE") 
  QSOIL       <- get.var.ncdf(Data.clm, "QSOIL") 
  QRUNOFF     <- get.var.ncdf(Data.clm, "QRUNOFF") 
  GPP         <- get.var.ncdf(Data.clm, "GPP")          
  RAIN        <- get.var.ncdf(Data.clm, "RAIN")         
  SNOW        <- get.var.ncdf(Data.clm, "SNOW")         
  TSOI        <- get.var.ncdf(Data.clm, "TSOI")         #soil temperature (K)
  TG          <- get.var.ncdf(Data.clm, "TG")           #ground temperature (K)
  TV          <- get.var.ncdf(Data.clm, "TV")           #veg temperature (K)
  H2OSOI      <- get.var.ncdf(Data.clm, "H2OSOI")       #Volumetric soil moisture
  SOILLIQ     <- get.var.ncdf(Data.clm, "SOILLIQ")      #soil liquid water 
  SNOW        <- get.var.ncdf(Data.clm, "SNOW") 
  SNOW_DEPTH  <- get.var.ncdf(Data.clm, "SNOW_DEPTH")   #snow height of snow covered area
  SOILLIQunits<- att.get.ncdf(Data.clm, "SOILLIQ", "units")
  H2OSOIunits <- att.get.ncdf(Data.clm, "H2OSOI", "units")
  TLAI        <- get.var.ncdf(Data.clm, "TLAI")         
  SOILPSI     <- get.var.ncdf(Data.clm, "SOILPSI")      #
  FSH            = get.var.ncdf(Data.clm, "FSH")        #sensible heat
  LH             = get.var.ncdf(Data.clm, "EFLX_LH_TOT")#latent heat

#  print(H2OSOIunits)
  TSOI   <- TSOI - 273.15
  TV     <- TV   - 273.15
  TG     <- TG   - 273.15
  nsteps <- length(MCDATE)
  MCDATE[nsteps-48]
  n1 <- 365*48  #for no-leap   data
  n2 <- n1 + 48 #for leap year data (leap)
  
  shift      <- 15  #adjust CLM to MST  
  mcdate     <- MCDATE[1:(nsteps-shift+1)]
  mcsec      <- MCSEC[1:(nsteps-shift+1)]
  mchour     <- mcsec/3600
  gpp        <- GPP[1+shift:nsteps-1]   
  lh         <- LH[1+shift:nsteps-1]  
  fsh        <- FSH[1+shift:nsteps-1]  
  #lh[gpp<=0] <- NA
  #fsh[gpp<=0]<- NA
  tg         <- TG[1+shift:nsteps-1]  
  tv         <- TV[1+shift:nsteps-1]   
  tsoi       <- TSOI[   ,1+shift:nsteps-1]   #select all soil layers
  soilliq    <- SOILLIQ[,1+shift:nsteps-1]
  h2osoi     <- H2OSOI[ ,1+shift:nsteps-1]	
  snow       <- SNOW[    1+shift:nsteps-1]
  qvegt      <- QVEGT[1+shift:nsteps-1]
  qvege      <- QVEGE[1+shift:nsteps-1]
  qsoil      <- QSOIL[1+shift:nsteps-1]
  qrunoff    <- QRUNOFF[1+shift:nsteps-1]
  snow_depth <- SNOW_DEPTH[1+shift:nsteps-1]
  soilpsi    <- SOILPSI[,1+shift:nsteps-1]
  nstep2     <- length(mcdate)
  clmdate  <- as.Date(as.character(mcdate),format='%Y%m%d')
  clmday      <- (as.character(clmdate, format = "%m-%d")) 
  
  day        <- as.factor(clmday) #to calculate DOY average
#  day        <- as.factor(mcdate) #to maintain inter annual variability
  
#  DATE_mean[e,]    <- tapply(mcdate, day, mean)
  GPP_mean[e,]     <- tapply(gpp,    day, mean)
  LH_mean[e,]      <- tapply(lh,     day, sum, na.rm=T)
  FSH_mean[e,]     <- tapply(fsh,    day, sum, na.rm=T)
  TG_min[e,]       <- tapply(tg,     day, min)
  TV_max[e,]       <- tapply(tv,     day, max)
  TV_min[e,]       <- tapply(tv,     day, min)
  TV_mean[e,]      <- tapply(tv,     day, mean)
  SNOW_mean[e,]    <- tapply(snow,   day, mean)
  QVEGT_mean[e,]   <- tapply(qvegt,  day, mean)
  QVEGE_mean[e,]   <- tapply(qvege,  day, mean)
  QSOIL_mean[e,]   <- tapply(qsoil,  day, mean)
  QRUNOFF_mean[e,] <- tapply(qrunoff,day, mean)
  SNOW_DEPTH_mean[e,] <- tapply(snow_depth,day, mean)
  for (z in 1:nlevgrnd) {
    TSOI_mean[e,z,]    <- tapply(tsoi[z,],    day, mean)
    SOILLIQ_mean[e,z,] <- tapply(soilliq[z,],day, mean)
    H2OSOI_mean[e,z,]  <- tapply(h2osoi[z,], day, mean)           
    SOILPSI_mean[e,z,] <- tapply(soilpsi[z,],day, mean)

    TSOI_sd[e,z,]    <- tapply(tsoi[z,],    day, sd)
    SOILLIQ_sd[e,z,] <- tapply(soilliq[z,],day, sd)
    H2OSOI_sd[e,z,]  <- tapply(h2osoi[z,], day, sd)           
    SOILPSI_sd[e,z,] <- tapply(soilpsi[z,],day, sd)
    
  }
 RUNOFFunits <- att.get.ncdf(Data.clm, "QRUNOFF", "units")
 QSOILunits  <- att.get.ncdf(Data.clm, "QSOIL",   "units")
 
 close.ncdf(Data.clm)  
}
FSH_mean[FSH_mean=='NaN'] <- NA
LH_mean[LH_mean=='NaN'] <- NA
clm_B <- FSH_mean / LH_mean
clmEF <- 1/(1+clm_B)



# ---------- make plots---------------
month   <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","")
days    <- c( 31  ,  28 ,  31 ,  30 ,   31,   30,   31,   31,   30,   31,   30,  31 ) 
totdays <- rep(0, 13)
for (m in 2:13) { totdays[m] <- sum(days[1:m-1]) }

fout <- paste(path,'FigSI_1_Tvan_WM_extendedSummer.pdf', sep='')
pdf(fout, width=7, height=6, compress = FALSE)
par(mfrow=c(3,1),mar=c(0,0,0,0), oma=c(5,6,3,2))

plot(GPP_mean[1,]*3600*24, type='l',xlim=c(0,360),ylim=range(GPP_mean*3600*24),
     lwd=3, xlab=NA, ylab=NA,xaxt='n', col=4) 
abline(h=0, lty=2) 
mtext('GPP', side = 2, line = 4)
mtext(expression(paste("(gC ",m^-2," ",d^-1,")")), 
      side = 2, line = 2, cex=0.8) 
text(0,5,'(a)', cex=1.5, adj=0)
lines(GPP_mean[2,]*3600*24, col=2, lwd=3)
lines(GPP_mean[3,]*3600*24, col=1, lwd=3)
axis(1, tck = 0.05, at=totdays, labels=NA)

plot(TSOI_mean[1,3,], type='l',xlim=c(0,360),ylim=range(TSOI_mean[,3,]),
     xaxt='n', lwd=3, ylab=NA, xlab=NA,xaxt='n', col=4) 
text(0,15,'(b)', cex=1.5, adj=0)
mtext('Soil Temp', side = 2, line = 4)
mtext(expression(paste('(',~degree~C,', 10 cm)'), sep=''), 
      side = 2, line = 2, cex=0.8) 
abline(h=0, lty=2) 
lines(TSOI_mean[2,3,], col=2, lwd=3)
lines(TSOI_mean[3,3,], col=1, lwd=3)
axis(1, tck = 0.05, at=totdays, labels=NA)

plot(SOILLIQ_mean[1,3,], type='l',xlim=c(0,360),ylim=range(SOILLIQ_mean[,3,]),
     xaxt='n', lwd=3, ylab=NA, xlab=NA,xaxt='n', col=4) 
text(0,21,'(c)', cex=1.5, adj=0)
abline(h=0, lty=2) 
lines(SOILLIQ_mean[2,3,], col=2, lwd=3)
lines(SOILLIQ_mean[3,3,], col=1, lwd=3)
mtext('Soil Moisture', side = 2, line = 4)
mtext(expression(paste("(kg ",m^-2,")")), 
      side = 2, line = 2, cex=0.8) 
axis(1, tck = 0.05, at=totdays, labels=month, cex=1.2)
dev.off()

par(mfrow=c(1,1),mar=c(0,0,0,0), oma=c(5,6,3,2))
TV_min[GPP_mean<=0] <- NA
TV_max[GPP_mean<=0] <- NA
plot(TV_min[1,], type='l',xlim=c(145,180),ylim=range(TV_min, na.rm=T),
     xaxt='n', lwd=3, ylab=NA, xlab=NA, col=4) 
text(0,21,'(c)', cex=1.5, adj=0)
abline(h=0, lty=2) 
lines(TV_min[2,], col=2, lwd=3)
lines(TV_min[3,], col=1, lwd=3)
lines(TV_min[1,], col=4, lwd=3)
mtext('Min TV', side = 2, line = 4)
mtext(expression(paste('(',~degree~C,')'), sep=''), 
      side = 2, line = 2, cex=0.8) 
axis(1, tck = 0.05,  cex=1.2)

TV_range <- TV_max-TV_min
plot(TV_range[1,], type='l',xlim=c(145,180),ylim=range(TV_range, na.rm=T),
     xaxt='n', lwd=3, ylab=NA, xlab=NA, col=4) 
text(0,21,'(c)', cex=1.5, adj=0)
abline(h=0, lty=2) 
lines(TV_range[2,], col=2, lwd=3)
lines(TV_range[3,], col=1, lwd=3)
lines(TV_range[1,], col=4, lwd=3)
mtext('range TV', side = 2, line = 4)
mtext(expression(paste('(',~degree~C,')'), sep=''), 
      side = 2, line = 2, cex=0.8) 
axis(1, tck = 0.05,  cex=1.2)



plot(clmEF[1,], type='l',xlim=c(0,360),ylim=c(-0.1,1.1),
     lwd=3, xlab=NA, ylab=NA,xaxt='n', col=4) 
mtext('Evap. Frac.', side = 2, line = 4)
text(0,0.8,'(d)', cex=1.5, adj=0)
lines(clmEF[2,], col=2, lwd=3)
lines(clmEF[3,], col=1, lwd=3)
axis(1, tck = 0.05, at=totdays, labels=NA)






#---------Other plots to see--------------
QET_mean <- QSOIL_mean + QVEGE_mean + QVEGT_mean

plot( QET_mean[1,], type='l', col=2)
lines(QET_mean[2,], col=4)

plot( QSOIL_mean[2,], type='l', col=4)
lines(QSOIL_mean[1,], type='l', col=2)

rowMeans(GPP_mean)     * 3600 * 24 * 365
rowMeans(QET_mean)     * 3600 * 24 * 365
rowMeans(QVEGT_mean)   * 3600 * 24 * 365
rowMeans(QSOIL_mean)   * 3600 * 24 * 365
rowMeans(QRUNOFF_mean) * 3600 * 24 * 365

temp <- 48*365*1
plot(TV_mean[1,(365:700)], type="l")
lines(TV_mean[2,(365:700)], col=2)

plot(SOILPSI[1,1:temp], type='l', ylim=c(-5,0))
lines(SOILPSI[2,1:temp], col=2)
lines(SOILPSI[3,1:temp], col=3)
lines(SOILPSI[4,1:temp], col=4)

TV_diff      <- TV_mean[2,]       -TV_mean[1,]
QRUNOFF_diff <- QRUNOFF_mean[2,]  - QRUNOFF_mean[1,]
QSOIL_diff   <- QSOIL_mean[2,]    - QSOIL_mean[1,]
TSOI_diff    <- TSOI_mean[2,,]    - TSOI_mean[1,,]
SOILPSI_diff <- SOILPSI_mean[2,,] - SOILPSI_mean[1,,]
SOILLIQ_diff <- SOILLIQ_mean[2,,] - SOILLIQ_mean[1,,]
plot(QSOIL_diff, type='l')
plot(QRUNOFF_diff, type='l')
abline(h=0)
plot(SOILPSI_diff[4,], type='l', col=4)
lines(SOILPSI_diff[3,], col=3)
lines(SOILPSI_diff[2,], col=2)

H2OSOI_diff  <- H2OSOI_mean[2,,]  - H2OSOI_mean[1,,]


# calculate annual cycle of soil temp and moisture
clmdate  <- as.Date(as.character(mcdate),format='%Y%m%d')
clmday      <- (as.character(clmdate, format = "%m-%d")) 
clmyear     <- as.numeric(as.character(clmdate, format = "%Y")) 
years       <- seq(2008,2013,1)

temp <- TSOI_mean
tempOBS <- Data.flx$GPP_f[1:length(temp)]
temp[tempOBS==NA] <- NA
GPP_mean <- tapply(temp, clmday, mean, na.rm=T) 




for (e in 1:length(exp)) {
  CLM_date <- as.Date(as.character(DATE_mean[2,]),format="%Y%m%d")  
  fout <- paste(path[e],'VRsoils',exp[e],'_',case[1],'.pdf', sep='')
  pdf(fout, width=9, height=3, compress = FALSE)
  dim(TSOI_mean)
  par(mar=c(2,4,2,2))
  fday  <- 365 * 4
  lday  <- 365 * 6
  nbot  <- 12
  x <- CLM_date[fday:lday]
  y <- -rev(levgrnd[1:nbot])
  z <- TSOI_mean[e,1:nbot,fday:lday]
  z <- t(z[nrow(z):1,]) 
  filled.contour(x,y, z, color = function(z)rev(rainbow(z)), nlevels=20,
                 main = paste(exp[e], 'TSOI (C)'), zlim=range(TSOI_mean, finite = TRUE)*0.8)
  dev.next()

  z <- TSOI_diff[1:nbot,fday:lday]
  z <- t(z[nrow(z):1,]) 
  filled.contour(x,y, z, color = function(z)rev(rainbow(z)), nlevels=10,
                 main = paste(e, 'TSOI diff (2-1)'), zlim=range(z, finite = TRUE)*0.8)
  dev.next()

  z <- SOILPSI_mean[e,1:nbot,fday:lday]
  z <- t(z[nrow(z):1,]) 
  filled.contour(x,y,z , color = function(z)rev(topo.colors(z)), nlevels=20, 
                 main = paste(path[e], 'SOIPSI'), zlim=range(SOILPSI_mean, finite = TRUE))
  dev.next()
  
  z <- SOILPSI_diff[1:nbot,fday:lday]
  z <- t(z[nrow(z):1,]) 
  filled.contour(x,y, z, color = function(z)(heat.colors(z)), nlevels=10,
                 main = paste(e, 'SOILPSI diff (2-1)'), zlim=range(z, finite = TRUE)/2)
  dev.next()
e<- 1
  z <-  H2OSOI_mean[e,1:nbot,fday:lday]
  z <- t(z[nrow(z):1,]) 
  filled.contour(x,y,z , color=function(z)rev(topo.colors(z)), nlevels=20, 
                 main = paste(path[e],'H2OSOI'), zlim=range(H2OSOI_mean, finite = TRUE)*0.8)
  dev.next()
plot(H2OSOI_diff[1,fday:lday], type='l')
lines(H2OSOI_diff[2,fday:lday], col=2)
lines(H2OSOI_diff[3,fday:lday], col=3)
abline(h=0)
  z <- H2OSOI_diff[1:nbot,fday:lday]
  z <- t(z[nrow(z):1,]) 
  filled.contour(x,y, z, color = function(z)(heat.colors(z)), nlevels=10,
                 main = paste(e, 'H2OSOI diff (2-1)'), zlim=range(z, finite = TRUE)*0.9)
  dev.next()
  

  plot(GPP_mean[1,]~CLM_date, type='l', col=2)
  lines(GPP_mean[2,]~CLM_date, col=1)
  dev.next()
  
  plot(SNOW_mean[1,fday:lday]~CLM_date[fday:lday], type='l', col=2)
  lines(SNOW_mean[2,fday:lday]~CLM_date[fday:lday], col=1)
  dev.next()
  
  plot(SNOW_DEPTH_mean[1,]~CLM_date, type='l', col=2)
  lines(SNOW_DEPTH_mean[2,]~CLM_date, col=1)
  dev.next()
  
  plot(TSOI_mean[1,2,]~CLM_date, type='l', col=2)
  lines(TSOI_mean[2,2,]~CLM_date, col=1)
  dev.next()
  
  plot(SOILPSI_mean[1,2,]~CLM_date, type='l', col=2)
  lines(SOILPSI_mean[2,2,]~CLM_date, col=1)
  
  dev.off()
  #, color= color.bar, asp=0.98,  xlim=c(min(clay), max(clay)), ylim=c(0.06, max(met)), nlevels = 20,
  #               plot.title= title(xlab="clay fraction",ylab="metabolic fracion", cex.lab=1.4, cex.axis=1.4))
  print(paste('wrote',fout))
  
}
plot(SOILLIQ_mean[1,3,]~CLM_date, type='l',col=2)
plot(SOILLIQ_mean[2,3,]~CLM_date, type='l',col=1)

# ----------critical growing degree days ------------------
# ----------critical soil water potential  ------------------
plot(TV[1:(48*365)], type='l', ylim=c(260,295))
abline(v=180*48)
par(mfrow=c(4,1), mar=c(0,4,0,1))
plot(TSOI[3,], type='l', ylim=c(272,278), xaxt='n')
lines(TSOI[2,1:(48*365)], col=2)
lines(TSOI[3,1:(48*365)], col=3)
lines(TSOI[4,1:(48*365)], col=4)
abline(v=180*48)
abline(h = 276.5)
abline(h = 273.15)
plot(GPP[1:(48*365)], type='l',xaxt='n')
abline(v=180*48)
plot(TLAI[1:(48*365)], type='l',xaxt='n')
abline(v=180*48)
plot(SNOW_DEPTH[1:(48*365)], type='l')
abline(v=180*48)

par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(SOILPSI[3,1:(48*365)], type='l')
abline(h=soilpsi_off)

fracday         <- 1./48.
annavg_t2m      <- 271.957553457408
SHR_CONST_TKFRZ <- 273.15   # freezing point of water
crit_onset_gdd  <- exp(4.8 + 0.13 *(annavg_t2m - SHR_CONST_TKFRZ))
crit_onset_gdd2 <- exp(4.8 + 0.5 *(annavg_t2m - SHR_CONST_TKFRZ))
crit_onset_fdd  <- 15
crit_offset_swi <- 15 
crit_onset_swi  <- 15
soilpsi_off     <- -2
soilpsi_on      <- -2
dims            <- dim(TSOI)
Onset_gdd       <- array(0,dims) 
Onset_fdd       <- array(0,dims) 
Onset_gdd       <- array(0,dims) 
Onset_flag      <- array(0,dims) 
Dormant_flag    <- array(0,dims) 
Onset_gddflag   <- array(0,dims) 
Onset_swi       <- array(0,dims) 
dormant_flag    <- rep(1, dims[1])
onset_flag      <- rep(0, dims[1])
onset_gddflag   <- rep(0, dims[1])
onset_gdd       <- rep(0, dims[1])
onset_fdd       <- rep(0, dims[1])
onset_swi       <- rep(0, dims[1])
days_active     <- rep(0, dims[1])
for (t in 2:(48*365)) {
  for (l in 1:dims[1]){
    if (dormant_flag[l] == 1) {
      if (onset_gddflag[l] == 0 && TSOI[l,t] < SHR_CONST_TKFRZ) {
        onset_fdd[l] = onset_fdd[l] + fracday
        onset_gdd[l] = 0.
      } 
      
      # if the number of freezing degree days exceeds a critical value,
      # then onset will require both wet soils and a critical soil
      # temperature sum.  If this case is triggered, reset any previously
      # accumulated value in onset_swi, so that onset now depends on
      # the accumulated soil water index following the freeze trigger
      if (onset_fdd[l] > crit_onset_fdd) {
        onset_gddflag[l] <- 1.
        onset_fdd[l]     <- 0.
        onset_swi[l]     <- 0.
      }
      
      if (onset_gddflag[l] == 1 && TSOI[l,t] > SHR_CONST_TKFRZ) {
        onset_gdd[l] = onset_gdd[l] + (TSOI[l,t]-SHR_CONST_TKFRZ)*fracday
      }
      
      if (SOILPSI[l,t] >= soilpsi_on) {
        onset_swi[l] = onset_swi[l] + fracday
      }
      
      # only check soil temperature criteria if freeze flag set since
      # beginning of last dormancy.  If freeze flag set and growing
      # degree day sum (since freeze trigger) is lower than critical
      # value, then override the onset_flag set from soil water.
      
      if (onset_swi[l] > crit_onset_swi) {
        onset_flag[l] = 1
        if (onset_gddflag[l] == 1 && onset_gdd[l] < crit_onset_gdd) {
          onset_flag[l] = 0
        }
      }
      
      if (onset_flag[l] == 1) {
        dormant_flag[l] = 0
        days_active[l] = 0
        onset_gddflag[l] = 0
        onset_fdd[l] = 0
        onset_gdd[l] = 0
        onset_swi[l] = 0
        #        onset_counter[l] = ndays_on * secspday
      }  
      
      
    } else {
      
    }
    
    Onset_swi[l,t] <- onset_swi[l]
    Onset_gdd[l,t] <- onset_gdd[l]
    Onset_fdd[l,t] <- onset_fdd[l]
    Onset_gddflag[l,t] <- onset_gddflag[l]
    Onset_flag[l,t]    <- onset_flag[l]    
    Dormant_flag[l,t]  <- dormant_flag[l]    
  }  #close layer loop
} # close time loop


plot(Onset_fdd[3,1:(48*365)], type='l')
plot(Onset_gddflag[3,1:(48*365)], type='l')
plot(Onset_flag[3,1:(48*365)], type='l')
abline(v=180*48)

plot(Onset_gdd[1,1:(48*365)], type='l', ylim=c(-1,crit_onset_gdd*1.05))
lines(Onset_gdd[2,1:(48*365)], col=2)
lines(Onset_gdd[3,1:(48*365)], col=3)
lines(Onset_gdd[4,1:(48*365)], col=4)
abline (h=crit_onset_gdd)
abline(v=180*48)
abline (h=crit_onset_gdd2)

plot(Onset_swi[1,1:(48*365)], type='l', ylim=c(0,crit_onset_swi*1.05))
lines(Onset_swi[2,1:(48*365)], col=2)
lines(Onset_swi[3,1:(48*365)], col=3)
lines(Onset_swi[4,1:(48*365)], col=4)
abline (h=crit_onset_swi)

lines(Onset_swi[2,1:(48*365)], col=1, lwd=2)
lines(Onset_swi[2,1:(48*365)], col=2, lwd=2)
lines(Onset_swi[3,1:(48*365)], col=3, lwd=2)
lines(Onset_swi[4,1:(48*365)], col=4, lwd=2)