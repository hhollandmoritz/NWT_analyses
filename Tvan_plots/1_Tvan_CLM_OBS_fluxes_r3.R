# Compares observed and simulated fluxes, soil temperature 

# Uses ReddyProc to gap fill Eddy-Flux data
# https://www.bgc-jena.mpg.de/bgi/index.php/Services/REddyProcWeb
# install.packages("raster", repos="http://R-Forge.R-project.org", type="source")

remove(list=ls())

library(REddyProc)
library(ncdf)
library(lattice)

case <- c('010_run','010_noRHpsn_run',
          '010dry_noRHpsn_run','010dry_noRHpsn_LOWbetarun',
          '010dry_LOWbetarun','010dry_dm_noRHpsn_run', 
          '010dry_dm_noRHpsn_LOWbetarun', '015dry_noRHpsn_LOWbetarun',
          '010dry_run','010dry_noRHpsn_LOWbetarun_FF',
          '010dry_noRHpsn_LOWbeta_Sand_run','010dry_noRHpsn_LOWbeta_Sand_FF_run',
          '010dry_noRHpsn_LOWbeta_Sand_lowORG_run','010dry_noRHpsn_LOWbeta_Sand_lowORG_lowClay_run',
          '010dry_noRHpsn_LOWbeta_lowWATSAT_run','010dry_noRHpsn_LOWbeta_lowWATSAT_run1',
          '010dry_noRHpsn_LOWbeta_lowWATSAT_lowORG_run','010dry_noRHpsn_lowWATSAT_run',
          '010dry_noRHpsn_lowWATSAT_newTEXTURE_run','010dry_noRHpsn_lowWATSAT_upFSAT_run',
          '010dry_noRHpsn_lowWATSAT_upRESIST_run','010dry_noRHpsn_lowWATSAT_dsl05_run',
          '010dry_noRHpsn_lowWATSAT_dsl10_run')

# 1 = 70 cm soils, increased foliar C:N (32), increased froot:leaf (2.2), 
#     modified gdd calculation (CNPhenologyMod),
# 2 = remove RH control on gs (stomata conductance) in PhotosynthesisMod
# 3 = increase MAM snow to 25% of observed (improves net ratiation)
# 4 = Shallower rooting depth, 80% of fine roots in top 10 cm BEST ONE
# 5 = increased MAM snow and shallow roots, but w/ orig gs calculation
# 6 = as #3, but with deeper (100 cm) soils
# 7 = as #4, but with deeper (100 cm) soils 
# 8 = 15% snow, 30% in MAM- otherwise same as #4
# 9 = 35% bare ground, added back in RH control on gs, CHECK PARAMETER FILES
#10 = as #4, but w/ 35% bare ground  GPP goes way up, Btran decreases (likely averaging 0 from bare soil?)
#11 = as #4, but with 2x sand, to better approximate rocky soils
#12 = as 11, but with 35% base ground (to try and bring up gpp...)
#13 = as 11, but with half ORGANIC
#14 = as 13, with low clay content, to hopefully further dry out soils
#15 = as #4, but with 66% lower watsat (water holding capacity at porosity) to slimulate rocky soils?
#16 = as #4, but reduced watsat 50%
#17 = as 15, also reduced organic
#18 = as 15, but w/ default rooting depth
#19 = as 18, but w/ coarser soil texture to reduce summer soil moisture bias.
#20 = as 18, but w/ increases runnoff (from fsat)
#21 = as 18, but w/ increases soil evaporation (from modifying resistance *0.5)
#22 = as 18, but modified dry soil layer constant to 5 (from 15)
#23 = as 18, but modified dry soil layer constant to 10(from 15)

e   <- 21

#---------------read in CLM variables----------------------------------
dir      <- ('/Users/wwieder/Desktop/Working_files/Niwot/NR_fluxes/TVan/Tvan_PTCLM45/')
path     <- paste(dir,'CLM_nc_files/PTrespmods/lowWATSAT/',sep='')
pre      <- ('Tvan_PTrespmods_allPHYS_noLW_lowVCmax_')
suf      <- ('.clm2.h1.2008-01-01-00000.nc')

infile   <- paste(path,pre,case[e],suf,sep='')
Data.clm <- open.ncdf(infile)   
nsteps   <- 48 * (365*6 + 2)
  print(paste("The file has",Data.clm$nvars,"variables"))
  summary(Data.clm)
  print(Data.clm)
  MCDATE         = get.var.ncdf(Data.clm, "mcdate")[1:nsteps] # getting/naming the variable
	MCSEC          = get.var.ncdf(Data.clm, "mcsec")[1:nsteps] 
  TV             = get.var.ncdf(Data.clm, "TV")[1:nsteps]           #veg temperature
  T10            = get.var.ncdf(Data.clm, "T10")[1:nsteps]          #mean 2m air temp
  RH             = get.var.ncdf(Data.clm, "RH")[1:nsteps]           #relative humidity
  RH_LEAF        = get.var.ncdf(Data.clm, "RH_LEAF")[1:nsteps]      #leaf relative humidity
  BTRAN          = get.var.ncdf(Data.clm, "BTRAN")[1:nsteps]       
  FSH            = get.var.ncdf(Data.clm, "FSH")[1:nsteps]          #sensible heat
  LH             = get.var.ncdf(Data.clm, "EFLX_LH_TOT")[1:nsteps]  #latent heat
	FGR            = get.var.ncdf(Data.clm, "FGR")[1:nsteps]          #ground heat flux (inc. snowmelt)[1:nsteps]
	FGR12          = get.var.ncdf(Data.clm, "FGR12")[1:nsteps]        #ground heat flux (between layers 1-2)[1:nsteps]
	FSDS           = get.var.ncdf(Data.clm, "FSDS")[1:nsteps]         #incident solar
	FIRA           = get.var.ncdf(Data.clm, "FIRA")[1:nsteps]         #net longwave
  FPSN           = get.var.ncdf(Data.clm, "FPSN")[1:nsteps]         #photosynthesis
  GPP            = get.var.ncdf(Data.clm, "GPP")[1:nsteps]          #gross primiary productivity
  NEE            = get.var.ncdf(Data.clm, "NEE")[1:nsteps]          #net ecosystem exchange
  ELAI           = get.var.ncdf(Data.clm, "ELAI")[1:nsteps]         #leaf area
  FSR            = get.var.ncdf(Data.clm, "FSR")[1:nsteps]          #reflected solar
	FSA            = get.var.ncdf(Data.clm, "FSA")[1:nsteps]          #absorbed solar
  TSOI           = get.var.ncdf(Data.clm, "TSOI")[,1:nsteps]         #soil temperature (K)[1:nsteps]
  SOILLIQ        = get.var.ncdf(Data.clm, "SOILLIQ")[,1:nsteps]       	#Volumetric soil moisture
  H2OSOI         = get.var.ncdf(Data.clm, "H2OSOI")[,1:nsteps]	     	#Volumetric soil moisture
  SOILPSI         = get.var.ncdf(Data.clm, "SOILPSI")[,1:nsteps]       	#Volumetric soil moisture
  SNOW   		     = get.var.ncdf(Data.clm, "SNOW")[1:nsteps]       
  RAIN     	     = get.var.ncdf(Data.clm, "RAIN")[1:nsteps]       
  SNOW_DEPTH     = get.var.ncdf(Data.clm, "SNOW_DEPTH")[1:nsteps]	#snow height of snow covered area
  QVEGT          = get.var.ncdf(Data.clm, "QVEGT")[1:nsteps]      # mm/s      
  QVEGE          = get.var.ncdf(Data.clm, "QVEGE")[1:nsteps]       
  QSOIL          = get.var.ncdf(Data.clm, "QSOIL")[1:nsteps]       
  FCEV           = get.var.ncdf(Data.clm, "FCEV")[1:nsteps]       #canopy evaporation (w/m2)[1:nsteps]
  FCTR           = get.var.ncdf(Data.clm, "FCTR")[1:nsteps]       #canopy transpiration (w/m2)[1:nsteps]
  FGEV           = get.var.ncdf(Data.clm, "FGEV")[1:nsteps]       #ground evaporation (w/m2)
  att.get.ncdf(Data.clm, "RH", "units")
  att.get.ncdf(Data.clm, "RH_LEAF", "units")
  GPPunits <- att.get.ncdf(Data.clm, "GPP", "units")
  GPPunits
  mean(GPP) * 3600 * 24 * 365    #gC/m^s/y
  SOILPSIunits <- att.get.ncdf(Data.clm, "SOILPSI", "units")
  SOILLIQunits <- att.get.ncdf(Data.clm, "SOILLIQ", "units")
  H2OSOIunits  <- att.get.ncdf(Data.clm, "H2OSOI",  "units")
  SOILLIQunits <- att.get.ncdf(Data.clm, "SOILLIQ", "units")
  SOILPSIunits
  SOILLIQunits
  mean(H2OSOI)
  BTRAN[GPP<=0] <- NA      #so only + values count in the calculation
  FPG[  GPP<=0] <- NA
  Rn            <- FSA - FIRA		
  TSOI          <- TSOI - 273.15
	nsteps        <- length(MCDATE)
	MCDATE[nsteps]
	MCSEC[2]
	MCDATE[nsteps-48]
  n1 <- 365*48  #for no-leap   data
  n2 <- n1 + 48 #for leap year data (leap)

  Rn[Rn < -9000 ] <- NA
  plot(Rn ~ SNOW_DEPTH)
  
  shift      <- 15  #adjust CLM to MST  
  mcdate     <- MCDATE[1:(nsteps-shift+1)]
  mcsec      <- MCSEC[1:(nsteps-shift+1)]
  mchour     <- mcsec/3600
  rn         <- Rn[1+shift:nsteps-1]       #omit last day (Jan 1, 2013)
  gpp        <- GPP[1+shift:nsteps-1]      #omit last day (Jan 1, 2013)
  nee        <- NEE[1+shift:nsteps-1]      #omit last day (Jan 1, 2013)
  tsoi       <- TSOI[3,1+shift:nsteps-1]   #select 3rd soil layer
  snow       <- SNOW[1+shift:nsteps-1]
  rain       <- RAIN[1+shift:nsteps-1]
  snow_depth <- SNOW_DEPTH[1+shift:nsteps-1]
  h2osoi     <- H2OSOI[,1+shift:nsteps-1]	
  soilliq    <- SOILLIQ[,1+shift:nsteps-1]  
  soilpsi    <- SOILPSI[,1+shift:nsteps-1]  
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
  SNOW_DEPTH_mean <- tapply(snow_depth,day, mean)
  nlevsoi      <- 5
  ndaily       <- length(DATE_mean)

  H2OSOI_mean  <- array(NA, c(nlevsoi,ndaily))  
  SOILLIQ_mean <- array(NA, c(nlevsoi,ndaily))	
  SOILPSI_mean <- array(NA, c(nlevsoi,ndaily))  
  for (i in 1:nlevsoi) {
  	H2OSOI_mean[i,]   <- tapply(h2osoi[i,],day, mean)     	
  	SOILLIQ_mean[i,]  <- tapply(soilliq[i,],day, mean)     	
  	SOILPSI_mean[i,]  <- tapply(soilpsi[i,],day, mean)     	
  }
  PPT_sum      <- SNOW_sum + RAIN_sum
#  close.ncdf(Data.clm)
print(mean(GPP)*3600*24*365)
print(mean(LH))
case[e]
  #-------------------------------------------------------------------------------
  #---------------read in Observations -------------------------------------------
  #-------------------------------------------------------------------------------
  #              DAILY SNOW DEPTH MEASURESMENTS
  #-------------------------------------------------------------------------------
  
  data.sno <- read.csv(paste(dir,'SnowDepth_daily.csv',sep=''))
  names(data.sno)
  #-------------------------------------------------------------------------------
  #              FLUX TOWER MEASURESMENTS
  #-------------------------------------------------------------------------------
  #  Load data with one header and one unit row from (tab-delimited) text file
  Data.flx <- fLoadTXTIntoDataframe(paste(dir,'Tvan_flux_OBS_2008-2013b.txt',sep=''))
  names(Data.flx)
  #  units(Data.flx)
  maxYear <- max(Data.flx$Year)

  #read in ground flux measurements:
  dir2 <- '/Users/wwieder/Desktop/Working_files/Niwot/NR_fluxes/TVan/G/'
  Data.G  <- fLoadTXTIntoDataframe(paste(dir2,'Wieder_Niwot_CLM_G.txt',sep=''))
  Gobs    <- Data.G$G[Data.G$Year <= maxYear]

  #check that data are same dimensions
  length(Gobs) / (365*48) 
  length(Data.flx$Rn) / (365 * 48)
  length(Data.flx$Rn)

  plot(Data.flx$Rn[1:(nsteps-1)], Rn[1:(nsteps-1)], pch=16, cex=0.2)
  n <- 24
  plot(Data.flx$SoilMoisture[Data.flx$Year == 2010] ~ Data.flx$DecimalDate[Data.flx$Year == 2010], 
       ylim=c(0,35),type='l')
  lines(SOILLIQ[2,Data.flx$Year == 2010]~ Data.flx$DecimalDate[Data.flx$Year == 2010], col=2)
  lines(SOILLIQ[4,Data.flx$Year == 2010]~ Data.flx$DecimalDate[Data.flx$Year == 2010], col=4)
  plot(Data.flx$SoilMoisture[Data.flx$Year == 2011] ~ Data.flx$DecimalDate[Data.flx$Year == 2011], 
       ylim=c(0,30),type='l')
  lines(SOILLIQ[3,Data.flx$Year == 2011]~ Data.flx$DecimalDate[Data.flx$Year == 2011], col=4)
  lines(SOILLIQ[2,Data.flx$Year == 2011]~ Data.flx$DecimalDate[Data.flx$Year == 2011], col=2)
  min(Data.flx$Rn)

#_____HERE____  
  plot(Data.flx$Reco[(160*48):(175*48)])
  plot(Data.flx$GPP_f[(160*48):(175*48)], col=2)
  lines(GPP[(160*48):(175*48)]*1000, col=4)
  abline(h=0)

plot(Data.flx$SoilMoisture, type='l',col=1, ylim=c(0,35), 
     ylab=paste('SOILLIQ[3]',  SOILLIQunits))
#lines(SOILLIQ[3,], col=4 )
#lines(SOILLIQ[5,], col=5 )
lines(SOILLIQ[3,], col=2)
  
#Remove negative GPP values, add to Reco
  Reco  <- Data.flx$Reco
  GPP_f <- Data.flx$GPP_f 
  range(Reco , na.rm=T)
  range(GPP_f, na.rm=T)

  Reco[!is.na(GPP_f) & GPP_f<0]  <- (Reco[!is.na(GPP_f) & GPP_f<0] - GPP_f[!is.na(GPP_f) & GPP_f<0])
  GPP_f[!is.na(Data.flx$Reco) & Data.flx$Reco <0] <- (GPP_f[!is.na(Data.flx$Reco) & Data.flx$Reco<0] - 
                                                           Data.flx$Reco[!is.na(Data.flx$Reco) & Data.flx$Reco<0])
  GPP_f[Data.flx$GPP_f<0] <- 0. 
  Reco[Data.flx$Reco<0]   <- 0. 
  mean(GPP_f, na.rm=T)
  mean(Data.flx$GPP_f, na.rm=T)
  range(GPP_f, na.rm=T)
  range(Data.flx$GPP_f , na.rm=T)
  y  <- 365*48*1
  di <- y+200*48
  df <- y+203*48
  plot(GPP_f[di:df] )
  lines(GPP[di:df]*1000, col=2)
  lines(Data.flx$Reco[di:df], col=4)
#  Data.flx$DD
#  Data.flx$GPP_f[Data.flx$GPP_fqc == 1] <- NA 
  cor.test(GPP[1:(nsteps-1)], GPP_f[1:(nsteps-1)])
  Data.flx$HOUR <- Data.flx$HR*100 + Data.flx$MM
  flxHOUR <- levels(as.factor(Data.flx$HOUR))  
mean(Data.flx$GPP_f[Data.flx$GPP_fqc <= 1], na.rm=T) * 1e-3 * 3600 * 24 * 365 * 12/44
mean(Data.flx$GPP_f, na.rm=T)* 1e-3 * 3600 * 24 * 365 * 12/44
length(Data.flx$GPP_f[Data.flx$GPP_fqc < 1])

#-----------------------------------------------------------------------
# plot seasonal energy fluxes
#----------------------------------------------------------------------
nCLM <- length(Rn)

#Summary Stats for results
cor.test(Data.flx$Rn_MDS , Rn)[[4]][[1]]
cor.test(Data.flx$LE , LH)[[4]][[1]]
cor.test(Data.flx$H , FSH)[[4]][[1]]
cor.test(Gobs, FGR)[[4]][[1]]
clmdate   <- as.Date(as.character(mcdate),format='%Y%m%d')
clmMONTH  <- as.numeric(format.Date(clmdate, "%m"))
clmdate[1]

season <- c('DJF', 'MAM', 'JJA', 'SON')
fout   <- paste(path,'Fig0_Seasonal_energy_fluxes_Tvan',case[e],'.pdf', sep='')
pdf(fout, width=7, height=7, compress = FALSE)

par(mfrow=c(4,4),mar=c(0,0,0,4), oma=c(5,5,3,2))
length(Rn)  
for (i in 1:4) {  
  if (i == 1) { m <- c(1,2,12)   }
  if (i == 2) { m <- c(3,4,5)    }
  if (i == 3) { m <- c(6,7,8)    }
  if (i == 4) { m <- seq(9,10,11) }
  
  tempRNobs  <- Data.flx$Rn_MDS[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3]] 
  tempFGRobs <- Gobs[           Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3]] 
  tempLHobs  <- Data.flx$LE[    Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3]] 
  tempFSHobs <- Data.flx$H[     Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3]] 
  tempGPPobs <- Data.flx$GPP_f[ Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3]] 
  tempGPPobs <- tempGPPobs *  12/44    #mgC/m2/s from mgCO2/m2/s
  range(tempGPPobs, na.rm=T)
  tempRNclm  <- Rn[ Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3]]
  tempLHclm  <- LH[ Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3]]
  tempFSHclm <- FSH[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3]]
  tempFGRclm <- FGR[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3]]
  tempGPPclm <- GPP[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3]]
  tempDATclm <- MCDATE[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3]]
  nall <- length(tempRNobs)  
  tempRNclm  <- tempRNclm[1:nall]   #omit Jan 1 2014
  tempLHclm  <- tempLHclm[1:nall]
  tempFSHclm <- tempFSHclm[1:nall]
  tempFGRclm <- tempFGRclm[1:nall]
  tempGPPclm <- tempGPPclm[1:nall] * 1e3 #mgC/m2/s
  range(tempRNobs, na.rm=T)
  range(tempFSHclm, na.rm=T)
  
  lim  <- c(-350,1200)   #c(-150,1150)
  lim2 <- c(-350,700 )
  lim3 <- c(-0.2,0.3 )
  
#  par(mar=c(0,0,0,4))
  plot(tempRNclm[1:(nall-1)] , tempRNobs[1:(nall-1)], 
       ylim=lim, xlim=lim, pch=".",
       xaxt = "n", yaxt = "n") 
  abline(0,1, lty=2)
  axis(2)
  text(500,1100, season[i], cex=1.5)
  if (i == 1) {text(-350,1100,'(a)', cex=1.5, adj=0)} 
  if (i == 4) {axis(1)}
  if (i == 4) {mtext ("CLM Rnet", side = 1, line = 3)}
  if (i == 3) {mtext (expression(paste("Observed Net Radiation (W ",m^-2,")")), side = 2, line = 2.1, adj = 0)}
  tempcor <- cor.test(tempRNobs , tempRNclm)[[4]][[1]]
  text(900,-200, paste("r =",signif(tempcor[1], digits=2)))

#  par(mar=c(0,0,0,4), oma=c(5,5,3,2))
  plot(tempFSHclm, tempFSHobs, 
       ylim=lim2 , xlim=lim2 , pch=".",
       xaxt = 'n', yaxt='n')
  abline(0,1, lty=2)
  axis(2)
  if (i == 1) {text(-350,625,'(b)', cex=1.5, adj=0)} 
  if (i == 4) {axis(1)}
  if (i == 4) {mtext ("CLM H", side = 1, line = 3)}
  if (i == 3) {mtext (expression(paste("Observed Sensible Heat (W ",m^-2,")")), side = 2, line = 2.1, adj = 0)}
  tempcor <- cor.test(tempFSHobs , tempFSHclm)[[4]][[1]]
  text(500,-250, paste("r =",signif(tempcor[1], digits=2)))
  
  #  par(mar=c(4,0,0,0), oma=c(5,5,3,2))
  plot(tempLHclm, tempLHobs, 
       ylim=lim2, xlim=lim2, pch=".",
       xaxt = 'n', yaxt='n')
  abline(0,1, lty=2)
  axis(2)
  if (i == 1) {text(-350,625,'(c)', cex=1.5, adj=0)} 
  if (i == 4) {axis(1)}
  if (i == 4) {mtext ("CLM LH", side = 1, line = 3)}
  if (i == 3) {mtext (expression(paste("Observed Latent Heat (W ",m^-2,")")), side = 2, line = 2.1, adj = 0)}
  tempcor <- cor.test(tempLHobs , tempLHclm)[[4]][[1]]
  text(500,-250, paste("r =",signif(tempcor[1], digits=2)))

plot(tempFGRclm, tempFGRobs, 
     ylim=lim2, xlim=lim2, pch=".",
     xaxt = 'n', yaxt='n')
abline(0,1, lty=2)
axis(2)
if (i == 1) {text(-350,625,'(c)', cex=1.5, adj=0)} 
if (i == 4) {axis(1)}
if (i == 4) {mtext ("CLM FGR", side = 1, line = 3)}
if (i == 3) {mtext (expression(paste("Observed Ground Heat (W ",m^-2,")")), side = 2, line = 2.1, adj = 0)}
tempcor <- cor.test(tempFGRobs , tempFGRclm)[[4]][[1]]
text(500,-250, paste("r =",signif(tempcor[1], digits=2)))



#  par(mar=c(4,0,0,0), oma=c(5,5,3,2))
#  plot(tempGPPclm, tempGPPobs, 
#       ylim=lim3, xlim=lim3, pch=".",
#       xaxt = 'n', yaxt='n')
#  abline(0,1, lty=2)
#  axis(2)
#  if (i == 4) {axis(1)}
#  if (i == 4) {mtext ("CLM GPP", side = 1, line = 3)}
#  if (i == 3) {mtext (expression(paste("Observed GPP (gC ",m^-2," ",s^-1,")")), 
#                      side = 2, line = 2.1, adj = 0)}
#  tempcor <- cor.test(tempGPPobs , tempGPPclm)[[4]][[1]]
#  text(0.2,-0.1, paste("r =",signif(tempcor[1], digits=2)))

}   # close seasonal loop  
dev.off()
  

#-----------------------------------------------------------------------
# look at Bowen ratio (Qh/Qe) and evaporative fraction (1/(1+B))
# here monthly sum of obs & modeled fluxes 
#-----------------------------------------------------------------------
year   <- seq(2008,2013,1) 
nyear  <- length(year)
dim    <- c(nyear,12)
obsRN_mon  <- array(NA, dim)
obsH_mon   <- array(NA, dim)
obsE_mon   <- array(NA, dim)
obsNEE_mon <- array(NA, dim)
obsGPP_mon <- array(NA, dim)
clmRN_mon  <- array(NA, dim)
clmH_mon   <- array(NA, dim)
clmE_mon   <- array(NA, dim)
clmGPP_mon <- array(NA, dim)

for (i in 1:nyear) {
  for (j in 1:12) {
    obsRN_mon[i,j]  <- sum( Data.flx$Rn_MDS[  Data.flx$MO == j & Data.flx$Year == year[i] ], na.rm=T)
    obsH_mon[i,j]   <- sum( Data.flx$H[  Data.flx$MO == j & Data.flx$Year == year[i] ], na.rm=T)
    obsE_mon[i,j]   <- sum( Data.flx$LE[ Data.flx$MO == j & Data.flx$Year == year[i]], na.rm=T)
    obsNEE_mon[i,j] <- mean(Data.flx$NEE[Data.flx$MO == j & Data.flx$Year == year[i]], na.rm=T)
    obsGPP_mon[i,j] <- mean(Data.flx$GPP_f[Data.flx$MO == j & Data.flx$Year == year[i]], na.rm=T)
    clmH_mon[i,j]   <- sum( FSH[Data.flx$MO == j & Data.flx$Year == year[i]], na.rm=T)
    clmG_mon[i,j]   <- sum( FGR[Data.flx$MO == j & Data.flx$Year == year[i]], na.rm=T)
    clmRN_mon[i,j]  <- sum(  Rn[Data.flx$MO == j & Data.flx$Year == year[i]], na.rm=T)
    clmE_mon[i,j]   <- sum( LH[Data.flx$MO == j & Data.flx$Year == year[i]], na.rm=T)
    clmGPP_mon[i,j] <- mean(GPP[Data.flx$MO == j & Data.flx$Year == year[i]], na.rm=T)
  }
}


obs_B <- obsH_mon / obsE_mon
obsEF <- 1/(1+obs_B)
clm_B <- clmH_mon / clmE_mon
clmEF <- 1/(1+clm_B)


plot(colMeans(clm_B), type="l",col=2,lwd=3)
lines(colMeans(obs_B), col=1, lwd=3)
plot(colMeans(obsEF), type="l",col=1,lwd=3)
lines(colMeans(clmEF), col=2, lwd=3)

data.frame(colMeans(obsGPP_mon), colMeans(obsEF), apply(obsEF, 2, sd)) 
data.frame(colMeans(clmGPP_mon),colMeans(clmEF), apply(obsEF, 2, sd))

#growing season EF
data.frame(mean(obsEF[,5:9]), sd(obsEF[,5:9])) 
data.frame(mean(clmEF[,5:9]), sd(clmEF[,5:9])) 



  

#---------------------------------------------------------
#---------------------------------------------------------
#calculate diurnal cycle seasonally
#---------------------------------------------------------
#---------------------------------------------------------

# ---- NA model where no obs ---- 
  GPP[is.na(Data.flx$GPP_f)] <- NA
  LH[is.na(Data.flx$LE)]     <- NA
  Rn[is.na(Data.flx$Rn_MDS)] <- NA
  FSH[is.na(Data.flx$H)]     <- NA
  LH[is.na(Data.flx$LE)]     <- NA
  FGR[is.na(Gobs)]           <- NA

  for (i in 1:nlevsoi) {
    H2OSOI[i,is.na(Data.flx$SoilMoisture)] <-  NA
    SOILLIQ[i,is.na(Data.flx$SoilMoisture)] <-  NA
  }
  min(H2OSOI[1,])


  Data.flx$hour<- as.numeric(paste(Data.flx$HR,".",Data.flx$MM, sep=""))
  Data.flx$HOUR<- Data.flx$hour - 7.3 
  Data.flx$HOUR[ Data.flx$HOUR < 0] <-   Data.flx$HOUR[Data.flx$HOUR < 0] + 24
  clmmon      <- as.numeric(as.character(clmdate, format = "%m")) 

  OBS_TIM_mean_mo   <- array(NA, c(4,48))
  OBS_FSH_mean_mo   <- array(NA, c(4,48))
  OBS_FGR_mean_mo   <- array(NA, c(4,48))
  OBS_RN_mean_mo    <- array(NA, c(4,48))
  OBS_LH_mean_mo    <- array(NA, c(4,48))
  OBS_GPP_mean_mo   <- array(NA, c(4,48))
  OBS_FSH_sd_mo   <- array(NA, c(4,48))
  OBS_FGR_sd_mo   <- array(NA, c(4,48))
  OBS_RN_sd_mo    <- array(NA, c(4,48))
  OBS_LH_sd_mo    <- array(NA, c(4,48))
  OBS_GPP_sd_mo     <- array(NA, c(4,48))
  OBS_TAIR_mean_mo  <- array(NA, c(4,48))
  OBS_SOILLIQ_mean_mo<- array(NA, c(4,48))

  CLM_TIM_mean_mo   <- array(NA, c(4,48))
  CLM_FPSN_mean_mo  <- array(NA, c(4,48))
  CLM_GPP_mean_mo   <- array(NA, c(4,48))
  CLM_RN_mean_mo    <- array(NA, c(4,48))
  CLM_FSH_mean_mo   <- array(NA, c(4,48))
  CLM_FGR_mean_mo   <- array(NA, c(4,48))
  CLM_FGR12_mean_mo <- array(NA, c(4,48))
  CLM_LH_mean_mo    <- array(NA, c(4,48))
  CLM_RN_sd_mo    <- array(NA, c(4,48))
  CLM_FSH_sd_mo   <- array(NA, c(4,48))
  CLM_FGR_sd_mo   <- array(NA, c(4,48))
  CLM_FGR12_sd_mo <- array(NA, c(4,48))
  CLM_LH_sd_mo    <- array(NA, c(4,48))
  CLM_GPP_sd_mo     <- array(NA, c(4,48))
  CLM_QVEGT_mean_mo <- array(NA, c(4,48))
  CLM_QVEGE_mean_mo <- array(NA, c(4,48))
  CLM_QSOIL_mean_mo <- array(NA, c(4,48))
  CLM_FCTR_mean_mo  <- array(NA, c(4,48))
  CLM_FCEV_mean_mo  <- array(NA, c(4,48))
  CLM_FGEV_mean_mo  <- array(NA, c(4,48))
  CLM_TV_mean_mo    <- array(NA, c(4,48))
  CLM_T10_mean_mo   <- array(NA, c(4,48))
  CLM_RH_mean_mo    <- array(NA, c(4,48))
  CLM_RH_LEAFmean_mo<- array(NA, c(4,48))
  CLM_BTRAN_mean_mo <- array(NA, c(4,48))
  CLM_SOILPSI1_mean_mo <- array(NA, c(4,48))
  CLM_SOILPSI2_mean_mo <- array(NA, c(4,48))
  CLM_SOILPSI3_mean_mo <- array(NA, c(4,48))
  CLM_SOILPSI4_mean_mo <- array(NA, c(4,48))
  CLM_SOILLIQ1_mean_mo <- array(NA, c(4,48))
  CLM_SOILLIQ2_mean_mo <- array(NA, c(4,48))
  CLM_SOILLIQ3_mean_mo <- array(NA, c(4,48))
  CLM_SOILLIQ4_mean_mo <- array(NA, c(4,48))

  SOILPSI1 <- SOILPSI[1,]
  SOILPSI2 <- SOILPSI[2,]
  SOILPSI3 <- SOILPSI[3,]
  SOILPSI4 <- SOILPSI[4,]

  SOILLIQ1 <- SOILLIQ[1,]
  SOILLIQ2 <- SOILLIQ[2,]
  SOILLIQ3 <- SOILLIQ[3,]
  SOILLIQ4 <- SOILLIQ[4,]

  for (i in 1:4) {  
    if (i == 1) { m <- c(1,2,12)   }
    if (i == 2) { m <- c(3,4,5)    }
    if (i == 3) { m <- c(6,7,8)    }
    if (i == 4) { m <- c(9,10,11) }
      
    
    OBS_TIM_mean_mo[i,] <- tapply(Data.flx$hour[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3] ], 
                                  Data.flx$hour[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3] ], 
                                  mean, na.rm=T)
    OBS_SOILLIQ_mean_mo[i,] <- tapply(Data.flx$SoilMoisture[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3] ],
                                      Data.flx$hour[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3] ], 
                                      mean, na.rm=T) 
    OBS_TAIR_mean_mo[i,] <- tapply(Data.flx$Tair[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3] ], 
                                   Data.flx$hour[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3] ], 
                                   mean, na.rm=T) + 273
    OBS_RN_mean_mo[i,]  <- tapply(Data.flx$Rn_MDS[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3] ], 
                                  Data.flx$hour[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3] ], 
                                  mean, na.rm=T) 
    OBS_FSH_mean_mo[i,] <- tapply(Data.flx$H[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3] ], 
                                  Data.flx$hour[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3] ], 
                                  mean, na.rm=T) 
    OBS_LH_mean_mo[i,]  <- tapply(Data.flx$LE[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3] ], 
                                  Data.flx$hour[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3] ], 
                                  mean, na.rm=T) 
    OBS_FGR_mean_mo[i,]  <- tapply(Gobs[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3] ], 
                                 Data.flx$hour[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3] ], 
                                 mean, na.rm=T) 
    OBS_GPP_mean_mo[i,] <- tapply(Data.flx$GPP_f[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3] ], 
                                  Data.flx$hour[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3] ], 
                                  mean, na.rm=T) *  12/44 #mgC/m2/s
    OBS_RN_sd_mo[i,]  <- tapply(Data.flx$Rn_MDS[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3] ], 
                                Data.flx$hour[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3] ], 
                                sd, na.rm=T) 
    OBS_FSH_sd_mo[i,] <- tapply(Data.flx$H[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3] ], 
                                  Data.flx$hour[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3] ], 
                                  sd, na.rm=T) 
    OBS_LH_sd_mo[i,]  <- tapply(Data.flx$LE[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3] ], 
                                Data.flx$hour[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3] ], 
                                sd, na.rm=T) 
    OBS_FGR_sd_mo[i,]  <- tapply(Gobs[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3] ], 
                                 Data.flx$hour[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3] ], 
                                 sd, na.rm=T) 
    OBS_GPP_sd_mo[i,]   <- tapply(Data.flx$GPP_f[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3] ], 
                                  Data.flx$hour[Data.flx$MO == m[1] | Data.flx$MO == m[2] | Data.flx$MO == m[3] ], 
                                  sd,   na.rm=T) *  12/44 #from mgCO2/m2/s

    CLM_TIM_mean_mo[i,]  <- tapply(MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                   MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T) 
    CLM_FPSN_mean_mo[i,] <- tapply(FPSN[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                   MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T)   
    CLM_GPP_mean_mo[i,]  <- tapply(GPP[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                   MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T)   * 1e3 #mg/m2/s 
    CLM_GPP_sd_mo[i,]    <- tapply(GPP[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                   MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], sd, na.rm=T)     * 1e3
    CLM_RN_mean_mo[i,]   <- tapply(Rn[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                   MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T)   
    CLM_RN_sd_mo[i,]   <- tapply(Rn[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                   MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], sd, na.rm=T)   
    CLM_RH_mean_mo[i,]   <- tapply(RH[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                   MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T)   
    CLM_RH_LEAFmean_mo[i,]<-tapply(RH_LEAF[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                   MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T)   
    CLM_BTRAN_mean_mo[i,]<- tapply(BTRAN[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                   MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T)   

    CLM_SOILLIQ1_mean_mo[i,]<- tapply(SOILLIQ1[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                      MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T)   
    CLM_SOILLIQ2_mean_mo[i,]<- tapply(SOILLIQ2[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                      MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T)   
    CLM_SOILLIQ3_mean_mo[i,]<- tapply(SOILLIQ3[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                      MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T)   
    CLM_SOILLIQ4_mean_mo[i,]<- tapply(SOILLIQ4[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                      MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T)   
    CLM_SOILPSI1_mean_mo[i,]<- tapply(SOILPSI1[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                      MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T)   
    CLM_SOILPSI2_mean_mo[i,]<- tapply(SOILPSI2[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                      MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T)   
    CLM_SOILPSI3_mean_mo[i,]<- tapply(SOILPSI3[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                      MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T)   
    CLM_SOILPSI4_mean_mo[i,]<- tapply(SOILPSI4[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                      MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T)   
    CLM_TV_mean_mo[i,]   <- tapply(TV[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                   MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T)   
    CLM_T10_mean_mo[i,]  <- tapply(T10[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                   MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T)   
    CLM_FSH_mean_mo[i,]  <- tapply(FSH[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                   MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T)   
    CLM_FGR_mean_mo[i,]  <- tapply(FGR[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                   MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T)   
    CLM_FGR12_mean_mo[i,]<- tapply(FGR12[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                   MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T)   
    CLM_LH_mean_mo[i,]   <- tapply(LH[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                   MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T)   
    CLM_FSH_sd_mo[i,]  <- tapply(FSH[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                   MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], sd, na.rm=T)   
    CLM_FGR_sd_mo[i,]  <- tapply(FGR[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                 MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], sd, na.rm=T)   
    CLM_FGR12_sd_mo[i,]<- tapply(FGR12[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                 MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], sd, na.rm=T)   
    CLM_LH_sd_mo[i,]   <- tapply(LH[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                   MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], sd, na.rm=T)   
    CLM_QVEGT_mean_mo[i,]<- tapply(QVEGT[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                   MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T)   
    CLM_QVEGE_mean_mo[i,]<- tapply(QVEGE[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                   MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T)   
    CLM_QSOIL_mean_mo[i,]<- tapply(QSOIL[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                   MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T)   
    CLM_FCTR_mean_mo[i,] <- tapply(FCTR[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                   MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T)   
    CLM_FCEV_mean_mo[i,] <- tapply(FCEV[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                   MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T)   
    CLM_FGEV_mean_mo[i,] <- tapply(FGEV[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], 
                                   MCSEC[clmmon== m[1] | clmmon== m[2] | clmmon == m[3] ], mean, na.rm=T)   
    print(paste(i,'OBS, CLM Rnet=',max(OBS_RN_sd_mo[i,]),max(CLM_RN_sd_mo[i,]), sep=' '))
    
  }

  OBS_TIM_mean_mo[3,]
  OBS_GPP_mean_mo[3,]
  daily <- seq(0,23.5,0.5)

CLM_GPP_sd_mo

# === Plot dirunal fluxes  ===
fout   <- paste(path,'Fig1_Daily_energy_fluxes_Tvan',case[e],'.pdf', sep='')
pdf(fout, width=7, height=7, compress = FALSE)

par(mfrow=c(4,4),mar=c(0,0,0,4), oma=c(5,5,3,2))

for (i in 1:4) {
  # ---- RN ----
  OBSrange <- range(c((OBS_RN_mean_mo + OBS_RN_sd_mo),
                      (OBS_RN_mean_mo - OBS_RN_sd_mo)), na.rm=T)
  CLMrange <- range(c((CLM_RN_mean_mo + CLM_RN_sd_mo),
                      (CLM_RN_mean_mo - CLM_RN_sd_mo)), na.rm=T)
  ylim     <- c(min(c(OBSrange,CLMrange))*0.9, max(c(OBSrange,CLMrange))*1.1)
  OBS_var  <- OBS_RN_mean_mo
  CLM_var  <- CLM_RN_mean_mo
  OBS_var2  <- OBS_RN_sd_mo
  CLM_var2  <- CLM_RN_sd_mo
  yval1   <-c(OBS_var[i,16:48] , OBS_var[i,1:15]) 
  yval1sd <-c(OBS_var2[i,16:48] , OBS_var2[i,1:15]) 
  yval2   <-c(CLM_var[i,16:48] , CLM_var[i,1:15]) 
  yval2sd <-c(CLM_var2[i,16:48] , CLM_var2[i,1:15]) 
  plot(yval1~daily, type = 'l', ylim=ylim,
       yaxt='n', xaxt='n', lwd=4)
  lines(yval2~daily, col=2, lwd=4)
  xx       <- c(daily, rev(daily))
  yy1      <- c((yval1+yval1sd),rev(yval1-yval1sd))
  yy2      <- c((yval2+yval2sd),rev(yval2-yval2sd))
  polygon(xx,yy1,col = rgb(0, 0, 0,0.3), border = NA)
  polygon(xx,yy2,col = rgb(1, 0, 0,0.3), border = NA)
  
  abline(h=0, lty=2)  
  axis(2)
  axTicks(2)  
  text(18,ylim[2]*.90,season[i], cex=1.4) 
  if (i == 1) {text(1,ylim[2]*.90,'(a)', cex=1.5, adj=0)} 
  if (i == 3) {
    mtext (expression(paste("Net Radiation (W ",m^-2,")")), 
           side = 2, line = 2.1, adj = 0)
  }
  if (i == 4) {
    axis(1, at=seq(0,24,6))
    mtext ("Hour", side = 1, line = 2.3)
  }

  
  # ---- FSH ----
  OBSrange <- range(c((OBS_FSH_mean_mo + OBS_FSH_sd_mo),
                      (OBS_FSH_mean_mo - OBS_FSH_sd_mo)), na.rm=T)
  CLMrange <- range(c((CLM_FSH_mean_mo + CLM_FSH_sd_mo),
                      (CLM_FSH_mean_mo - CLM_FSH_sd_mo)), na.rm=T)
  ylim     <- c(min(c(OBSrange,CLMrange)), max(c(OBSrange,CLMrange)))
  OBS_var  <- OBS_FSH_mean_mo
  CLM_var  <- CLM_FSH_mean_mo
  OBS_var2  <- OBS_FSH_sd_mo
  CLM_var2  <- CLM_FSH_sd_mo
  yval1   <-c(OBS_var[i,16:48] , OBS_var[i,1:15]) 
  yval1sd <-c(OBS_var2[i,16:48] , OBS_var2[i,1:15]) 
  yval2   <-c(CLM_var[i,16:48] , CLM_var[i,1:15]) 
  yval2sd <-c(CLM_var2[i,16:48] , CLM_var2[i,1:15]) 
  plot(yval1~daily, type = 'l', ylim=ylim,
       yaxt='n', xaxt='n', lwd=4)
  lines(yval2~daily, col=2, lwd=4)
  yy1      <- c((yval1+yval1sd),rev(yval1-yval1sd))
  yy2      <- c((yval2+yval2sd),rev(yval2-yval2sd))
  polygon(xx,yy1,col = rgb(0, 0, 0,0.3), border = NA)
  polygon(xx,yy2,col = rgb(1, 0, 0,0.3), border = NA)
  abline(h=0, lty=2)
  axis(2)
  axTicks(2)
  if (i == 1) {text(1,ylim[2]*.90,'(b)', cex=1.5, adj=0)} 
  if (i == 3) {    
    mtext(expression(paste("Sensible Heat Flux (W ",m^-2,")")), 
          side = 2, line = 2.1, adj = 0)
  }
  if (i == 4) {
    axis(1, at=seq(0,24,6))
    mtext('Hour', side = 1, line =2.3)
  }

  #---- LH ----
  OBSrange <- range(c((OBS_LH_mean_mo + OBS_LH_sd_mo),
                      (OBS_LH_mean_mo - OBS_LH_sd_mo)), na.rm=T)
  CLMrange <- range(c((CLM_LH_mean_mo + CLM_LH_sd_mo),
                      (CLM_LH_mean_mo - CLM_LH_sd_mo)), na.rm=T)
  ylim     <- c(min(c(OBSrange,CLMrange)), max(c(OBSrange,CLMrange)))
  OBS_var  <- OBS_LH_mean_mo
  CLM_var  <- CLM_LH_mean_mo
  OBS_var2  <- OBS_LH_sd_mo
  CLM_var2  <- CLM_LH_sd_mo
  yval1   <-c(OBS_var[i,16:48] , OBS_var[i,1:15]) 
  yval1sd <-c(OBS_var2[i,16:48] , OBS_var2[i,1:15]) 
  yval2   <-c(CLM_var[i,16:48] , CLM_var[i,1:15]) 
  yval2sd <-c(CLM_var2[i,16:48] , CLM_var2[i,1:15]) 
  plot(yval1~daily, type = 'l', ylim=ylim,
       yaxt='n', xaxt='n', lwd=4)
  lines(yval2~daily, col=2, lwd=4)
  yy1      <- c((yval1+yval1sd),rev(yval1-yval1sd))
  yy2      <- c((yval2+yval2sd),rev(yval2-yval2sd))
  polygon(xx,yy1,col = rgb(0, 0, 0,0.3), border = NA)
  polygon(xx,yy2,col = rgb(1, 0, 0,0.3), border = NA)
  abline(h=0, lty=2)
  axis(2)
  axTicks(2)
  if (i == 1) {text(1,ylim[2]*.90,'(c)', cex=1.5, adj=0)} 
  if (i == 3) {
    mtext(expression(paste("Latent Heat Flux (W ",m^-2,")")), 
          side = 2, line = 2.1, adj = 0)
  }
  if (i == 4) {
    axis(1, at=seq(0,24,6))
    mtext('Hour', side = 1, line =2.3)
  }
  
  #---- FGR ----
#  OBSrange <- range(c((OBS_FGR_mean_mo + OBS_FGR_sd_mo),
#                      (OBS_FGR_mean_mo - OBS_FGR_sd_mo)), na.rm=T)
#  CLMrange <- range(c((CLM_FGR_mean_mo + CLM_FGR_sd_mo),
#                      (CLM_FGR_mean_mo - CLM_FGR_sd_mo)), na.rm=T)
#  ylim     <- c(min(c(OBSrange,CLMrange)), max(c(OBSrange,CLMrange)))
#  OBS_var  <- OBS_FGR_mean_mo
#  CLM_var  <- CLM_FGR12_mean_mo
#  OBS_var2  <- OBS_FGR_sd_mo
#  CLM_var2  <- CLM_FGR12_sd_mo
#  yval1   <-OBS_var[i,] 
#  yval1sd <-OBS_var2[i,] 
#  yval2   <-c(CLM_var[i,16:48] , CLM_var[i,1:15]) 
#  yval2sd <-c(CLM_var2[i,16:48] , CLM_var2[i,1:15]) 
#  yval3   <-c(CLM_FGR_mean_mo[i,16:48], CLM_FGR_mean_mo[i,1:15])
#  plot(yval1~daily, type = 'l', ylim=ylim,
#       yaxt='n', xaxt='n', lwd=4)
#  lines(yval2~daily, col=2, lwd=4)
##  lines(yval3~daily, col=3, lwd=2)
#  yy1      <- c((yval1+yval1sd),rev(yval1-yval1sd))
#  yy2      <- c((yval2+yval2sd),rev(yval2-yval2sd))
#  polygon(xx,yy1,col = rgb(0, 0, 0,0.3), border = NA)
#  polygon(xx,yy2,col = rgb(1, 0, 0,0.3), border = NA)
#  abline(h=0, lty=2)
#  axis(2)
#  axTicks(2)
#  if (i == 1) {text(1,ylim[2]*.90,'(d)', cex=1.5, adj=0)} 
#  if (i == 3) {
#    mtext(expression(paste("Ground Heat Flux_1-2 (W ",m^-2,")")), 
#          side = 2, line = 2.1, adj = 0)
#  }
#  if (i == 4) {
#    axis(1, at=seq(0,24,6))
#    mtext('Hour', side = 1, line =2.3)
#  }
  
  # ---- GPP ----
  OBSrange <- range(c((OBS_GPP_mean_mo + OBS_GPP_sd_mo),
                      (OBS_GPP_mean_mo - OBS_GPP_sd_mo)), na.rm=T)
  CLMrange <- range(c((CLM_GPP_mean_mo + CLM_GPP_sd_mo),
                      (CLM_GPP_mean_mo - CLM_GPP_sd_mo)), na.rm=T)
  ylim     <- c(min(c(OBSrange,CLMrange)), max(c(OBSrange,CLMrange)))
  OBS_var  <- OBS_GPP_mean_mo
  CLM_var  <- CLM_GPP_mean_mo
  OBS_var2  <- OBS_GPP_sd_mo
  CLM_var2  <- CLM_GPP_sd_mo
  yval1   <-c(OBS_var[i,16:48] , OBS_var[i,1:15]) 
  yval1sd <-c(OBS_var2[i,16:48] , OBS_var2[i,1:15]) 
  yval2   <-c(CLM_var[i,16:48] , CLM_var[i,1:15]) 
  yval2sd <-c(CLM_var2[i,16:48] , CLM_var2[i,1:15]) 
  plot(yval1~daily, type = 'l', ylim=ylim,
       yaxt='n', xaxt='n', lwd=4)
  lines(yval2~daily, col=2, lwd=4)
  yy1      <- c((yval1+yval1sd),rev(yval1-yval1sd))
  yy2      <- c((yval2+yval2sd),rev(yval2-yval2sd))
  polygon(xx,yy1,col = rgb(0, 0, 0,0.3), border = NA)
  polygon(xx,yy2,col = rgb(1, 0, 0,0.3), border = NA)
  abline(h=0, lty=2)
  axis(2, at=seq(-0.03,0.06,0.03))
  axTicks(2)
  if (i == 1) {text(1,ylim[2]*.90,'(d)', cex=1.5, adj=0)} 
  if ( i == 3 ) {
    mtext(expression(paste("GPP (gC ",m^-2," ",s^-1,")")), 
          side = 2, line = 2.1, adj = 0)
  }  
  if (i == 4) {
    axis(1, at=seq(0,24,6))
    mtext('Hour', side = 1, line =2.3)
  }
}   # close seasonal loop  


dev.off()


#-----------------------------------------------------------------------
# plot annual GPP values
#-----------------------------------------------------------------------
hour <- seq(0,23.5,0.5)
month  <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
days   <- c( 31  ,  28 ,  31 ,  30 ,   31,   30,   31,   31,   30,   31,   30,  31 ) 
DOY    <- c(15,40,60,75,90,115,150,215,315) # random days to plot below
MONTH   <- c('01','02','03','04','05','06','07','08','09','10','11','12')
nMONTH  <- length(MONTH)
paste(year[1],'-',MONTH[1],'-01',sep='')
clmdate  <- as.Date(as.character(mcdate),format='%Y%m%d')
clmdate[1]
mcdate[1]
s2y          <- 365 * 24 * 3600

ann_obs_GPP    <- tapply(Data.flx$GPP_f, Data.flx$Year, mean, na.rm=T) * s2y * 1e-3 * 12/44    #gC/m2/y from mgCO2/m2/s 
ann_clm_GPP    <- tapply(GPP, Data.flx$Year, mean, na.rm=T) * s2y
cbind(mean(ann_obs_GPP[1:6]),sd(ann_obs_GPP[1:6]))
cbind(mean(ann_clm_GPP[1:6]),sd(ann_clm_GPP[1:6]))

Data.flx$day <- as.character(paste(Data.flx$MO,"-",Data.flx$DD, sep=""), format ="%m-%d")
Data.flx$day <- as.Date(Data.flx$day, format ="%m-%d")
OBS_DAY      <- tapply(as.numeric(Data.flx$day), Data.flx$day, mean, na.rm=T)

OBS_GPP_mean <- tapply(Data.flx$GPP_f, Data.flx$day, mean, na.rm=T)
OBS_GPP_sd   <- tapply(Data.flx$GPP_f, Data.flx$day, sd,   na.rm=T)

OBS_SM_mean  <- tapply(Data.flx$SoilMoisture/100, Data.flx$day, mean, na.rm=T) # what are normal units
OBS_SM_sd    <- tapply(Data.flx$SoilMoisture/100, Data.flx$day, sd, na.rm=T) 

OBS_GPP_mean <- OBS_GPP_mean * 1e-3 * 3600 * 24 * 12/44    #gC/m2/d from mgCO2/m2/s
OBS_GPP_sd   <- OBS_GPP_sd   * 1e-3 * 3600 * 24 * 12/44   #gC/m2/d
OBS_Tair_mean    <- tapply(Data.flx$Tair, Data.flx$day, mean, na.rm=T)

obsJuly <- Data.flx$GPP_f[Data.flx$MO==7] 
obsJuly <- GPP_f[Data.flx$MO==7] 
mean(obsJuly, na.rm=T) * 1e-3 * 3600 * 24 * 12/44
plot(OBS_Tair_mean, type="l")
abline(h=0)

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
month   <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","")
days    <- c( 31  ,  28 ,  31 ,  30 ,   31,   30,   31,   31,   30,   31,   30,  31 ) 
totdays <- rep(0, 13)
for (m in 2:13) { totdays[m] <- sum(days[1:m-1]) }

fout <- paste(path,'Fig2_Tvan',case[e],'_Annual_GPP.pdf', sep='')
pdf(fout, width=7, height=6, compress = FALSE)
#par(mfrow=c(1,1),mar=c(4,4,1,1), oma=c(0,0,0,0),mgp=c(2.5,1,0), cex=1.3)
par(mfrow=c(3,1),mar=c(0,0,0,0), oma=c(5,6,3,2))
plot(OBS_GPP_mean, type="l", ylim=c(-3,6),lwd=3, xaxt='n')
mtext('GPP', side = 2, line = 4)
mtext(expression(paste("(gC ",m^-2," ",d^-1,")")), 
      side = 2, line = 2, cex=0.8) 
text(0,5,'(a)', cex=1.5, adj=0)
x  <- seq(1,length(OBS_DAY),1)
xx <- c(x, rev(x))
yy <- c((OBS_GPP_mean+OBS_GPP_sd),rev(OBS_GPP_mean-OBS_GPP_sd))
polygon(xx,yy,col = rgb(0, 0, 0,0.5), border = NA)
abline(h=0, lty=2) 
clmday      <- (as.character(clmdate, format = "%m-%d")) 
clmyear     <- as.numeric(as.character(clmdate, format = "%Y")) 
years       <- seq(2008,2013,1)

temp <- gpp
tempOBS <- Data.flx$GPP_f[1:length(temp)]
temp[tempOBS==NA] <- NA
GPP_mean <- tapply(temp, clmday, mean, na.rm=T) * 3600 * 24    #gC/m2/d
GPP_sd   <- tapply(temp, clmday, sd, na.rm=T)   * 3600 * 24    #gC/m2/d
x        <- seq(1,length(GPP_mean),1)
xx       <- c(x, rev(x))
yy       <- c((GPP_mean+GPP_sd),rev(GPP_mean-GPP_sd))
polygon(xx,yy,col = rgb(1, 0, 0,0.3), border = NA)
lines(GPP_mean, col=2,lwd=3)
lines(OBS_GPP_mean, col=1,lwd=3)
axis(1, tck = 0.05, at=totdays, labels=NA)

#legend('topleft',legend=c("CLM", "TVan obs"), lwd=4, col=c(2,1),bty='n',cex=1.1)        

print(mean(GPP_mean)*365)
print(mean(OBS_GPP_mean)*365)
# ---- Soil Temperature from Tvan ----

Temp.Path    <- '/Users/wwieder/Desktop/Working_files/Niwot/NR_fluxes/TVan/Soil_temp/'
Temp.FF      <- read.csv(paste(Temp.Path, 'SoilTemp_FF_2.csv' , sep=''))
TEMP_DATE_FF <- as.Date(Temp.FF$Time, "%m/%d/%y %H:%M")
TEMP_YEAR_FF <- as.numeric(as.character(TEMP_DATE_FF, format = "%Y")) 
TEMP_DAY_FF  <- as.character(TEMP_DATE_FF, format = "%m-%d") 

Temp.FF$Temp_avg[Temp.FF$Temp_avg =='NaN']  <- NA
OBS_FF_TSOI   <- Temp.FF$Temp_avg[TEMP_YEAR_FF<=max(years)]
OBS_FF_DAY    <- TEMP_DAY_FF[TEMP_YEAR_FF<=max(years)]
OBS_TSOI_mean <- tapply(OBS_FF_TSOI, OBS_FF_DAY, mean, na.rm=T)
OBS_TSOI_sd   <- tapply(OBS_FF_TSOI, OBS_FF_DAY, sd  , na.rm=T)
TSOI_DAY      <- tapply(OBS_FF_TSOI, OBS_FF_DAY, mean, na.rm=T)

#DAILY TSOI HERE
#par(mar=c(3,4,1,2))
#plot(TEMP_DATE_FF[TEMP_YEAR_FF=='2011'],Temp.FF$Temp_avg[TEMP_YEAR_FF=='2011'], type='l',
#     ylab='soil temp', xlab=NA)
#lines(clmdate[clmyear=='2011'], tsoi[clmyear=='2011'], col=2)
#abline(h=0, lty=2)
mean(ELAI[ELAI>0])

plot(OBS_TSOI_mean, type="l", xlim=c(0,360),ylim=c(-15,25),lwd=3, , xaxt='n')
mtext('Soil Temp', side = 2, line = 4)
mtext(expression(paste('(',~degree~C,', 10 cm)'), sep=''), 
      side = 2, line = 2, cex=0.8) 
text(0,21,'(b)', cex=1.5, adj=0)
x  <- seq(1,length(OBS_DAY),1)
xx <- c(x, rev(x))
yy <- c((OBS_TSOI_mean+OBS_TSOI_sd),rev(OBS_TSOI_mean-OBS_TSOI_sd))
polygon(xx,yy,col = rgb(0, 0, 0,0.5), border = NA)
abline(h=0, lty=2) 

TSOI_mean <- tapply(tsoi, clmday, mean) 
TSOI_sd   <- tapply(tsoi, clmday, sd)   
x         <- seq(1,length(TSOI_mean),1)
xx        <- c(x, rev(x))
yy        <- c((TSOI_mean+TSOI_sd),rev(TSOI_mean-TSOI_sd))
polygon(xx,yy,col = rgb(1, 0, 0,0.3), border = NA)
lines(TSOI_mean, col=2,lwd=3)
lines(OBS_TSOI_mean, col=1,lwd=3)
axis(1, tck = 0.05, at=totdays, labels=NA)

# ---- Soil moisture ----
OBS_SM_mean[OBS_TSOI_mean <= 0] <- NA # mask out missing values
OBS_SM_sd[  OBS_TSOI_mean <= 0] <- NA
plot(OBS_SM_mean, type="l", xlim=c(0,360),ylim=c(0.05,0.25),lwd=3, xaxt='n')
mtext('Soil Moisture', side = 2, line = 4)
mtext(expression(paste("(",mm^3," ",mm^-3,")")), 
      side = 2, line = 2, cex=0.8) 
text(0,0.23,'(c)', cex=1.5, adj=0)
axis(1, tck = 0.05, at=totdays, labels=month, cex=1.2)
#mtext('Day', side = 1, line =2.3)
x  <- seq(1,length(OBS_DAY),1)
xx <- c(x, rev(x))
yy <- c((OBS_SM_mean+OBS_SM_sd),rev(OBS_SM_mean-OBS_SM_sd))
polygon(na.omit(cbind(x = xx, y = yy)),col = rgb(0, 0, 0,0.3), border = NA)
clmday      <- (as.character(clmdate, format = "%m-%d")) 
clmyear     <- as.numeric(as.character(clmdate, format = "%Y")) 
years       <- seq(2008,2013,1)
nyears      <- length(years)

temp <- h2osoi[3,]
tempOBS <- Data.flx$SoilMoisture[1:length(temp)]
temp[tempOBS==NA] <- NA

SM_mean <- tapply(temp, clmday, mean,na.rm=T) 
SM_sd   <- tapply(temp, clmday, sd  ,na.rm=T)   
SM_mean[TSOI_mean <= 0] <- NA 
SM_sd[  TSOI_mean <= 0] <- NA

x       <- seq(1,length(SM_mean),1)
xx      <- c(x, rev(x))
yy      <- c((SM_mean+SM_sd),rev(SM_mean-SM_sd))

polygon(na.omit(cbind(x = xx, y = yy)),col = rgb(1, 0, 0,0.3), border = NA)
lines(SM_mean, col=2,lwd=3)
lines(OBS_SM_mean, col=1,lwd=3)

cor.test(GPP_mean, OBS_GPP_mean)
cor.test(TSOI_mean, OBS_TSOI_mean)
cor.test(SM_mean, OBS_SM_mean)

dev.off()
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

# Calculate annual GPP for obs & model
M_Ogpp <- rep(NA,nyear)
ANNgpp_clm <- rep(NA,nyear)
ANNrain    <- rep(NA,nyear)
for(i in 1:nyear) {
  ytemp <- tapply(Data.flx$GPP_f[Data.flx$Year == year[i]], 
                  Data.flx$day[Data.flx$Year == year[i]], 
                  mean, na.rm=T) * 1e-3 * 3600 * 24 * 12/44 
  
  ztemp <- tapply(Data.flx$Reco[Data.flx$Year == year[i]], 
                  Data.flx$day[Data.flx$Year == year[i]], 
                  mean, na.rm=T) * 1e-3 * 3600 * 24 * 12/44 
  
  ytemp_CLM <- tapply(gpp[clmyear == year[i]], 
                      clmday[clmyear == year[i]], 
                      mean, na.rm=T) * 3600 * 24  # daily CLM GPP flux
  
  
  print(paste(year[i],sum(ytemp), sum(ztemp)))
  may1  <- sum(days[1:4])
  oct1 <- sum(days[1:9])
  ytemp2 <- ytemp[may1:oct1]
  print(paste('MAY-SEPT',year[i],sum(ytemp2[ytemp2>0])))
#  M_Ogpp[i] <- sum(ytemp2[ytemp2>0])
  M_Ogpp[i] <- sum(ytemp2)
  ANNgpp_clm[i] <- sum(ytemp_CLM)
  ANNrain[i]    <- sum(tapply(rain[clmyear == year[i]], 
                              clmday[clmyear == year[i]], 
                              mean, na.rm=T) * 3600 * 24)
  #  if (year[i] <= 2013) { lines(ytemp, lwd=1, col=1) }
}
# WW HERE
print(paste('OBS GPP',mean(M_Ogpp),sd(M_Ogpp)))
print(paste('CLM GPP',mean(ANNgpp_clm),sd(ANNgpp_clm)))
OBS_GPP_mean2 <-OBS_GPP_mean[may1:oct1]
print(sum(OBS_GPP_mean2[OBS_GPP_mean2>0]))

plot(M_Ogpp~ANNgpp_clm)
cbind(M_Ogpp,ANNgpp_clm,ANNrain)
mean(ANNrain[1:3])
mean(ANNrain[4:6])
cor.test(M_Ogpp,ANNgpp_clm)
cor.test(ANNrain,ANNgpp_clm)

par(mfrow=c(1,1),mar=c(4,4,1,1), oma=c(0,0,0,0),mgp=c(2.5,1,0), cex=1.3)
plot(OBS_SOILLIQ_mean, type='l', lwd=3, xlab="Day",
     ylab='soilliq', ylim=c(0,35))
lines(SOILLIQ_mean[1,], lwd=3, col=1, lty=2)
lines(SOILLIQ_mean[3,], lwd=3, col=3)
#lines(SOILLIQ_mean[2,], lwd=3, col=2)
#lines(SOILLIQ_mean[4,], lwd=3, col=4)
lines(SOILLIQ_mean[5,], lwd=3, col=5)
lines(SOILLIQ_mean[6,], lwd=3, col=6)

#lines(H2OSOI_mean[1,]*100, lwd=3, col=1, lty=2)
#lines(H2OSOI_mean[2,]*100, lwd=3, col=2)
#lines(H2OSOI_mean[3,]*100, lwd=3, col=3)
#lines(H2OSOI_mean[4,]*100, lwd=3, col=4)
dev.off()

sum(GPP_mean)
sum(OBS_GPP_mean)

GPP_mean2    <- array(NA, dim=c(nyear,366))
GPP_sd2      <- array(NA, dim=c(nyear,366))
length(clmyear)
names(Data.flx)

# Plot annual cycle or soil moisture


  #----------------  END  ----------------------------------
  #---------------------------------------------------------
  
  
