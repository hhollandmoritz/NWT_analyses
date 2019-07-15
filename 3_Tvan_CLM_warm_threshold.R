# looks at GPP, snow depth, lai, & ecosystem response w/ changes in precip.
# default 1:1 froot:leaf ratio provided HIGH LAI in wet condition
# Chose precipitation level based on max snow depth reported in saddle grid survey for each veg community
remove(list=ls())

  library(ncdf)
  library(ggplot2)
  library(grid)
  library(ggthemes)
  library(boot)
  library(REddyProc)
  library(lattice)
  
  
  # Multiple plot function
  #
  # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
  # - cols:   Number of columns in layout
  # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
  #
  # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  # then plot 1 will go in the upper left, 2 will go in the upper right, and
  # 3 will go all the way across the bottom.
  #
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
      
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }

#used to define decimal places in plots
fmt_dcimals <- function(decimals=0){
  function(x) format(x,nsmall = decimals,scientific = FALSE)
}


#------------------------------------------------------------------------------------  
# --- to use PTrespmods --- 
dir      <- ('/Users/wwieder/Desktop/Working_files/Niwot/NR_fluxes/TVan/Tvan_PTCLM45/')
path     <- paste(dir,'CLM_nc_files/PTrespmods/lowWATSAT/',sep='')
pre      <- ('Tvan_PTrespmods_allPHYS_noLW_lowVCmax_')
suf      <- ('.clm2.h1.2008-01-01-00000.nc')
case     <- c('010dry_ff_noRHpsn_lowWATSAT_lowRESIST_run','010dry_dm_noRHpsn_lowWATSAT_lowRESIST_run',
              '100_mm_noRHpsn_lowRESIST_run','075_wm_noRHpsn_lowRESIST_run',
              '200_sb_noRHpsn_lowWATSAT_lowRESIST_run',
              '010dry_ff_noRHpsn_lowWATSAT_lowRESIST_M-Awarm','010dry_dm_noRHpsn_lowWATSAT_lowRESIST_M-Awarm',
              '100_mm_noRHpsn_lowRESIST_M-Awarm','075_wm_noRHpsn_lowRESIST_M-Awarm',
              '200_sb_noRHpsn_lowWATSAT_lowRESIST_M-Awarm',
              '010dry_ff_noRHpsn_lowWATSAT_lowRESIST_BlackSand','010dry_dm_noRHpsn_lowWATSAT_lowRESIST_BlackSand',
              '100_mm_noRHpsn_lowRESIST_BlackSand','075_wm_noRHpsn_lowRESIST_BlackSand',
              '200_sb_noRHpsn_lowWATSAT_lowRESIST_BlackSand')

  nrows    <- length(case)
  years    <- seq(2008,2013,1)
  nyears   <- length(years)

#create arrays to store annual data  
  rows     <- case
  nrows    <- length(rows)
  cols     <- c(as.character(years))
  dims     <- c(nrows,nyears)
#  NEP_a    <- array(NA,dim=dims,dimnames=list(rows,cols)) # total NEP (gC/m2/s) pos. for sink
  NEE_a    <- array(NA,dim=dims,dimnames=list(rows,cols)) # total NEE (gC/m2/s)
  NPP_a    <- array(NA,dim=dims,dimnames=list(rows,cols)) # total NPP (gC/m2/s)
  ANPP_a   <- array(NA,dim=dims,dimnames=list(rows,cols)) # Aboveground NPP (gC/m2/s)
  BNPP_a   <- array(NA,dim=dims,dimnames=list(rows,cols)) # Belowground NPP (gC/m2/s)
  GPP_a    <- array(NA,dim=dims,dimnames=list(rows,cols)) # total GPP (gC/m2/s)
  FPSN_a   <- array(NA,dim=dims,dimnames=list(rows,cols)) # photosynthesis
  TLAI_a   <- array(NA,dim=dims,dimnames=list(rows,cols)) # median LAI for each year
  LEAF_f   <- array(NA,dim=dims,dimnames=list(rows,cols)) # first day of growing season
  LEAF_l   <- array(NA,dim=dims,dimnames=list(rows,cols)) # last  day of growing season
  TLAI_f   <- array(NA,dim=dims,dimnames=list(rows,cols)) # first day when LAI > 0
  TLAI_l   <- array(NA,dim=dims,dimnames=list(rows,cols)) # last  day when LAI > 0
  VEGC_a   <- array(NA,dim=dims,dimnames=list(rows,cols)) # TOTVEGC
  TOTC_a   <- array(NA,dim=dims,dimnames=list(rows,cols)) # TOTECOSYSC
  PPT_a    <- array(NA,dim=dims,dimnames=list(rows,cols)) # 'precipitation mm/y'
  RAIN_a   <- array(NA,dim=dims,dimnames=list(rows,cols))
  SNOW_a   <- array(NA,dim=dims,dimnames=list(rows,cols))
  BTRAN_a  <- array(NA,dim=dims,dimnames=list(rows,cols))
  FPG_a    <- array(NA,dim=dims,dimnames=list(rows,cols))
  FPI_a    <- array(NA,dim=dims,dimnames=list(rows,cols))
  TV_g    <- array(NA,dim=dims,dimnames=list(rows,cols)) #growing season TV
  TSOI_g   <- array(NA,dim=dims,dimnames=list(rows,cols)) #growing season TSOI
  TSOI_a   <- array(NA,dim=dims,dimnames=list(rows,cols)) # annual TSOI
  TMIN_a   <- array(NA,dim=dims,dimnames=list(rows,cols)) # min TSOI
  TMAX_a   <- array(NA,dim=dims,dimnames=list(rows,cols)) # max TSOI
  H2OSOI_g <- array(NA,dim=dims,dimnames=list(rows,cols)) #growing season H2OSOI
  SOILLIQ_g<- array(NA,dim=dims,dimnames=list(rows,cols)) #growing season liquid water (gravametric)
  HR_a     <- array(NA,dim=dims,dimnames=list(rows,cols))
  AR_a     <- array(NA,dim=dims,dimnames=list(rows,cols))
  LH_a     <- array(NA,dim=dims,dimnames=list(rows,cols))
  FSH_a    <- array(NA,dim=dims,dimnames=list(rows,cols))
  Snow_Depth_a <- array(NA,dim=dims,dimnames=list(rows,cols))
  Snow_Depth_l <- array(NA,dim=dims,dimnames=list(rows,cols))
  Snow_Depth_f <- array(NA,dim=dims,dimnames=list(rows,cols))
  QRUNOFF_a<- array(NA,dim=dims,dimnames=list(rows,cols))
  YEAR_a   <- array(NA,dim=dims,dimnames=list(rows,cols))
  CASE_a    <- array(NA,dim=dims,dimnames=list(rows,cols))
  
  for (i in 1:nyears) {YEAR_a[,i] <- years[i]}
  for (i in 1:nrows)  {CASE_a[i,] <- case[i]}
#-----------------------------------------------------------------------
#---------------read in CLM variables----------------------------------
#-----------------------------------------------------------------------
for (e in 1:15) {
  infile   <- paste(path,pre,case[e],suf,sep='')
  Data.clm <- open.ncdf(infile)   

#  print(paste("The file has",Data.clm$nvars,"variables"))
#  summary(Data.clm)
#  print(Data.clm)
  MCDATE         = get.var.ncdf(Data.clm, "mcdate") # getting/naming the variable
  MCSEC          = get.var.ncdf(Data.clm, "mcsec") 
  BTRAN          = get.var.ncdf(Data.clm, "BTRAN") 
  FPSN           = get.var.ncdf(Data.clm, "FPSN") 
  TLAI           = get.var.ncdf(Data.clm, "TLAI") 
  FPG            = get.var.ncdf(Data.clm, "FPG")      
  FPI            = get.var.ncdf(Data.clm, "FPI")      
  HR             = get.var.ncdf(Data.clm, "HR")      
  AR             = get.var.ncdf(Data.clm, "AR")      
  GPP            = get.var.ncdf(Data.clm, "GPP")      
  NPP            = get.var.ncdf(Data.clm, "NPP")      
  ANPP           = get.var.ncdf(Data.clm, "AGNPP")      
  BNPP           = get.var.ncdf(Data.clm, "BGNPP")      
  NEE            = get.var.ncdf(Data.clm, "NEE")      
  FSH            = get.var.ncdf(Data.clm, "FSH")         #sensible heat
  LH             = get.var.ncdf(Data.clm, "EFLX_LH_TOT") #latent heat
#  NEP            = get.var.ncdf(Data.clm, "NEP")      

  TOTVEGC        = get.var.ncdf(Data.clm, "TOTVEGC")     	
  TOTECOSYSC     = get.var.ncdf(Data.clm, "TOTECOSYSC")       
  TSOI           = get.var.ncdf(Data.clm, "TSOI")	- 273     	#Soil Temp.
  TV             = get.var.ncdf(Data.clm, "TV")   - 273     	#VEG Temp.
  H2OSOI         = get.var.ncdf(Data.clm, "H2OSOI")       #Volumetric soil moisture
  QRUNOFF        = get.var.ncdf(Data.clm, "QRUNOFF")    	#liquid runoff (does not include QSNWCPICE)
  SNOW   		     = get.var.ncdf(Data.clm, "SNOW")       
  RAIN     	     = get.var.ncdf(Data.clm, "RAIN")       
  SNOWDEPTH      = get.var.ncdf(Data.clm, "SNOW_DEPTH")   #gridcell mean snow
  SOILLIQ	       = get.var.ncdf(Data.clm, "SOILLIQ")      #soil liquid water 
  SNOWICE	       = get.var.ncdf(Data.clm, "SNOWICE")      #snow ice water 
  SNOW_DEPTH     = get.var.ncdf(Data.clm, "SNOW_DEPTH")	#snow height of snow covered area
  QRUNOFFunits<- att.get.ncdf(Data.clm, "QRUNOFF", "units")
  QRUNOFFunits
  H2OSOIunits <- att.get.ncdf(Data.clm, "H2OSOI", "units")
  H2OSOIunits
  NEEunits    <- att.get.ncdf(Data.clm, "NEE", "units")
  NEEunits
#  NEPunits    <- att.get.ncdf(Data.clm, "NEP", "units")
#  NEPunits
  PRECIP        <- RAIN + SNOW
  BTRAN[GPP<=0] <- NA      #so only + values count in the calculation
  FPG[  GPP<=0] <- NA

  paste('mean Runoff', mean(QRUNOFF)* 3600 * 24 * 365, 'mm/y')
  paste('mean Precip', mean(PRECIP) * 3600 * 24 * 365, 'mm/y')
  paste('mean GPP',    mean(GPP)    * 3600 * 24 * 365, 'gC/m2/y')
  paste('mean NPP',    mean(NPP)    * 3600 * 24 * 365, 'gC/m2/y')
  paste('mean ANPP',   mean(ANPP)   * 3600 * 24 * 365, 'gC/m2/y')
  paste('mean BNPP',   mean(BNPP)   * 3600 * 24 * 365, 'gC/m2/y')
  paste('mean NEE',    mean(NEE)    * 3600 * 24 * 365, 'gC/m2/y')
#  paste('mean NEP',    mean(NEP)    * 3600 * 24 * 365, 'gC/m2/y')
  paste('mean BTRAN',  mean(BTRAN[GPP>0]) )
  paste('mean FPG',    mean(FPG[GPP>0])   )
  print(paste('mean TSOIg',   mean(TSOI[3,][GPP>0])))
  print(paste('mean H2OSOI', mean(H2OSOI[3,][GPP>0])))
  paste('max SnowDepth' , max(SNOW_DEPTH), 'm')
  
  day2 <- as.Date.character(MCDATE, format='%Y%m%d')
  year <- format(day2,format='%Y')
  mo   <- as.numeric(format(day2,format='%m'))
  doy  <- as.numeric(strftime(day2,format='%j'))
  snowyear <- as.numeric(year)
  snowyear[mo >= 10] <- as.numeric(year[mo >= 10]) + 1
  ascalar      <- 3600*24*365  # converto from gC/m2/s to gC/m2/y
  FPSN_a[e,]   <- (tapply(FPSN     ,   year, mean) )[1:nyears]  
  GPP_a[e,]    <- (tapply(GPP*ascalar,   year, mean) )[1:nyears]       #'gC/m2/y'
  NPP_a[e,]    <- (tapply(NPP*ascalar,   year, mean) )[1:nyears]       #'gC/m2/y'
  ANPP_a[e,]   <- (tapply(ANPP*ascalar,  year, mean) )[1:nyears]       #'gC/m2/y'
  BNPP_a[e,]   <- (tapply(BNPP*ascalar,  year, mean) )[1:nyears]       #'gC/m2/y'
  HR_a[e,]     <- (tapply(HR *ascalar,   year, mean) )[1:nyears]       #'gC/m2/y'
  AR_a[e,]     <- (tapply(AR *ascalar,   year, mean) )[1:nyears]       #'gC/m2/y'
  NEE_a[e,]    <- (tapply(NEE*ascalar,   year, mean) )[1:nyears]       #'gC/m2/y'
#  NEP_a[e,]    <- (tapply(NEP*ascalar,   year, mean) )[1:nyears]       #'gC/m2/y'
  TLAI_a[e,]   <- (tapply(TLAI[TLAI>0], year[TLAI>0], median))[1:nyears]  
  TLAI_f[e,]   <- (tapply(doy[TLAI>0],  year[TLAI>0],min))[1:nyears]  # first day with LAI > 0
  TLAI_l[e,]   <- (tapply(doy[TLAI>0],  year[TLAI>0],max))[1:nyears]  # last day with LAI > 0
  LEAF_f[e,]   <- (tapply(doy[GPP>0],   year[GPP>0],min))[1:nyears]  #first day of growing season
  LEAF_l[e,]   <- (tapply(doy[GPP>0],   year[GPP>0],max))[1:nyears]  # last day of growing season  
  VEGC_a[e,]   <- (tapply(TOTVEGC,        year, max  ))[1:nyears]            
  TOTC_a[e,]   <- (tapply(TOTECOSYSC,     year, max  ))[1:nyears]            
  PPT_a[e,]    <- (tapply(PRECIP*ascalar, year, mean) )[1:nyears]       #'mm/y'
  SNOW_a[e,]   <- (tapply(SNOW*ascalar,   snowyear, mean) )[1:nyears]       #'mm/y'
  RAIN_a[e,]   <- (tapply(RAIN*ascalar,   year, mean) )[1:nyears]       #'mm/y'
  BTRAN_a[e,]  <- (tapply(1-BTRAN[GPP>0], year[GPP>0], mean, na.rm=T))[1:nyears] # reverse so larger values = > limitation
  FPG_a[e,]    <- (tapply(1-FPG[GPP>0],   year[GPP>0], mean, na.rm=T))[1:nyears] 
  FPI_a[e,]    <- (tapply(1-FPI[GPP>0],   year[GPP>0], mean))[1:nyears] 
  TSOI_g[e,]   <- (tapply(TSOI[3,][GPP>0],year[GPP>0],        mean))[1:nyears] 
  TSOI_a[e,]   <- (tapply(TSOI[3,],    year,        mean))[1:nyears] 
  TMAX_a[e,]   <- (tapply(TSOI[3,],    year,        max))[1:nyears]  # max tsoi 
  TMIN_a[e,]   <- (tapply(TSOI[3,],    year,        min))[1:nyears]  # min tsoi 
  Snow_Depth_a[e,] <- (tapply(SNOW_DEPTH, snowyear, max))[1:nyears]         #'m'
  Snow_Depth_l[e,] <- (tapply(doy[SNOW_DEPTH==0], year[SNOW_DEPTH==0], min))[1:nyears]         #'m'
  Snow_Depth_f[e,] <- (tapply(doy[SNOW_DEPTH==0], year[SNOW_DEPTH==0], max))[1:nyears]         #'m'
  QRUNOFF_a[e,]<- (tapply(QRUNOFF*ascalar, snowyear, mean) )[1:nyears] #'mm/y'

  LH_a[e,]    <- (tapply(LH[GPP>0],  year[GPP>0], mean))[1:nyears] 
  FSH_a[e,]   <- (tapply(FSH[GPP>0],  year[GPP>0], mean))[1:nyears] 


  close(Data.clm)
  remove(MCDATE,NEE,NPP,ANPP,BNPP,GPP,PRECIP,SNOW,RAIN,BTRAN,FPG,FPI,AR,HR,SNOW_DEPTH,QRUNOFF,TSOI,H2OSOI,TLAI,FPSN,SOILLIQ)
  print(paste('finished',rows[e]))
}

#-------------------------------------------------------------------------  
# Finished reading and analizing data
#-------------------------------------------------------------------------  
GROW_a  <- LEAF_l - LEAF_f    #growing season length GPP > 0
GROW_b  <- TLAI_l - TLAI_f    #growing season, TLAAI > 0
mean(TSOI_g[6:10, ]/TSOI_g[1:5,])
mean(TSOI_g[11:15,]/TSOI_g[1:5,])
print(c(mean(GPP_a[9,]/GPP_a[4,]-1), mean(GPP_a[14,]/GPP_a[4,]-1)))
print(c(sd(GPP_a[9,]/GPP_a[4,]-1), sd(GPP_a[14,]/GPP_a[4,]-1)))
print(c(mean(AR_a[9,]/AR_a[4,]-1), mean(AR_a[14,]/AR_a[4,]-1)))
print(c(mean(TLAI_a[9,]/TLAI_a[4,]-1), mean(TLAI_a[14,]/TLAI_a[4,]-1)))
print(c(mean(VEGC_a[9,]/VEGC_a[4,]-1), mean(VEGC_a[14,]/VEGC_a[4,]-1)))
print(c(mean(LH_a[9,]/LH_a[4,]),   mean(LH_a[14,]/LH_a[4,])  ))
print(c(mean(FSH_a[9,]/FSH_a[4,]), mean(FSH_a[14,]/FSH_a[4,])))
print(c(mean(GROW_a[9,]-GROW_a[4,]), mean(GROW_a[14,]-GROW_a[4,])))
print(c(sd(GROW_a[9,]-GROW_a[4,]), sd(GROW_a[14,]-GROW_a[4,])))

mean(PPT_a[3,])
sd(PPT_a[3,])
mean(RAIN_a[3,])
mean(SNOW_a[3,])

LEAF_f
TLAI_f
TLAI_l
Snow_Depth_f
LEAF_l
relRAIN   <- RAIN_a/SNOW_a
  veg       <- c(1:5)            #'best' Veg_classification 
  M_Awarm   <- c(6:10)           # spring & summer warmed plots
  BlackSand    <- c(11:15)           # spring & summer warmed plots
#  DJFwarm   <- c(11:15)           # spring & summer warmed plots
#  ndep   <- c(14:18)          # Ndep plots (+5gN/m2/y)
  vegpch <- c(16,17,2,4,1) 
  COL        <- array(NA,dim=dims) 
  COL[veg,]  <- 4
  COL[M_Awarm,] <- 'darkred'
  COL[BlackSand,]  <- 1
  plot(relRAIN~PPT_a, col=COL)
  plot(GPP_a~  PPT_a, col=COL)
  rowMeans(GPP_a)
  vegLegend <- c('FF','DM','MM','WM','SB')

#-------------------------------------------------------------------------  
# Read in OBS data 
#-------------------------------------------------------------------------  
  # NPP estimate from observations (Table 9.1 in Niwot book)
  # different estimates for productivity are for ANPP reported in different papers
  # NOTE OMITTED data from Walker et al. 1994, included data from medicine bow (Scott & Billings)
  # BNPP estimated by mean BNPP/ANPP ratio across all studies & communities
  ANPP_OBS <- c(171,129,
               164,136,134,155,143,
               196,167,262,232,
               139,299,162,291,248,
               124,97,74)
  BNPP_obs <- c(NA  ,NA  ,198 ,198 ,198 ,198 ,198 ,230 ,230 ,230 ,230 ,364 ,364 ,364 ,364 ,364 ,NA  ,NA  , NA )
  VEG_obs  <- c('FF','FF','DM','DM','DM','DM','DM','MM','MM','MM','MM','WM','WM','WM','WM','WM','SB','SB','SB')
  BNPP_OBS <- BNPP_obs
  VEG_OBS  <- factor(VEG_obs, levels=vegLegend)
  
  multiplier <- mean(BNPP_obs/ANPP_OBS, na.rm=T)
  BNPP_OBS[is.na(BNPP_OBS)] <- ANPP_OBS[is.na(BNPP_OBS)] * multiplier
  NPP_OBS <- (ANPP_OBS + BNPP_OBS) / 2                        #convert from gDW/m2/y to gC/m2/y 
  tapply(NPP_OBS, VEG_OBS, mean)
  npp_obs_data     <- data.frame(VEG_OBS,NPP_OBS,ANPP_OBS,BNPP_OBS)
  npp_obs_data$NPP <- npp_obs_data$NPP_OBS
  npp_obs_data$VEG <- npp_obs_data$VEG_OBS
  npp_obs_data$TX  <- rep('Control',length(npp_obs_data$VEG_OBS))
  
  # root biomass estimate from observations (Table 9.2 in Niwot book)
  # unigs g/m2/y?
  # different estimates for reported in different papers
  # Biomass estimates sum of Root & shoot from below:above ground ratios
  ROOT_obs    <- c('FF','DM','DM','DM','MM','MM','WM','WM','WM','SB')
  ROOT_TOTAL  <- c(1858,2929,2772, NA ,2563,3039,3583,5674, NA ,1502) / 2        #convert to gC/m2
  ROOT_LIVE   <- c(1534,1739,2300, 575,2090,1055,2927,3300,1232,1208) / 2        #convert to gC/m2
  ROOT_SHOOTt <- c(2.5 ,3.0 ,4.8 , NA ,2.6 ,3.4 ,5.6 ,8.0 , NA ,5.6 )
  ROOT_SHOOTl <- c(9.0 ,11. ,12. , 3.7,11. ,4.0 ,21. ,9.5 ,4.2 ,10. )  # no Bowman data fro WM, replace w/ mean
  VEG_TOTAL  <- (ROOT_TOTAL + (ROOT_TOTAL/ ROOT_SHOOTt)) 
  VEG_LIVE   <- (ROOT_LIVE  + (ROOT_LIVE / ROOT_SHOOTl)) 
  ROOT_OBS   <- factor(ROOT_obs, levels=vegLegend)
  
  tapply(VEG_LIVE, ROOT_OBS, mean)
  biomass_obs_data      <- data.frame(ROOT_OBS, VEG_TOTAL, VEG_LIVE)
  biomass_obs_data$VEGC <- biomass_obs_data$VEG_LIVE
  biomass_obs_data$VEG  <- biomass_obs_data$ROOT_OBS
  biomass_obs_data$TX   <- rep('Control',length(biomass_obs_data$ROOT_OBS))
  
  #-------------------------------------------------------------------------  
  # ANPP from Saddle Grid
  #-------------------------------------------------------------------------  
  # Plots community veg. data from Jane via Emily Farrer 
  # uses long term community analysis 'class_3'
  data.veg         <- read.csv(paste(dir,'NWT_SnowXprod.csv',sep=''))
  DATA.VEG         <- subset(data.veg,class_3!='rock' & class_3!='SF' & class_3!='ST',select=names(data.veg))
  DATA.VEG         <- subset(DATA.VEG,DATA.VEG$year>=2008 & DATA.VEG$year<=2013,select=names(DATA.VEG))
  good             <- list('FF','DM','MM','WM','SB')
  DATA.VEG$class   <- factor(DATA.VEG$class_3, levels=good)
  VEG.PLOTS        <- as.numeric(levels(as.factor(DATA.VEG$plot)))
  VEG.TYPE         <- DATA.VEG$class[1:length(VEG.PLOTS)]
  VEG.ANPP         <- DATA.VEG$anpp / 2                            #convert biomass to gC/m2
  VEG.YEAR         <- DATA.VEG$year  
  VEG.ANPP_ALLmean <- tapply(VEG.ANPP, INDEX=list(DATA.VEG$class), FUN=mean, na.rm=T)                             #convert biomass to gC/m2
  VEG.ANPP_ALLsd   <- tapply(VEG.ANPP, INDEX=list(DATA.VEG$class), FUN=sd, na.rm=T)                             #convert biomass to gC/m2
  VEG.ANPP_mean    <- tapply(VEG.ANPP, INDEX=list(DATA.VEG$class,VEG.YEAR), FUN=mean, na.rm=T)                             #convert biomass to gC/m2
  VEG.ANPP_sd      <- tapply(VEG.ANPP, INDEX=list(DATA.VEG$class,VEG.YEAR), FUN=sd  , na.rm=T)                             #convert biomass to gC/m2
  VEG.ANPP[VEG.YEAR==2009]

  cor.test(as.vector(VEG.ANPP_mean[1:2,]), as.vector(ANPP_a[1:2,]))
  cor.test(VEG.ANPP_mean[1,], RAIN_a[1,])
  cor.test(NPP_a[1,], RAIN_a[1,])
  
  plot(RAIN_a[1,], VEG.ANPP_mean[1,], pch=16, ylim=c(0,150))
  points(RAIN_a[2,], VEG.ANPP_mean[2,], pch=17)
  points(RAIN_a[1,], ANPP_a[1,], pch=16, col=2)
  points(RAIN_a[2,], ANPP_a[2,], pch=17, col=2)
  points(RAIN_a[1,], NPP_a[1,], pch=1, col=2)
  points(RAIN_a[2,], NPP_a[2,], pch=2, col=2)

  VEG.ANPP_mean[,2] <- NA
  VEG.ANPP_sd[,2]   <- NA
  VEG_ANPP_ALLmean  <- as.vector(VEG.ANPP_ALLmean)
  VEG_ANPP_ALLsd    <- as.vector(VEG.ANPP_ALLsd)
  VEG_ANPP_mean     <- as.vector(VEG.ANPP_mean)
  VEG_ANPP_sd       <- as.vector(VEG.ANPP_sd)
  VEG_CLASS         <- rep(factor(vegLegend, as.ordered(vegLegend)), length(VEG_ANPP_mean)/length(good))
  anpp_obs_data     <- data.frame(VEG_ANPP_mean, VEG_CLASS)
  anpp_obs_data$TX  <- rep('control', length(VEG_ANPP_mean))
  anpp_obs_data$ANPP<- VEG_ANPP_mean
  anpp_obs_data$VEG <- VEG_CLASS
  
  length(VEG.TYPE)
  length(VEG_ANPP_mean)

#-------------------------------------------------------------------------------
#              FLUX TOWER MEASURESMENTS
#-------------------------------------------------------------------------------
#  Load data with one header and one unit row from (tab-delimited) text file
Data.flx <- fLoadTXTIntoDataframe(paste(dir,'Tvan_flux_OBS_2008-2013b.txt',sep=''))
names(Data.flx)
GPP_obs <- Data.flx$GPP_f
s2y     <- 3600 * 24 * 365
ann_obs_GPP    <- tapply(Data.flx$GPP_f, Data.flx$Year, mean, na.rm=T) * s2y * 1e-3 * 12/44    #gC/m2/y from mgCO2/m2/s 


# Calculate MAY-OCT annual GPP for obs & model
Data.flx$day <- as.character(paste(Data.flx$MO,"-",Data.flx$DD, sep=""), format ="%m-%d")
Data.flx$day <- as.Date(Data.flx$day, format ="%m-%d")
days   <- c( 31  ,  28 ,  31 ,  30 ,   31,   30,   31,   31,   30,   31,   30,  31 ) 
may1   <- sum(days[1:4])
oct1   <- sum(days[1:9])
M_Ogpp <- rep(NA,nyears)
for(i in 1:nyears) {
  ytemp <- tapply(Data.flx$GPP_f[Data.flx$Year == years[i]], 
                  Data.flx$day[Data.flx$Year == years[i]], 
                  mean, na.rm=T) * 1e-3 * 3600 * 24 * 12/44 
  
  ytemp2 <- ytemp[may1:oct1]
#  print(paste('MAY-SEPT',years[i],sum(ytemp2[ytemp2>0])))
#  M_Ogpp[i] <- sum(ytemp2[ytemp2>0])
  print(paste('MAY1-Oct1',years[i],sum(ytemp2)))
  M_Ogpp[i] <- sum(ytemp2)
}
print(paste('OBS GPP',mean(M_Ogpp),sd(M_Ogpp)))

gpp_obs_data     <- data.frame(M_Ogpp, rep(VEG_CLASS[1], nyears))
gpp_obs_data$TX  <- rep('control', length(M_Ogpp))
gpp_obs_data$GPP <- M_Ogpp
gpp_obs_data$VEG <- rep(VEG_CLASS[1], nyears)

#----------------------------------------------------------
#reshape data into vectors for ggplot  
#----------------------------------------------------------
TRANGE_a<- abs(TMAX_a) + abs(TMIN_a) 
  SNOW  <- c(as.vector(SNOW_a[veg,]) )#,  as.vector(SNOW_a[M_Awarm,]) )
  RAIN  <- c(as.vector(RAIN_a[veg,]) )#,  as.vector(RAIN_a[M_Awarm,]) )
  PPT   <- c(as.vector(PPT_a[veg,]) )#,   as.vector(PPT_a[M_Awarm,])  )
  BTRAN <- c(as.vector(BTRAN_a[veg,]) )#, as.vector(BTRAN_a[M_Awarm,]))
  FPG   <- c(as.vector(FPG_a[veg,]) )#,   as.vector(FPG_a[M_Awarm,])  )
  FPI   <- c(as.vector(FPI_a[veg,]) )#,   as.vector(FPI_a[M_Awarm,])  )
  HR    <- c(as.vector( HR_a[veg,]) )#,   as.vector( HR_a[M_Awarm,])  )
  AR    <- c(as.vector( AR_a[veg,]) )#,   as.vector( AR_a[M_Awarm,])  )
  GPP   <- c(as.vector(GPP_a[veg,]) )#,   as.vector(GPP_a[M_Awarm,])  )
  NPP   <- c(as.vector(NPP_a[veg,]) )#,   as.vector(NPP_a[M_Awarm,])  )
  ANPP  <- c(as.vector(ANPP_a[veg,]) )#,  as.vector(ANPP_a[M_Awarm,]) )
  BNPP  <- c(as.vector(BNPP_a[veg,]) )#,  as.vector(BNPP_a[M_Awarm,]) )
#  NEP   <- c(as.vector(NEP_a[veg,]) )#,   as.vector(NEP_a[M_Awarm,])  )
  NEE   <- c(as.vector(NEE_a[veg,]) )#,   as.vector(NEE_a[M_Awarm,])  )
  FPSN  <- c(as.vector(FPSN_a[veg,]) )#,  as.vector(FPSN_a[M_Awarm,]) )
  TLAI  <- c(as.vector(TLAI_a[veg,]) )#,  as.vector(TLAI_a[M_Awarm,]) )
  TVg   <- c(as.vector(TV_g[veg,  ]) )#,  as.vector(TV_g[M_Awarm,  ]) )
  TSOIg <- c(as.vector(TSOI_g[veg,]) )#,  as.vector(TSOI_g[M_Awarm,]) )
  TSOIa <- c(as.vector(TSOI_a[veg,]) )#,  as.vector(TSOI_a[M_Awarm,]) )
  TRANGEa<-c(as.vector(TRANGE_a[veg,]) )#,as.vector(TRANGE_a[M_Awarm,]) )
  H2OSOI<- c(as.vector(H2OSOI_g[veg,]) )#,as.vector(H2OSOI_g[M_Awarm,]) )
  VEGC  <- c(as.vector(VEGC_a[veg,]) )#,  as.vector(VEGC_a[M_Awarm,]) )
  TOTC  <- c(as.vector(TOTC_a[veg,]) )#,  as.vector(TOTC_a[M_Awarm,]) )
  CASE  <- c(as.vector(CASE_a[veg,]) )#,  as.vector(CASE_a[M_Awarm,]) )
  YEAR  <- c(as.vector(YEAR_a[veg,]) )#,  as.vector(YEAR_a[M_Awarm,]) )
  GROW  <- c(as.vector(GROW_a[veg,]) )#,  as.vector(GROW_a[M_Awarm,]) )
  TX    <- c(rep("Control",length(CASE_a[veg,])) )#, rep("TMS",length(CASE_a[veg,])))
  VEG   <- c(rep(vegLegend,6))
  data  <- data.frame(YEAR,CASE,TX,VEG,GPP,NPP,ANPP,BNPP,NEE,FPSN,TLAI,TSOIg,TSOIa,TRANGEa,H2OSOI,
                      VEGC,TOTC,FPG,FPI,HR,AR,BTRAN,GROW,RAIN,SNOW,PPT)
  data$VEG  <- factor(data$VEG, levels=vegLegend)
#  data$TX   <- factor(data$TX, levels=c('TMS','Control'))
  data$TX   <- factor(data$TX, levels=c('Control'))
data$VEGno<- as.numeric(data$VEG)
  data$TXno <- as.numeric(data$TX)
  
  #--------make gg box plots of data-------------------
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  shortPalette <- c("#000000", "#D55E00")

  plottheme<-theme_classic()+theme(text = element_text(size=15),axis.title.x=element_blank(),
                                  panel.border = element_rect(colour = "black", fill=NA), 
                            #     axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                            #     axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
                                  legend.position="none")
 
  plottheme2<-theme_classic() + theme(text = element_text(size=15),axis.title.x=element_blank(),
                                      panel.border = element_blank(), axis.line = element_line())
  plot0 <-ggplot(data=data, aes(x=VEG, y=FPSN, colour=TX)) + geom_boxplot(lwd=1.0)+ plottheme + 
            ylab(paste('Photosynthesis'))
  plot1 <-ggplot(data=data, aes(x=VEG, y=GPP,  colour=TX)) + geom_boxplot(lwd=1.0)+ plottheme +
            ylab(expression(paste('GPP (gC ',m^-2,' ', y^-1,')',sep='')))  
  plot1a<-ggplot(data=data, aes(x=VEG, y=GPP,  colour=TX)) + geom_boxplot(lwd=1.0)+ plottheme +
            ylab(expression(paste('GPP (gC ',m^-2,' ', y^-1,')',sep='')))  +
            theme(axis.text.x=element_blank())
  plot1a<-plot1a +  geom_point(data=gpp_obs_data, shape=18, size=4,col=1) + ylim(0,600) 
  plot1a<-plot1a + annotate('text', x=0.6, y=(600*0.95), label='(a)',size=6) 
#+ ggtitle('(a)') + theme(plot.title=element_text(hjust=0,vjust=1))
  plot2 <-ggplot(data=data, aes(x=VEG, y=NPP,  colour=TX)) + geom_boxplot(lwd=1.0)+ plottheme + 
            ylab(expression(paste('NPP (gC ',m^-2,' ', y^-1,')',sep=''))) +
            theme(axis.text.x=element_blank())
  plot2 <-plot2 + geom_point(data=npp_obs_data, shape=18, size=4,col=1) + ylim(0,300)
  plot2 <-plot2  + annotate('text', x=0.6, y=(300*0.95), label='(b)',size=6) 

  plot2b<-ggplot(data=data, aes(x=VEG, y=ANPP,  colour=TX)) + geom_boxplot(lwd=1.0)+ plottheme + 
            ylab(expression(paste('ANPP (gC ',m^-2,' ', y^-1,')',sep=''))) #+
#            theme(axis.text.x=element_blank())
  plot2b <-plot2b + geom_point(data=anpp_obs_data, shape=18, size=4,col=1) + ylim(0,150)
  plot2b <-plot2b + annotate('text', x=0.6, y=(150*0.95), label='(c)',size=6) 
  plot3a<-ggplot(data=data, aes(x=VEG, y=NEP,  colour=TX)) + geom_boxplot(lwd=1.0)+ plottheme2+
            ylab(expression(paste('NEP (gC ',m^-2,' ', y^-1,')',sep='')))
  plot3 <-ggplot(data=data, aes(x=VEG, y=NEE,  colour=TX)) + geom_boxplot(lwd=1.0)+ plottheme2+
            ylab(expression(paste('NEE (gC ',m^-2,' ', y^-1,')',sep='')))
  plot4 <-ggplot(data=data, aes(x=VEG, y=HR,   colour=TX)) + geom_boxplot(lwd=1.0)+ plottheme +
            ylab(expression(paste('Heterotrophic Resp. (gC ',m^-2,' ', y^-1,')',sep=''))) +
            theme(axis.text.x=element_blank())
  plot5 <-ggplot(data=data, aes(x=VEG, y=AR,   colour=TX)) + geom_boxplot(lwd=1.0)+ plottheme +
            ylab(expression(paste('Autotrophic Resp. (gC ',m^-2,' ', y^-1,')',sep='')))
  plot6 <-ggplot(data=data, aes(x=VEG, y=BTRAN,colour=TX)) + geom_boxplot(lwd=1.0)+ plottheme +
            ylab('Soil moisture stress (unitless)')
  plot6 <-plot6 + ylim(0,0.8) + annotate('text', x=0.8, y=(0.8*0.98), label='(c)',size=6) 
  plot7 <-ggplot(data=data, aes(x=VEG, y=FPG,  colour=TX)) + geom_boxplot(lwd=1.0)+ plottheme +
            ylab('N limitation (unitless)')
  plot7 <-plot7 + ylim(0,0.8) + annotate('text', x=0.8, y=(0.8*0.98), label='(d)',size=6) 
  plot8 <-ggplot(data=data, aes(x=VEG, y=FPI,  colour=TX)) + geom_boxplot(lwd=1.0)+ plottheme2+
            ylab('Microbial N limitation (unitless)')
  plot9 <-ggplot(data=data, aes(x=VEG, y=TLAI, colour=TX)) + geom_boxplot(lwd=1.0)+ plottheme2+
            ylab(expression(paste('Leaf Area Index (',m^2,' ',m^-2,')',sep='')))
  plot10 <-ggplot(data=data,aes(x=VEG, y=TSOIg, colour=TX)) + geom_boxplot(lwd=1.0)+ plottheme +
            ylab(expression('Summer Soil Temp. (',~degree~,'C)'))+theme(axis.text.x=element_blank())
  plot10b<-ggplot(data=data,aes(x=VEG, y=TSOIa, colour=TX)) + geom_boxplot(lwd=1.0)+ plottheme +
            ylab(expression('Mean Soil Temp. ('*~degree*C*')'))+theme(axis.text.x=element_blank())
  plot10b<-plot10b + ylim(1.5,6) + annotate('text', x=0.8, y=(6*0.97), label='(a)',size=6) + scale_y_continuous(labels = fmt_dcimals(1))

  plot10c<-ggplot(data=data,aes(x=VEG, y=TRANGEa, colour=TX)) + geom_boxplot(lwd=1.0)+ plottheme +
            ylab('Mean Soil Temp. Range (C)')+ theme(axis.text.x=element_blank())
  plot11 <-ggplot(data=data,aes(x=VEG,y=H2OSOI,colour=TX)) + geom_boxplot(lwd=1.0)+ plottheme2+
            ylab(expression(paste('Soil Moisture (g ', g^-1,')',sep='')))
  plot12<-ggplot(data=data, aes(x=VEG, y=VEGC, colour=TX)) + geom_boxplot(lwd=1.0)+ plottheme+
            ylab(expression(paste('Plant Biomass (gC ',m^-2,')',sep='')))
  plot12 <-plot12 + geom_point(data=biomass_obs_data, shape=18, size=4,col=1)
  plot12b<-ggplot(data=data, aes(x=VEG, y=TOTC, colour=TX)) + geom_boxplot(lwd=1.0)+ plottheme+
           ylab(expression(paste('Ecosystem C (gC ',m^-2,')',sep='')))
  plot13<-ggplot(data=data, aes(x=VEG, y=GROW, colour=TX)) + geom_boxplot(lwd=1.0)+ plottheme +
            ylab('Growing season (d)')+ theme(axis.text.x=element_blank())
  plot13<-plot13 + ylim(100,200) + annotate('text', x=0.8, y=(200*0.99), label='(b)',size=6) 
  plot14<-ggplot(data=data, aes(x=VEG, y=(FPG*BTRAN), colour=TX)) + geom_boxplot(lwd=1.0)+ plottheme
  
  fout <- paste(path,'Fig4_productivity.pdf',sep='')
  pdf(fout, compress = FALSE)
  multiplot(plot1a,plot2, plot2b,cols=1)
  dev.off()

  fout <- paste(path,'Fig5_site stats.pdf',sep='')
  pdf(fout, compress = FALSE)
  multiplot(plot10b, plot6, plot13, plot7, cols=2)
  dev.off()

  fout <- paste(path,'noRHphn_ggplots_summary.pdf',sep='')
  pdf(fout, compress = FALSE)
  multiplot(plot10b, plot6, plot13, plot7, cols=2)
  dev.next()
  multiplot(plot0, plot1,cols=1)
  dev.next()
  multiplot(plot2, plot2b,cols=1)
  dev.next()
  multiplot(plot2, plot12,cols=1)
  dev.next()
  multiplot(plot4, plot5,cols=1)
  dev.next()
  multiplot(plot10b, plot10c,cols=1)
  dev.next()
# multiplot(plot3a, plot3,cols=1)
#  dev.next()
  plot8
  dev.next()
  plot9
  dev.next()
# plot11
  par(mfrow=c(1,1), mar=c(5,5,2,2), oma=c(0,0,0,0))
  lim = c(0,100)
  plot(data$ANPP, anpp_obs_data$ANPP, pch=as.numeric(data$VEG), cex=1.4, 
       ylim=lim, ylab=NA,
       xlim=lim, xlab=NA )
  abline(a=0, b=1, lty=2)
  legend('bottomright', pch=as.numeric(data$VEG), legend=data$VEG[1:5], 
         bty = "n", cex=1.4)
  mtext(expression(paste('(gC ',m^-2,' ',y^-1,')', sep='')), side=2, line=2)
  mtext(expression(paste('(gC ',m^-2,' ',y^-1,')', sep='')), side=1, line=3.5)
  mtext('Simulated ANPP', side=1, line=2, font=2, cex=1.3)
  mtext('Observed ANPP' , side=2, line=3.5, font=2, cex=1.3)
dev.off()

cor.test(data$ANPP, anpp_obs_data$ANPP)
error <- na.omit(data$ANPP - anpp_obs_data$ANPP)
n     <- length(error)
bias  <- sum(error)/n  
mean(data$ANPP)

cbind(rowMeans(GPP_a[veg,]),rowMeans(NPP_a[veg,]),rowMeans(ANPP_a[veg,]),rowMeans(BNPP_a[veg,]))
cbind(NPP_a[veg,], ANPP_a[veg,])
cbind(rowMeans(NPP_a[veg,1:5]),(rowMeans(ANPP_a[veg,])+rowMeans(BNPP_a[veg,2:6])))
cbind(rowMeans(ANPP_a[veg,]),VEG_ANPP_ALLmean)
rowMeans(ANPP_a[veg,])-VEG_ANPP_ALLmean

  #---------------------------------------------------------
  #------ Change in Biomass plots -------------------------------
  #---------------------------------------------------------
  library(scales)     # Need the scales package

  dVEGc_M_Awarm <- as.vector(VEGC_a[M_Awarm,] - VEGC_a[veg,])
  pVEGc_M_Awarm <- as.vector(log(VEGC_a[M_Awarm,]/ VEGC_a[veg,] ))
  dBTRAN_M_Awarm <- as.vector(BTRAN_a[M_Awarm,] - BTRAN_a[veg,])
  dFPG_M_Awarm  <- as.vector(FPG_a[M_Awarm,]  - FPG_a[veg,])
  dTSOI_M_Awarm <- as.vector(TSOI_a[M_Awarm,] - TSOI_a[veg,])
  dGROW_M_Awarm <- as.vector(GROW_a[M_Awarm,] - GROW_a[veg,])
  pGROW_M_Awarm <- as.vector(log(GROW_a[M_Awarm,] /  GROW_a[veg,] ))
  dGPP_M_Awarm  <- as.vector(GPP_a[M_Awarm,] -  GPP_a[veg,] )
  pGPP_M_Awarm  <- as.vector(log(GPP_a[M_Awarm,] /  GPP_a[veg,] ))
  dNPP_M_Awarm  <- as.vector(NPP_a[M_Awarm,] -  NPP_a[veg,] )
  pNPP_M_Awarm  <- as.vector(log(NPP_a[M_Awarm,] /  NPP_a[veg,] ))
  pAR_M_Awarm   <- as.vector(log(AR_a[M_Awarm,]  / AR_a[veg,] ))

  dVEGc_BlackSand <- as.vector(VEGC_a[BlackSand,] - VEGC_a[veg,])
  pVEGc_BlackSand <- as.vector(log(VEGC_a[BlackSand,]/ VEGC_a[veg,] ))
  dBTRAN_BlackSand <- as.vector(BTRAN_a[BlackSand,] - BTRAN_a[veg,])
  dFPG_BlackSand <- as.vector(FPG_a[BlackSand,] - FPG_a[veg,])
  dGROW_BlackSand <- as.vector(GROW_a[BlackSand,] - GROW_a[veg,])
  pGROW_BlackSand <- as.vector(log(GROW_a[BlackSand,] /  GROW_a[veg,] ))
  dGPP_BlackSand  <- as.vector(GPP_a[BlackSand,] -  GPP_a[veg,] )
  pGPP_BlackSand  <- as.vector(log(GPP_a[BlackSand,] /  GPP_a[veg,] ))
  dNPP_BlackSand  <- as.vector(NPP_a[BlackSand,] -  NPP_a[veg,] )
  pNPP_BlackSand  <- as.vector(log(NPP_a[BlackSand,] /  NPP_a[veg,] ))
  pAR_BlackSand   <- as.vector(log(AR_a[BlackSand,]  / AR_a[veg,] ))
  pGROW_BlackSand <- as.vector(log(GROW_a[BlackSand,] /  GROW_a[veg,] ))


#  dVEGc_ndep <- as.vector(VEGC_a[ndep,] - VEGC_a[veg,])
#  pVEGc_ndep <- as.vector(log(VEGC_a[ndep,] / VEGC_a[veg,] ))
  
  dBTRAN<- c(dBTRAN_M_Awarm,  dBTRAN_BlackSand) 
  dFPG<- c(dFPG_M_Awarm,  dFPG_BlackSand) 
  dVEGC <- c(dVEGc_M_Awarm,  dVEGc_BlackSand) 
  pVEGC <- c(pVEGc_M_Awarm,  pVEGc_BlackSand) 
  dGROW <- c(dGROW_M_Awarm,  dGROW_BlackSand) 
  pGROW <- c(pGROW_M_Awarm,  pGROW_BlackSand) 
  dGPP  <- c(dGPP_M_Awarm,  dGPP_BlackSand) 
  pGPP  <- c(pGPP_M_Awarm,  pGPP_BlackSand) 
  dNPP  <- c(dNPP_M_Awarm,  dNPP_BlackSand) 
  pNPP  <- c(pNPP_M_Awarm,  pNPP_BlackSand) 
  pAR   <- c(pAR_M_Awarm,  pAR_BlackSand) 
  TX    <- c(rep("TMS",length(dVEGc_M_Awarm)), rep("BS",length(dVEGc_BlackSand)))
  VEG   <- factor(c(rep(vegLegend,length(TX)/length(vegLegend))), as.ordered(vegLegend))
  year_out  <- YEAR[1:30]
  maxSNOW   <- as.vector(Snow_Depth_a[M_Awarm,])
  data2  <- data.frame(dVEGC, pVEGC, dNPP, pNPP,pGPP,dGPP,pAR,TX, VEG,maxSNOW,year_out,dGROW,dBTRAN,dFPG)
  write.csv(data2,paste(path,'clm_Tvan_TMS.csv', sep=''))


  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  shortPalette <- c("#000000", "#D55E00")

#  plottheme<-theme_classic() + theme(text = element_text(size=15),axis.title.x=element_blank(),legend.position="none")
#  plottheme2<-theme_classic() + theme(text = element_text(size=15),axis.title.x=element_blank())
plotdGROW <-ggplot(data=data2, aes(x=VEG, y=dGROW, colour=TX)) + 
  geom_boxplot(lwd=1.0)+ plottheme + 
  geom_hline(yintercept=0, linetype=2) + 
  theme(axis.text.x=element_blank()) +
  ylab(expression(paste(Delta,' Growing season (d)'))) +
  scale_color_manual(values=shortPalette) +
  annotate('text', x=0.8, y=34, label='(a)',size=6) 

plotdFPG <-ggplot(data=data2, aes(x=VEG, y=dFPG, colour=TX)) + 
  geom_boxplot(lwd=1.0)+ plottheme + 
  geom_hline(yintercept=0, linetype=2) + 
  ylab(expression(paste(Delta,' N limitation (unitless)'))) +
  scale_color_manual(values=shortPalette) +
  annotate('text', x=0.8, y=0.1, label='(c)',size=6) 

plotdBTRAN <-ggplot(data=data2, aes(x=VEG, y=dBTRAN, colour=TX)) + 
  geom_boxplot(lwd=1.0)+ plottheme + 
  geom_hline(yintercept=0, linetype=2) + 
  theme(axis.text.x=element_blank()) +
  ylab(expression(paste(Delta,' soil moisture stress (unitless)'))) +
  scale_color_manual(values=shortPalette) +
  annotate('text', x=0.8, y=0.53, label='(b)',size=6) 

plotBIOMASS <-ggplot(data=data2, aes(x=VEG, y=pVEGC, colour=TX)) + 
  geom_boxplot(lwd=1.0)+ plottheme + 
  geom_hline(yintercept=0, linetype=2) + 
  ylab(expression(paste('Biomass (effect size)'))) +
  scale_color_manual(values=shortPalette) +
  annotate('text', x=0.8, y=(0.12), label='(d)',size=6) 

fout <- paste(path,'Fig6_Tx_effects.pdf',sep='')
pdf(fout, compress = FALSE)
par(mar=c(5,5,2,2),oma=c(0,0,0,0))
multiplot(plotdGROW, plotdFPG,plotdBTRAN, plotBIOMASS, cols=2)
dev.off()

yearOUT <- c(year_out,year_out)
tapply(dGROW, list(TX,VEG,yearOUT), mean)
tapply(pGPP, list(TX,VEG,yearOUT), mean)

plot1 <-ggplot(data=data2, aes(x=VEG, y=pAR, colour=TX)) + 
  geom_boxplot(lwd=1.0)+ plottheme + 
  geom_hline(yintercept=0, linetype=2) + 
  ylab(expression(paste('AR (ln[Tx/control])'))) +
  theme(axis.text.x=element_blank()) +
  scale_color_manual(values=shortPalette)

plotGPP <-ggplot(data=data2, aes(x=VEG, y=pGPP, colour=TX)) + 
  geom_boxplot(lwd=1.0)+ plottheme + 
  geom_hline(yintercept=0, linetype=2) + 
  ylab(expression(paste('GPP (ln[Tx/control])'))) +
  theme(axis.text.x=element_blank()) +
  scale_color_manual(values=shortPalette)

plotdNPP <-ggplot(data=data2, aes(x=VEG, y=dNPP, colour=TX)) + 
  geom_boxplot(lwd=1.0)+ plottheme + 
  geom_hline(yintercept=0, linetype=2) + 
  ylab(expression(paste('dNPP (gC m-2 y-1)'))) +
  scale_color_manual(values=shortPalette)

plotNPP <-ggplot(data=data2, aes(x=VEG, y=pNPP, colour=TX)) + 
  geom_boxplot(lwd=1.0)+ plottheme + 
  geom_hline(yintercept=0, linetype=2) + 
  ylab(expression(paste('NPP (ln[Tx/control])'))) +
  scale_color_manual(values=shortPalette)

fout <- paste(path,'More_Tx_effects.pdf',sep='')
pdf(fout, compress = FALSE)
par(mar=c(5,5,2,2),oma=c(0,0,0,0))

multiplot(plotdGROW, plot1, plotGPP, plotNPP, cols=2)
dev.next()

plot1
  dev.next()
  
  plot(dVEGc_M_Awarm~as.vector(Snow_Depth_a[M_Awarm,]), ylab='Biomass Change (ln[Tx/control])',
       xlab='Max Snow Depth (m)')
  abline(h=0, lty=2) 

plot(dGPP~dGROW, pch=as.numeric(VEG), col=as.numeric(data2$TX))

  plot0 <-ggplot(data=data2, aes(x=VEG, y=dVEGC, colour=TX)) + 
    geom_boxplot(lwd=1.0)+ plottheme2 + 
    geom_hline(yintercept=0, linetype=2) + 
    ylab(expression(paste('Biomass Change (gC ',m^-2,')'))) +
  scale_color_manual(values=shortPalette)
  plot0

dev.off()
print(paste('wrote',fout))
  
  cbind(vegLegend,as.numeric(rowMeans(dVEGc)[1:5]))
  colMeans(dVEGc)
  cbind(rowMeans(VEGC_a[veg,]),    rowMeans(VEGC_a[M_Awarm,]),  rowMeans(VEGC_a[BlackSand,]))
  cbind(rowMeans(GPP_a[veg,]),    rowMeans(GPP_a[M_Awarm,]),  rowMeans(GPP_a[BlackSand,]))
  cbind(rowMeans(NPP_a[veg,]),    rowMeans(NPP_a[M_Awarm,]),  rowMeans(NPP_a[BlackSand,]))
  cbind(rowMeans(NEE_a[veg,]),    rowMeans(NEE_a[M_Awarm,]),  rowMeans(NEE_a[BlackSand,]))
  cbind(rowMeans(AR_a[veg,]),     rowMeans(AR_a[M_Awarm,]) ,  rowMeans(AR_a[BlackSand,]))
  cbind(rowMeans(HR_a[veg,]),     rowMeans(HR_a[M_Awarm,]) ,  rowMeans(HR_a[BlackSand,]))
  cbind(rowMeans(TSOI_a[veg,]),   rowMeans(TSOI_a[M_Awarm,]), rowMeans(TSOI_a[BlackSand,]))
  cbind(rowMeans(SOILLIQ_g[veg,]),rowMeans(SOILLIQ_g[M_Awarm,]),rowMeans(SOILLIQ_g[BlackSand,]))
  cbind(rowMeans(H2OSOI_g[veg,]), rowMeans(H2OSOI_g[M_Awarm,]), rowMeans(H2OSOI_g[BlackSand,]))
  
#---------------------------------------------------------
#------ Make GPP plots -----------------------------------
#---------------------------------------------------------
  #Model 1: Full model
  lnMAP  <- as.vector(log(PPT_a[1:5,]))
  MAP    <- as.vector(PPT_a[1:5,])
  gpp    <- as.vector(GPP_a[1:5,])
  gpp_norm1 <- glm(gpp ~ lnMAP + MAP ,family=gaussian(link="log"))
  gpp_norm1
  #Residuals look fine - scale location plot shows constant variance.
  #Nonlinearity is not apparent in the residual plot - good.
  #QQ plot indicates some problems. Maybe a little peakier than a Normal and
  #very slightly skewed left but I don't think this is an issue. See histogram. 
  #We conclude that the Normal distribution is good enough.
  par(mfrow=c(3,2), mar=c(4,4,2,1))
  plot(gpp_norm1)
  hist(gpp_norm1$residuals, xlab = "Residuals", main = "Histogram of residuals",freq=FALSE)
  rr <- seq(min(gpp_norm1$residuals),max(gpp_norm1$residuals),length.out=100)
  lines(rr,dnorm(rr,mean=0,sd=sd(gpp_norm1$residuals)),col="red")
  box()
  
  #Parameter estimates for Model 1
  gpp_norm1_coef <- matrix(coef(gpp_norm1),1,3,byrow=TRUE)
#  gpp_norm1_coef[2,] <- colSums(gpp_norm1_coef)
#  rownames(gpp_norm1_coef) <- c("Lowland","Montane")
  colnames(gpp_norm1_coef) <- c("a","b","(-)c")
  gpp_norm1_coef[1,1]<- exp(gpp_norm1_coef[1,1])
#  gpp_norm1_coef[2,1]<- exp(gpp_norm1_coef[2,1])
  gpp_norm1_coef
  
  fout <- paste(path,'Figs/warm_GPP_vs_PPT.pdf',sep='')
  pdf(fout, compress = FALSE)
  par(mfrow=c(1,1), mar=c(5,5,2,2))
  MAPseq <- seq(300,2100,20)
  plot(GPP_a[1:8,]~  PPT_a[1:8,], cex=1.5, pch=16, col=1,
    #pch=PCH[odd2,], col=COL[odd2,], 
    xlab=expression(paste('Precipitation (mm ',y^-1,')',sep='')),
    ylab=expression(paste('GPP (gC',m^-2," ",y^-1,')',sep='')),
    xlim=c(250,2600), ylim=c(0,600))
    #legend('topright', pch=pch,legend=ppt,bty='n')
  points(GPP_a[1:8,]~  PPT_a[1:8,], cex=1.5, pch=1)
  preds <- predict(gpp_norm1, data.frame(MAP=MAPseq, lnMAP=log(MAPseq)), type="response")
  preds2 <- predict(gpp_norm1, type="response")
  cor.test(preds2,gpp)
  lines(MAPseq, preds, col=1, lwd=3)

  lnMAP  <- as.vector(log(PPT_a[M_Awarm,]))
  MAP    <- as.vector(PPT_a[M_Awarm,])
  gpp    <- as.vector(GPP_a[M_Awarm,])
  gpp_norm2 <- glm(gpp ~ lnMAP + MAP ,family=gaussian(link="log"))
  gpp_norm2
  preds2 <- predict(gpp_norm2, data.frame(MAP=MAPseq, lnMAP=log(MAPseq)), type="response")
  points(GPP_a[M_Awarm,]~  PPT_a[M_Awarm,], cex=1.5, pch=16, col=2)
 # lines(MAPseq, preds2, col=2, lwd=3)  
  dev.next()
  
  boxplot.matrix(GPP_a[veg,],use.cols=F,ylab="GPP (gC/m2/y)",names=vegLegend)
  apply(GPP_a[veg,], MARGIN=1,FUN=mean,na.rm=T)
  apply(GPP_a[veg,], MARGIN=1,FUN=sd,na.rm=T)
  dev.next()
  
  GPPrr <- GPP_a[M_Awarm,] / GPP_a[veg,]  
  VEGCrr <- VEGC_a[M_Awarm,]/VEGC_a[veg,]

  par(mfrow=c(2,1), mar=c(0,5,3,2))
  boxplot.matrix(GPPrr,use.cols=F,ylab="GPP response (M_Awarm/control)",names=vegLegend)
    abline(h=1, lty=2)  
  par(mar=c(3,5,0,2))
  boxplot.matrix(VEGCrr,use.cols=F,ylab="VEGC response (M_Awarm/control)",names=vegLegend)
     abline(h=1, lty=2)

  dev.off()
  print(paste('wrote',fout))
  
  t.test(VEGC_a[veg,], VEGC_a[M_Awarm,], paired=T)
  t.test(VEGC_a[veg[1:2],], VEGC_a[M_Awarm[1:2],], paired=T)
  t.test(VEGC_a[veg[3:5],], VEGC_a[M_Awarm[3:5],], paired=T)
  t.test(GPP_a[veg,], GPP_a[M_Awarm,], paired=T)
  t.test(GPP_a[veg[1:2],], GPP_a[M_Awarm[1:2],], paired=T)
  t.test(GPP_a[veg[3:5],], GPP_a[M_Awarm[3:5],], paired=T)
  
  
#---------------------------------------------------------
#------ Make NEE plots -----------------------------------
#---------------------------------------------------------
  fout <- paste(path,'Figs/warm_NEE_vs_PPT.pdf',sep='')
  pdf(fout, compress = FALSE)
  par(mfrow=c(1,1), mar=c(5,5,2,2))
  plot(NEE_a[veg,]~  PPT_a[veg,], cex=1.5, pch=16, col=1,
       #pch=PCH[odd2,], col=COL[odd2,], 
       xlab="Precipitation (mm/y)",ylab='NEE (gC/m2/y)',
       xlim=c(250,2600), ylim=c(-150,150))
  abline(h=0,lty=2)
  dev.next()
  
  plot(NEE_a[veg,]~  PPT_a[veg,], pch=vegpch, col=1, cex=1.5,
       xlab="Precipitation (mm/y)",ylab='NEE (gC/m2/y)',
       xlim=c(200,3000), ylim=c(-200,200))
  legend('bottomright', pch=vegpch[1:5],legend=vegLegend,bty='n',cex=1.1)
  abline(h=0,lty=2)
#  points(NEE_a[(veg+1),]~  PPT_a[(veg+1),],pch=vegpch, col=2)
  dev.next()
  
  NEErr <- NEE_a[M_Awarm,] / NEE_a[veg,]  
  boxplot.matrix(NEE_a[veg,],use.cols=F,ylab="NEE (gC/m2/y)",names=vegLegend)
  boxplot.matrix(NEErr,use.cols=F,ylab="NEE response (warm/control)",names=vegLegend)
  abline(h=0, lty=2)
  dev.next()
  
  ylim=c(0.8,1.3)
  ARrr   <- AR_a[M_Awarm,] / AR_a[veg,]
  HRrr   <- HR_a[M_Awarm,] / HR_a[veg,]
  par(mfrow=c(2,1), mar=c(0,5,3,2))
  boxplot.matrix(ARrr,use.cols=F,ylab="ARrr",names=vegLegend, xaxt='n',
                 ylim=ylim, main="warm/control")
  abline(h=1, lty=2)
  par(mar=c(3,5,0,2))
  boxplot.matrix(HRrr,use.cols=F,ylab="HRrr",names=vegLegend,
                 ylim=ylim)
  abline(h=1, lty=2)
  dev.next()
  
  RECO   <- AR_a + HR_a
  RECOrr <- RECO[M_Awarm,] / RECO[veg,]
  par(mfrow=c(2,1), mar=c(0,5,3,2))
  boxplot.matrix(GPPrr,use.cols=F,ylab="GPPrr",xaxt='n', 
                 ylim=ylim, main="warm/control")
  abline(h=1, lty=2)
  par(mar=c(3,5,0,2))
  boxplot.matrix(RECOrr,use.cols=F,ylab="Reco_rr",names=vegLegend,
                 ylim=ylim)
  abline(h=1, lty=2)
  dev.off()
  

  H2OSOIrr  <- H2OSOI_g[M_Awarm,]  / H2OSOI_g[veg,]
  SOILLIQrr <- SOILLIQ_g[M_Awarm,] / SOILLIQ_g[veg,]
  rowMeans(H2OSOIrr)
  rowMeans(SOILLIQrr)
  
  #---------------------------------------------------------
  #------ Make NEE plots -----------------------------------
  #---------------------------------------------------------
  fout <- paste(path,'Figs/warm_Limitation.pdf',sep='')
  par(mfrow=c(2,1), mar=c(0,5,3,2))
  boxplot.matrix(BTRAN_a[veg,],use.cols=F,ylab="Water limitation", 
                 ylim=c(0.5,1), xaxt='n')
  par(mar=c(3,5,0,2))
  boxplot.matrix(FPG_a[veg,],use.cols=F,ylab="N limitation",names=vegLegend,
                 ylim=c(0.5,1))
  dev.off()

  
  plot(NEE_a[veg[1:5],]~SNOW_a[veg[1:5],], pch=c(1:5))
  points(NEE_a[veg[3:5]+1,]~SNOW_a[veg[3:5],], pch=c(1:3), col=2)
  plot(NEE_a[veg[3:5],]~SNOW_a[veg[3:5]+1,], pch=c(1:3))
  plot(NEE_a[1:6,]~SNOW_a[1:6,], pch=c(1:6))
  plot(NEE_a[1:6,]~RAIN_a[1:6,], pch=c(1:6))
  snowdom <- c(7:10,15,16)
  plot(NEE_a[snowdom,]~SNOW_a[snowdom,], pch=c(snowdom))
  plot(NEErr[3:4,]~SNOWrr[3:4,], pch=as.numeric(3:4)  )
  
  plot(BTRAN_a ~ PPT_a)
  boxplot.matrix(BTRAN_a[veg,],use.cols=F,ylab="median BTRAN",names=vegLegend)
  boxplot.matrix(BTRAN_a[veg+1,],use.cols=F,ylab="median BTRAN",names=vegLegend)
  boxplot.matrix(FPG_a[veg,],use.cols=F,ylab="median FPG",names=vegLegend)
  boxplot.matrix(FPG_a[veg+1,],use.cols=F,ylab="median FPG",names=vegLegend)
  
  
#------------------------------------------------------------------
#------------------------------------------------------------------
# Other analyses  
#------------------------------------------------------------------
#------------------------------------------------------------------
  
  
  plot(data$GPP,data$NEE, pch=vegpch)
  #No significant relationships w/ NEE & precip across veg communities
  cor.test(data$NEE ,data$RAIN)
  cor.test(data$NEE ,data$SNOW)
  cor.test(data$NEE ,data$PPT)
  plot(data$NEE~data$RAIN, pch=vegpch)
  plot(data$NEE~data$SNOW, pch=vegpch)
  
  NEEcorRAIN_p <- rep(NA,5) 
  NEEcorRAIN_r <- rep(NA,5) 
  NEEcorSNOW_p <- rep(NA,5) 
  NEEcorSNOW_r <- rep(NA,5) 
  
  for (i in 1:5) {
    NEEcorRAIN_p[i] <- cor.test(data$NEE[data$VEGno==i],data$RAIN[data$VEGno==i])[3]
    NEEcorRAIN_r[i] <- cor.test(data$NEE[data$VEGno==i],data$RAIN[data$VEGno==i])[4]
    
    NEEcorSNOW_p[i] <- cor.test(data$NEE[data$VEGno==i],data$SNOW[data$VEGno==i])[3]
    NEEcorSNOW_r[i] <- cor.test(data$NEE[data$VEGno==i],data$SNOW[data$VEGno==i])[4]
  }
  plot(data$NEE[data$VEGno == 4]~data$RAIN[data$VEGno == 4], col=data$TXno [data$VEGno == 4]+1)
  plot(NEE_a[15:16,]~RAIN_a[7:8,]) #actual rain received, rest is runoff
  RAIN_effects <- RAIN_a[15:16,] / RAIN_a[7:8,]
  GPP_effects <- GPP_a[15:16,] / GPP_a[7:8,]
  QRUNOFF_effects <- QRUNOFF_a[15:16,] / QRUNOFF_a[7:8,]
  BTRAN_effects <- BTRAN_a[15:16,] / BTRAN_a[7:8,]
  FPG_effects   <- FPG_a[15:16,]   / FPG_a[7:8,]
  plot(NEE_a[15:16,]~runoff) #actual rain received, rest is runoff
  plot(GPP_a[15:16,]~runoff) #actual rain received, rest is runoff
  plot(GPP_a[15:16,]~RAIN_a[15:16,]) #actual rain received, rest is runoff
  
  RAIN_a
  cor.test(data$NEE[data$VEGno > 2],data$RAIN[data$VEGno > 2])
  cor.test(data$NEE[data$VEGno > 2],data$SNOW[data$VEGno > 2])
  plot(data$NEE[data$VEGno >2 ]~data$SNOW[data$VEGno >2 ], pch=vegpch[3:5])
  plot(data$NEE[data$VEGno >2 ]~data$RAIN[data$VEGno >2 ], pch=vegpch[3:5]) #all snowbed!
  
  cor.test(data$NEE[data$VEGno <=2],data$SNOW[data$VEGno <=2])
  cor.test(data$NEE[data$VEGno <=2],data$RAIN[data$VEGno <=2])
  plot(data$NEE[data$VEGno <=2 ]~data$RAIN[data$VEGno <=2 ], pch=vegpch[1:2])
  plot(data$NEE[data$VEGno <=2 ]~data$SNOW[data$VEGno <=2 ], pch=vegpch[1:2])
  
  
  rainsites <- c(1)
  test <- NEE_a[rainsites,] - rowMeans(NEE_a)[rainsites]
  cor.test(GPP_a[rainsites,] , SNOW_a[rainsites,])
  cor.test(GPP_a[rainsites,] , RAIN_a[rainsites,])
  cor.test(NEE_a[rainsites,] , RAIN_a[rainsites,])
  cor.test(NEE_a[rainsites,] , SNOW_a[rainsites,])
  cor.test(NEE_a[rainsites,] , PPT_a[rainsites,])
  cor.test(BTRAN_a[rainsites,] , RAIN_a[rainsites,])
  cor.test(BTRAN_a[rainsites,] , SNOW_a[rainsites,])
  cor.test(BTRAN_a[rainsites,] , PPT_a[rainsites,])
  plot(GPP_a[rainsites,]~RAIN_a[rainsites,])
  plot(NEE_a[rainsites,]~RAIN_a[rainsites,])
  plot(BTRAN_a[rainsites,]~RAIN_a[rainsites,])
  plot(BTRAN_a[rainsites,]~PPT_a[rainsites,])
  rowMeans(NEE_a)
  
  snowsitesA <- c(7:10,15,16)
  snowsites <- c(9,15)
  SB        <- 13
  #snowsites <- c(7:10)
  #snowsites <- c(veg[3:4])#, (veg[3:5]+1))
  test <- NEE_a[snowsites,] - rowMeans(NEE_a)[snowsites]
  cor.test(GPP_a[snowsites,] , SNOW_a[snowsites,])
  cor.test(GPP_a[snowsites,] , RAIN_a[snowsites,])
  cor.test(NEE_a[snowsites,] , SNOW_a[snowsites,])
  cor.test(NEE_a[snowsites,] , RAIN_a[snowsites,])
  cor.test(NEE_a[snowsites,] , PPT_a[snowsites,])
  cor.test(FPG_a[snowsites,] , PPT_a[snowsites,])
  cor.test(FPG_a[snowsites,] , RAIN_a[snowsites,])
  cor.test(FPG_a[snowsites,] , SNOW_a[snowsites,])
  plot(NEE_a[snowsitesA,]~SNOW_a[snowsitesA,], col=0, xlim=c(200,1500))
  points(NEE_a[snowsites+1,]~SNOW_a[snowsites+1,], pch=as.numeric(snowsites), col=2)
  points(NEE_a[snowsites,]~SNOW_a[snowsites,], pch=as.numeric(snowsites), col=1)

  
  t.test(NEE_a[rainsites,], NEE_a[rainsites+1,], paired=T)
  t.test(NEE_a[snowsites,], NEE_a[snowsites+1,], paired=T)
  t.test(NEE_a[SB,], NEE_a[SB+1,], paired=T)
  
  t.test(TLAI_a[rainsites,], TLAI_a[rainsites+1,], paired=T)
  t.test(TLAI_a[snowsites,], TLAI_a[snowsites+1,], paired=T)
  t.test(TLAI_a[SB,], TLAI_a[SB+1,], paired=T)
  print(paste(mean(TLAI_a[rainsites,]),mean(TLAI_a[rainsites+1,])))
  print(paste(mean(TLAI_a[snowsites,]),mean(TLAI_a[snowsites+1,])))
  print(paste(mean(TLAI_a[SB,]),       mean(TLAI_a[SB+1,])))
  print(paste(mean(TLAI_a[odd,]),      mean(TLAI_a[even,])))
  
  t.test(BTRAN_a[rainsites,], BTRAN_a[rainsites+1,], paired=T)
  t.test(BTRAN_a[snowsites,], BTRAN_a[snowsites+1,], paired=T)
  t.test(BTRAN_a[SB,], BTRAN_a[SB+1,], paired=T)
  print(paste(mean(BTRAN_a[rainsites,]),mean(BTRAN_a[rainsites+1,])))
  print(paste(mean(BTRAN_a[snowsites,]),mean(BTRAN_a[snowsites+1,])))
  print(paste(mean(BTRAN_a[SB,]),       mean(BTRAN_a[SB+1,])))
  print(paste(mean(BTRAN_a[odd,]),      mean(BTRAN_a[even,])))

  t.test(FPG_a[rainsites,], FPG_a[rainsites+1,], paired=T)
  t.test(FPG_a[snowsites,], FPG_a[snowsites+1,], paired=T)
  t.test(FPG_a[SB,], FPG_a[SB+1,], paired=T)
  print(paste(mean(FPG_a[rainsites,]),mean(FPG_a[rainsites+1,])))
  print(paste(mean(FPG_a[snowsites,]),mean(FPG_a[snowsites+1,])))
  print(paste(mean(FPG_a[SB,]),       mean(FPG_a[SB+1,])))
  print(paste(mean(FPG_a[odd,]),      mean(FPG_a[even,])))

  t.test(GPP_a[rainsites,], GPP_a[rainsites+1,], paired=T)
  t.test(GPP_a[snowsites,], GPP_a[snowsites+1,], paired=T)
  t.test(GPP_a[SB,], GPP_a[SB+1,], paired=T)
  t.test(NEE_a[odd,], NEE_a[even,], paired=T)
  t.test(GPP_a[odd,], GPP_a[even,], paired=T)
  print(paste(mean(GPP_a[rainsites,]),mean(GPP_a[rainsites+1,])))
  print(paste(mean(GPP_a[snowsites,]),mean(GPP_a[snowsites+1,])))
  print(paste(mean(GPP_a[SB,]),       mean(GPP_a[SB+1,])))
  print(paste(mean(GPP_a[odd,]),      mean(GPP_a[even,])))

  plot(GPP_a[snowsitesA,]~SNOW_a[snowsitesA,], col=0, xlim=c(200,1500))
  points(GPP_a[snowsites+1,]~SNOW_a[snowsites+1,], pch=as.numeric(snowsites), col=2)
  points(GPP_a[snowsites,]~SNOW_a[snowsites,], pch=as.numeric(snowsites), col=1)
  
  
  plot(NEE_a[snowsites,]~PPT_a[snowsites,], pch=as.numeric(snowsites), xlim=c(200,1500))
  points(NEE_a[snowsites+1,]~PPT_a[snowsites+1,], pch=as.numeric(snowsites), col=2)
  rowMeans(NEE_a[c(9,10,15,16,13,14),])
  plot(NEE_a[veg[3:5],]~SNOW_a[veg[3:5],], pch=as.numeric(veg[3:5]))
  plot(GPP_a[snowsites,]~SNOW_a[snowsites,], pch=as.numeric(snowsites))
  plot(BTRAN_a[snowsites,]~SNOW_a[snowsites,], pch=as.numeric(snowsites))
    
  snowbed <- c(5:14)
  snowbed <- c(11:16)

  cor.test(GPP_a[snowbed,], SNOW_a[snowbed,])
  cor.test(GPP_a[snowbed,], grow_season[snowbed,])
  cor.test(NEE_a[snowbed,], grow_season[snowbed,])
  cor.test(SNOW_a[snowbed,],grow_season[snowbed,])
  
  plot(GPP_a[snowbed,] ~ grow_season[snowbed,])
  plot(NEE_a[snowbed,] ~ grow_season[snowbed,])

  
  plot(NEE_a[5:14,] ~ grow_season[5:14,])
  
  rowMeans(grow_season,na.rm=T)
  GROWrr <- grow_season[(veg+1),] / grow_season[veg,] 
  GROWlength <- grow_season[(veg+1),] - grow_season[veg,] 
  boxplot.matrix(GROWrr,use.cols=F,ylab="Growing Season response (dry/control)",names=vegLegend)
  abline(h=1, lty=2)
  boxplot.matrix(LEAF_f[veg,],use.cols=F,ylab="Green up ",names=vegLegend)
  boxplot.matrix(LEAF_l[veg,],use.cols=F,ylab="Green up ",names=vegLegend)
  boxplot.matrix(grow_season[veg,],use.cols=F,ylab="Growing Season (d) ",names=vegLegend)
  
  apply(LEAF_f, MARGIN=1,FUN=mean,na.rm=T)
  apply(LEAF_f, MARGIN=1,FUN=sd,na.rm=T)
  apply(grow_season[odd,], MARGIN=1,FUN=mean,na.rm=T)
  apply(grow_season[odd,], MARGIN=1,FUN=sd,na.rm=T)
  rowMeans(GROWrr,na.rm=T)
  rowMeans(GROWlength,na.rm=T)
  plot(grow_season~SNOW_a)
  plot(NEE_a~  RAIN_a, pch=PCH, col=COL)
  plot(NEE_a[odd,]~  SNOW_a[odd,], pch=PCH[odd,], col=COL[odd,])
  plot(NEE_a[odd,]~  RAIN_a[odd,], pch=PCH[odd,], col=COL[odd,])
  GPP_rr <- GPP_a[even,2:6]-GPP_a[odd,2:6]
  NEE_rr <- NEE_a[even,2:6]-NEE_a[odd,2:6]
  boxplot.matrix(GPP_rr[2:7,], use.cols=F,ylab="GPP response (tx-control)")
  abline(h=0, lty=2)
  boxplot.matrix(NEE_rr[2:7,], use.cols=F, ylab="NEE response (tx-control)")
  abline(h=0, lty=2)
  rowMeans(GPP_rr)
  rowMeans(NEE_rr)
  
  yr <- seq(2008,2013,1)
  TLAIrr <- TLAI_a[even,]/TLAI_a[odd,]
  plot(TLAIrr[1,]~yr, col=0, ylim=c(0.70,1.05))
  for (i in 1:length(odd)) {
    lines(TLAIrr[i,]~yr, col=i,lwd=2)
  }
  
  GPPrr <-  GPP_a[even,] /GPP_a[odd,]
  plot(GPPrr[1,]~yr, col=0, ylim=c(0.70,1.25))
  for (i in 1:length(odd)) {
    lines(GPPrr[i,]~yr, col=i)
  }
  abline(h=1,lty=2)

  plot(GPP_a[1,]~yr, col=0, ylim=c(200,900))
  for (i in 1:length(odd)) {
    lines(GPP_a[odd[i],]~yr, col=i,lwd=2)
    abline(h=mean(GPP_a[odd[i],]), col=i, lty=2)
  }
  rowMeans(GPP_a[odd,])  
  plot(TLAIrr,GPPrr)
  
  NEErr <-  NEE_a[even,] - NEE_a[odd,]
  plot(NEErr[1,]~yr, col=0, ylim=c(-80,80))
  for (i in 1:length(odd)) {
    lines(NEErr[i,]~yr, col=i)
  }
  abline(h=0,lty=2)

  plot(NEE_a[1,]~yr, col=0, ylim=c(-150,150))
  for (i in 1:length(odd)) {
    lines(NEE_a[odd[i],]~yr, col=i,lwd=2)
    abline(h=mean(NEE_a[odd[i],]), col=i, lty=2)
    
  }

  rowMeans(NEE_a[odd,])  
  FPGrr <-  FPG_a[even,] / FPG_a[odd,]
  plot(FPGrr[1,]~yr, col=0, ylim=c(0.75,1.15))
  for (i in 1:length(odd)) {
    lines(FPGrr[i,]~yr, col=i,lwd=2)
  }
  abline(h=1,lty=2)

  plot(TLAIrr,GPPrr)
  grow_season_vr    <- rep(NA,nrows)
  NEE_vr            <- rep(NA,nrows)
  GPP_vr            <- rep(NA,nrows)
  for (i in 1:nrows) {
    NEE_vr[i]         <- var(NEE_a[i,],       na.rm=T)
    GPP_vr[i]         <- var(GPP_a[i,],       na.rm=T)
    grow_season_vr[i] <- var(grow_season[i,], na.rm=T)
  }
  var(GPP_a[10,])
  rowMeans(GPP_a)
  rowMeans(TLAI_a)
  rowMeans(LEAF_f)
  rowMeans(LEAF_l)
  rowMeans(NPP_a)
  rowMeans(NEE_a)
  colMeans(NEE_a)
  colMeans(PPT_a)
  rowMeans(SnowDP_a)
  rowMeans(BTRAN_a)
  rowMeans(FPG_a)
  
  plot(rowMeans(NEE_a)[even]-rowMeans(NEE_a)[odd])
  plot(  colMeans(NEE_a[odd,]) ~colMeans(SNOW_a[odd,]),ylim=c(-100,80),xlim=c(300,900))
  points(colMeans(NEE_a[even,])~colMeans(SNOW_a[even,]),col=2)
  
  
  
  plot(NEE_a~PPT_a)
  plot(NEE_a~RAIN_a)
  plot(NEE_a~SNOW_a)

  plot(GPP_a~PPT_a)
  plot(GPP_a~RAIN_a)
  plot(GPP_a~SNOW_a)
  
  relRAIN
  cor.test(NEE_a[relRAIN>0.75],RAIN_a[relRAIN>0.75])
  cor.test(NEE_a[relRAIN>0.75],PPT_a[relRAIN>0.75])
  cor.test(GPP_a[relRAIN>0.75],RAIN_a[relRAIN>0.75])
  cor.test(NEE_a[1,],PPT_a[1,])
  cor.test(NEE_a[1,],RAIN_a[1,])
  cor.test(NEE_a[1,],SNOW_a[1,])
  cor.test(GPP_a[1,],RAIN_a[1,])
  cor.test(GPP_a[1,],SNOW_a[1,])
  cor.test(GPP_a[1,],PPT_a[1,])
  
  plot(NEE_a[1,]~PPT_a[1,], ylim=c(-150,150))
  plot(NEE_a[1,]~RAIN_a[1,], ylim=c(-150,150))
  points(NEE_a[2,]~RAIN_a[2,], pch=16)
  points(NEE_a[3,]~RAIN_a[3,], pch=2)
  points(NEE_a[4,]~RAIN_a[4,], pch=17)
#  points(NEE_a[5,]~RAIN_a[5,], pch=17,col=2)
  
  plot(GPP_a[1,]~PPT_a[1,], ylim=c(200,800))
  plot(GPP_a[1,]~RAIN_a[1,], ylim=c(200,800))
  points(GPP_a[2,]~RAIN_a[2,], pch=16)
  points(GPP_a[3,]~RAIN_a[3,], pch=2)
  points(GPP_a[4,]~RAIN_a[4,], pch=17)
  points(GPP_a[5,]~RAIN_a[5,], pch=2, col=2)
  points(GPP_a[6,]~RAIN_a[6,], pch=17, col=2)
  
  plot(GPP_a[1:5,]~PPT_a[1:5,], ylim=c(200,1000))
  plot(GPP_a[1:8,]~PPT_a[1:8,], ylim=c(200,1000),xlim=c(0,2000))
  points(GPP_a[7:12,]~PPT_a[7:12,],pch=16)
  plot(GPP_a~PPT_a)
  
  cor.test(NEE_a[5:10,],RAIN_a[5:10,])
  cor.test(NEE_a[5:10,],SNOW_a[5:10,])
  cor.test(NEE_a[5,],PPT_a[5,])
  plot(NEE_a[5,]~PPT_a[5,], ylim=c(-150,150))
  plot(NEE_a[5,]~SNOW_a[5,], ylim=c(-200,200),xlim=c(300,1500))
  points(NEE_a[6,]~SNOW_a[6,], pch=16)
  points(NEE_a[9,]~SNOW_a[9,], pch=3)
  points(NEE_a[10,]~SNOW_a[10,], pch=18)

  points(NEE_a[7,]~SNOW_a[7,], pch=2, col=2)
  points(NEE_a[8,]~SNOW_a[8,], pch=17,col=2)
  
  cor.test(GPP_a[5:10,],SNOW_a[5:10,])
  plot(GPP_a[5,]~SNOW_a[5,], ylim=c(300,900),xlim=c(300,1500))
  points(GPP_a[6,]~SNOW_a[6,], pch=16)
  points(GPP_a[9,]~SNOW_a[9,], pch=3)
  points(GPP_a[10,]~SNOW_a[10,], pch=18)
  points(GPP_a[7,]~SNOW_a[7,], pch=2,col=2)
  points(GPP_a[8,]~SNOW_a[8,], pch=2,col=2)
  
  
  
  diff1 <- GPP_a[1,] - rowMeans(GPP_a)[[1]] 
  diff3 <- GPP_a[3,] - rowMeans(GPP_a)[[3]] 
  
  veg        <- c('FF','DM','MM','WM','SB')
  dims       <- c(length(veg),nyears)
  GPP_rr_dry <- array(NA,dim=dims, dimnames=list(veg,cols)) 
  NPP_rr_dry <- array(NA,dim=dims, dimnames=list(veg,cols)) 
  reps <- c(1,3,5,7,9)
  
  for (i in 1:length(veg)) {
    r <- reps[i]
    GPP_rr_dry[i,] <- GPP_a[(r+1),]/GPP_a[r,]
  }
  
  rr25dry <- log(GPP_a[2,] /GPP_a[1,])
  rr25wet <- log(GPP_a[3,] /GPP_a[1,])
  
  rr50dry <- log(GPP_a[5,] /GPP_a[4,])
  rr50wet <- log(GPP_a[6,] /GPP_a[4,])
  
  rr75dry <- log(GPP_a[8,] /GPP_a[7,])
  rr75wet <- log(GPP_a[9,] /GPP_a[7,])
  
  rr100dry <- log(GPP_a[11,] /GPP_a[10,])
  rr100wet <- log(GPP_a[12,] /GPP_a[10,])

  
  rr75wet[4]  <- NA
  rr100dry[4] <- NA
  rr100wet[4] <- NA
  
  par(mfrow=c(1,1),mar=c(5,4,4,2))
  boxplot(data.frame(rr25dry,rr50dry,rr75dry,rr100dry), main="reduced snow", ylab="GPP response")
  abline(h=0, lty=2)  
  
  boxplot(data.frame(rr25wet,rr50wet,rr75wet,rr100wet), main="increased snow", ylab="GPP response")
  abline(h=0, lty=2)  

  pdf("mean_CLM_snow_depth.pdf", compress = FALSE)
  boxplot.matrix(SnowDP_a[c(1,4,7,10),1:5]*100,use.cols=F,
                 main='CLM 2008-2013', ylab='max snow depth (cm)')
  dev.off()    
  
  
  
  
  cor(diff1, RAIN_a[1,])
  cor(diff3, RAIN_a[3,])
  cor(diff3, SNOW_a[3,])
  pch = as.numeric(factor(row.names(GPP_a)))
  plot(GPP_a~RAIN_a, pch=pch, ylim=c(min(GPP_a)*0.8,max(GPP_a)*1.1))
  YEARS <- c("'08","'09","'10","'11","'12","'13")
  text(RAIN_a[1,],y=min(GPP_a)*0.9,as.character(YEARS))

  plot(GPP_a~SNOW_a, pch=pch, ylim=c(min(GPP_a)*0.8,max(GPP_a)*1.1))
  plot(GPP_a~PPT_a, pch=pch, ylim=c(min(GPP_a)*0.8,max(GPP_a)*1.1))
  plot(BTRAN_a,FPG_a, pch=pch)
  
  
  
  fout <- paste('Figs_Tvan',case[e],'/Tvan',case[e],'_QRUNOFF.csv', sep='')
  
  dout <- data.frame(MCDATE,PRECIP,QRUNOFF)
  write.csv(dout,file=fout)
  
  par(mfrow=c(2,1), mar=c(2,5,1,1))
  plot(QRUNOFF, ylim=c(0,0.002), type="l")
  plot(BTRAN,  type="l")

  plot(FPG,  type="l")
  plot((FPG*BTRAN),  type="l")

  plot(FPG,  type="l",xlim=c(8000,9300))
  plot(BTRAN,  type="l",xlim=c(8000,9300))
  plot(PRECIP,  type="l",xlim=c(8000,9300))
  GPPunits <- att.get.ncdf(Data.clm, "GPP", "units")
  GPPunits
  SNOW_DEPTHunits <- att.get.ncdf(Data.clm, "SNOW_DEPTH", "units")
  SNOW_DEPTHunits


  
  