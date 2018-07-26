# Plots sample niwot data
# here data from node #13
# downloaded from SD card Jan 19, 2017
# manually deleted data w/ bad time stamps
# manuall added header to each column for ease of plotting

remove(list=ls())
library(zoo)


dir <- paste("/Users/wwieder/Desktop/Working_files/Niwot/Sensors/Sensor_data/Sadd_temp_long/",sep="")
fin <- paste(dir,"sdlcr23x-cr1000.daily.ml.data.csv",sep="") 
dat <- read.csv(fin, header = TRUE,sep=',')
data <- read.csv(paste(dir,'temp.csv',sep=''),header=T)
MART <- read.table(paste(dir,'Barnes_13.txt',sep=''), header=T, sep=',')
MART_ww <- read.csv(paste(dir,'Mart_chem/MART_ann_flux.csv',sep='')) 
MART_LYS<- read.csv(paste(dir,'William_lys_chem/MART_lys_ann.csv',sep=''))
SAD_LYS <- read.csv(paste(dir,'William_lys_chem/SAD_lys_ann.csv' ,sep=''))
#NO3 & DON = g  FLOW = m^3 / day

names(dat)
names(MART)
# n-day running mean
n = 7
temp_ma <- rollapply(dat$soiltemp_avg, n , mean, na.rm=TRUE, align="center")
temp_ma[1:5]
dat$temp_ma <- rep(NA, length(dat$jday))
dat$temp_ma[3:(length(dat$jday)-4)] <- temp_ma
remove(temp_ma)

attach(dat)

estimate_mode <- function(x) {
  d <- density(x, na.rm=T)
  d$x[which.max(d$y)]
}

par(mfrow=c(1,1))
plot(jday[year==2012],soiltemp_avg[year==2012], xlim=c(0,365), type='l')
lines(jday[year==2012],temp_ma[year==2012], col=2)
abline(h=0, lty=2)

# vecors to store results
nyears <- max(year)-min(year)+1
fday   <- rep(NA, nyears)
lday   <- rep(NA, nyears)
nday   <- rep(NA, nyears)

for(y in (min(year)+1):max(year)) {
  i     <- y-min(year)+1
  x     <- jday[year==y]
  temp  <- round(soiltemp_avg[year==y],digits=1)
  temp2 <- round(temp_ma[year==y], digits=1)
  plot(x,temp, xlim=c(0,365), main=y,type='l')
  lines(x,temp2, col=2)
  
  abline(h=0, lty=2)
  #  lines(jday[year==y],soiltemp_avg[year==y])
  min <- round(estimate_mode(temp2[temp2>-1 & temp2<1]), digits=1)
  dt  <- temp2 - min
  tol <- 0.2
    
  fday[i] <- min(x[abs(dt)<=tol], na.rm=T)
  lday[i] <- max(x[abs(dt)<=tol & x < 200], na.rm=T)
  nday[i] <- lday[i]-fday[i]
  abline(h=min,lty=2, col=2)
  abline(v=fday[i])
  abline(v=lday[i])
}

y     <- seq(min(year),max(year),1)
labels = c('00','01','02','03','04','05','06','07','08','09',
           '10','11','12','13','14')
dout  <- data.frame(y,fday,nday,labels) 

fout <- paste(dir,"NO3_vs_Thaw.pdf", sep='')
pdf(fout, width=7, height=4)
par(mar=c(4,5,2,2))

plot(fday, nday, ylab='number of isothermic days', xlab='First soil thaw day, Saddle',
     type='n')
text(fday,nday,labels)
m0 <- lm(nday~fday)
abline(m0)
p0 <- signif(summary(m0)[[4]][8],1)
legend('topright',c(paste('p =',p0)),bty='n')
dev.next()

#---------------------------------
# Mark's soil lysemeter data
#---------------------------------
SAD_LYS$y             <- SAD_LYS$X
SAD_LYS$SAD_NO3_surf  <- rowMeans(SAD_LYS[,6:9], na.rm=T)
SAD_LYS$SAD_NO3_deep  <- rowMeans(SAD_LYS[,2:5], na.rm=T)
SAD_LYS$SAD_NO3_all   <- rowMeans(SAD_LYS[,2:9], na.rm=T)
SAD_LYS$MART_NO3_surf <- rowMeans(MART_LYS[,5:7], na.rm=T)
SAD_LYS$MART_NO3_deep <- rowMeans(MART_LYS[,2:4], na.rm=T)
SAD_LYS$MART_NO3_all <- rowMeans(MART_LYS[,2:7], na.rm=T)
mer   <- merge(x = dout, y = SAD_LYS, by = "y", all = TRUE)

plot(mer$fday, mer$SAD_NO3_all, ylab='[NO3]', 
#     ylim=c(250,750),
     xlab='First soil thaw day, Saddle',
     type='n', main="Saddle soil lysemeters")
text(mer$fday,mer$SAD_NO3_all,labels)
#text(mer$fday,mer$SAD_NO3_surf,labels)
#text(mer$fday,mer$SAD_NO3_deep,labels,col=2)
#text(mer$fday,mer$MART_NO3_all,labels,col=4)

m0 <- lm(mer$SAD_NO3_all~mer$fday)
m1 <- lm(mer$SAD_NO3_surf~mer$fday)
m2 <- lm(mer$SAD_NO3_deep~mer$fday)
m3 <- lm(mer$MART_NO3_all~mer$fday)
abline(m0)
#abline(m1,col=2)
#abline(m2,col=3)
#abline(m3,col=4)
p0 <- signif(summary(m0)[[4]][8],2)
p1 <- signif(summary(m1)[[4]][8],1)
p2 <- signif(summary(m2)[[4]][8],2)
p3 <- signif(summary(m3)[[4]][8],2)
legend('topright',c(paste('all lysemeters, p =',p0)),bty='n')
#legend('topright',c(paste('shallow, p =',p1), paste('    deep, p =',p2)),
#       text.col=c(1,2),bty='n')
summary(m0)
summary(m3)

dev.next()

#---------------------------------
# look at MART data I processed
#---------------------------------
plot(fday,MART_ww$NO3_vwm*1e3,type='n',
     main='MART NO3 export, flow weighted concentration/annual discharge',
     ylab="NO3_vwm",xlab='First soil thaw day, Saddle')
text(fday,MART_ww$NO3_vwm*1e3,labels)
m0 <- lm(MART_ww$NO3_vwm*1e3~fday)
abline(m0)
summary(m0)
p0 <- signif(summary(m0)[[4]][8],1)
legend('topright',c(paste('p =',p0)),bty='n')

dev.next()

plot(fday,MART_ww$NO3_vwm_1,type='n',
     main='raw VWM_MART NO3 export, flow weighted concentration',
     ylab="NO3_vwm_1 (ug/L)",xlab='First soil thaw day, Saddle')
text(fday,MART_ww$NO3_vwm_1,labels)
m0 <- lm(MART_ww$NO3_vwm_1~fday)
abline(m0)
summary(m0)
p0 <- signif(summary(m0)[[4]][8],1)
legend('topright',c(paste('p =',p0)),bty='n')
dev.next()

plot(NO3_vwm_1~cumFLOW, data=MART_ww, type='n',
     main='raw VWM_MART NO3 export vs. annual flow',
     ylab="NO3_vwm_1 (ug/L)",xlab='Annual flow')
text(NO3_vwm_1~cumFLOW,labels, data=MART_ww)
m0 <- lm(NO3_vwm_1~cumFLOW, data=MART_ww)
abline(m0)
summary(m0)
p0 <- signif(summary(m0)[[4]][8],1)
legend('topleft',c(paste('p =',p0)),bty='n')
dev.next()

plot(fday,MART_ww$cumFLOW,type='n',
     main='Mart annual Flow',
     ylab="Annual Flow",xlab='First soil thaw day, Saddle')
text(fday,MART_ww$cumFLOW,labels)
m0 <- lm(MART_ww$cumFLOW~fday)
abline(m0)
summary(m0)
p0 <- signif(summary(m0)[[4]][8],1)
legend('topleft',c(paste('p =',p0)),bty='n')
dev.next()


plot(fday,MART_ww$cumNO3/MART_ww$cumFLOW,type='n',
     main='NO3 yield / annual Flow',
     ylab="NO3 yield/ annual Flow",xlab='First soil thaw day, Saddle')
text(fday,MART_ww$cumNO3/MART_ww$cumFLOW,labels)
m0 <- lm(MART_ww$cumNO3/MART_ww$cumFLOW~fday)
abline(m0)
summary(m0)
p0 <- signif(summary(m0)[[4]][8],1)
legend('topright',c(paste('p =',p0)),bty='n')
dev.next()


#---------------------------------
# Lake data, here from GL4
# maybe try for lower lakes too?
#---------------------------------
data$Depth
plot(fday,data$AvgNO3[data$Depth == '0'],type='n',
     main='GL4 average NO3',
     ylab="NO3",xlab='First soil thaw day, Saddle')
text(fday,data$AvgNO3[data$Depth == '0'],labels)
m0 <- lm(data$AvgNO3[data$Depth == '0']~fday)
abline(m0)
summary(m0)
p0 <- signif(summary(m0)[[4]][8],1)
legend('topright',c(paste('p =',p0)),bty='n')

dev.off()




# ----- some other plots --------------
# --- MART stream data ---
plot(MART_ww$cumFLOW,MART_ww$cumNO3,type='n', ylab="MART NO3",xlab='Mart_flow')
text(MART_ww$cumFLOW,MART_ww$cumNO3,labels)

plot(MART_ww$cumDON,MART_ww$cumNO3,type='n', ylab="MART NO3",xlab='MART DON')
text(MART_ww$cumDON,MART_ww$cumNO3,labels)

plot(fday,MART_ww$cumNO3,type='n', ylab="MART NO3",xlab='First soil thaw day, Saddle')
text(fday,MART_ww$cumNO3,labels)

plot(fday,MART_ww$cumNO3/(MART_ww$cumFLOW*1e3),type='n', 
     ylab="MART NO3/MART FLOW",xlab='First soil thaw day, Saddle')
text(fday,MART_ww$cumNO3/(MART_ww$cumFLOW*1e3),labels)
m2 <- lm(MART_ww$cumNO3/(MART_ww$cumFLOW*1e3)~fday)
abline(m2)
summary(m2)
plot(MART_ww$cumNO3/(MART_ww$cumFLOW*1e3*62),MART_NO3, type='n')
text(MART_ww$cumNO3/(MART_ww$cumFLOW*1e3*62),MART_NO3,labels)


# plot with flow weighted concentrations from Barnes et al 2013
MART_NO3 <- rep(NA, length(labels))
nMART    <- length(MART$year)
MART_NO3[1:nMART] <- MART$MRT_NO3.uM.

plot(fday,MART_NO3,type='n', ylab="MART VWM NO3 (uM)",xlab='First soil thaw day, Saddle')
text(fday,MART_NO3,labels)
m1 <- lm(MART_NO3~fday)
abline(m1)

summary(m1)
hist(MART_NO3)
hist(fday)

#---------------------------------
# More of Mark's soil lysemeter data
#---------------------------------
plot(mer$nday, mer$SAD_NO3_deep, ylab='[NO3]', 
     ylim=c(250,750),
     xlab='First soil thaw day, Saddle',
     type='n', main="Saddle, deep lysemeters")
text(mer$nday,mer$SAD_NO3_deep,labels)
text(mer$nday,mer$MART_NO3_deep,labels,col=2)

m0 <- lm(mer$SAD_NO3_deep~mer$nday)
m1 <- lm(mer$MART_NO3_deep~mer$nday)
abline(m0)
abline(m1,col=2)
summary(m0)
summary(m1)


plot(mer$nday, mer$SAD_NO3_deep, ylab='[NO3]', xlab='First soil thaw day, Saddle',
     type='n', main="Saddle, deep lysemeters")
text(mer$nday,mer$SAD_NO3_deep,labels)
m0 <- lm(mer$SAD_NO3_deep~mer$nday)
abline(m0)
summary(m0)


plot(mer$fday, mer$MART_NO3, ylab='Mart soil lys NO3', xlab='First soil thaw day, Saddle',
     type='n')
text(mer$fday,mer$MART_NO3,labels)



plot(fday,lday)
hist(fday)
hist(lday)
doy <- seq(0,365,1)
mean_temp <- tapply(soiltemp_avg,jday, mean)
plot(doy, mean_temp)
abline(h=0, lty=2)
 max(year)
 