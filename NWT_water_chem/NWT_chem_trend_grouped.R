# Will Wieder
# Created April 2015
# Modified May & Oct 2018
# Plots DOC-NO3 relationship (See Taylor and Townsend, Nature 2010)
# Looks for trends in pH, SO4.., and DOC from Nel Ca..ine's Niwot water chemistry data 

# uses gamm function [generalized additive mixed modelling] 
# in mgcv (& nlme) package

# http://stackoverflow.com/questions/12623027/how-to-analyse-irregular-time-series-in-r
# using an additive model to "decompose" the seasonal and trend components. 
# As this is a regression-based approach you need to model the residuals as a time series process to
# account for lack of independence in the residuals.
#http://www.fromthebottomoftheheap.net/2011/06/12/additive-modelling-and-the-hadcrut3v-global-mean-temperature-series/

# !! code here is really messy !! 
# It could like be made more efficient & user friendly
# switching to tidyR and ggplot may be helpful...

remove(list=ls())

library(mgcv)
library(nlme)
library(dae)
#also explore rucm and ggplot?
library(rucm)
library(ggplot2)

#load custom functions
# from https://github.com/gavinsimpson/random_code/blob/master/tsDiagGamm.R

## For GAMM models from mgcv:::gamm

## Model Checking function

tsDiagGamm <- function(x, timevar, observed, f = 0.3, type = "normalized") {
    resi <- resid(x$lme, type = type)
    fits <- fitted(x$lme)
    on.exit(layout(1))
    layout(matrix(1:6, ncol = 3, byrow = TRUE))

    plot(resi ~ fits, ylab = "Normalized Residuals",
         xlab = "Fitted Values", main = "Fitted vs. Residuals")
    lines(lowess(x = fits, y = resi, f = f), col = "blue",
          lwd = 2)
    plot(resi ~ timevar, ylab = "Normalized Residuals",
         xlab = "Time", main = "Time series of residuals")
    lines(lowess(x = timevar, y = resi, f = f), col = "blue", lwd = 2)
    plot(observed ~ fits, ylab = "Observed",
         xlab = "Fitted Values", main = "Fitted vs. Observed",
         type = "n")
    abline(a = 0, b = 1, col = "red")
    points(observed ~ fits)
    lines(lowess(x = fits, y = observed, f = f), col = "blue",
          lwd = 2)
    hist(resi, freq = FALSE, xlab = "Normalized Residuals")

    qqnorm(resi)
    qqline(resi)
    acf(resi, main = "ACF of Residuals")

}

#---------------read in DATA----------------------
# REMOVED BEAVER POND, GL4_Waterfall, GL4_Tunnel, GL4_INLET, & TUNNEL_GROUNDWATER
# eventually will need to point to EDI master data and join sites
#-------------------------------------------------
Data.in <- read.csv('all_nc_chem.csv')
names(Data.in)

min(Data.in$DOC, na.rm =TRUE)
min(Data.in$DON, na.rm =TRUE)
Data.in$date[Data.in$local_site =='MART_PRECIP']
length(Data.in$date)
levels(Data.in$local_site)

OMIT <- list('FLUME', 'GREEN LAKE 4 WATERFALL')
data.in <- Data.in
for (i in 1:length(OMIT)) {
  data.in <- data.in[ which(data.in$local_site != OMIT[i]), ]
}
remove(Data.in)

data.in$date <- as.Date(as.character(data.in$date), "%Y-%m-%d")
data.in$date 
levels(data.in$local_site)

sites  <- levels(data.in$local_site)
nsites <- length(sites)  
# remove some absurd values

data.in$DON [data.in$DON <=0] = NA
data.in$NO3.[data.in$NO3.<=0] = NA
data.in$NO3.[data.in$NO3.>=200] = 200


# SET UP CATEGORICAL VARIABLES FOR GROUPING ANALYSES
data.in$ele <- rep(NA,length(data.in$local_site))
data.in$ele[data.in$local_site == "ALBION"] <- "Low"
data.in$ele[data.in$local_site == "ALBION INLET"] <- "Low"
data.in$ele[data.in$local_site == "ARIKAREE"]     <- "High"
data.in$ele[data.in$local_site == "GREEN LAKE 1"] <- "Low"
data.in$ele[data.in$local_site == "GREEN LAKE 4"] <- "High"
data.in$ele[data.in$local_site == "GREEN LAKE 5"] <- "High"
data.in$ele[data.in$local_site == "GREEN LAKE 5 ROCK GLACIER"]       <- "High"
data.in$ele[data.in$local_site == "MARTINELLI"]   <- "Low"
data.in$ele[data.in$local_site == "NAVAJO"]       <- "High"
data.in$ele[data.in$local_site == "SADDLE STREAM"]<- "Low"

data.in$ele2 <- data.in$ele 
data.in$ele2[data.in$local_site == "GREEN LAKE 5"]    <- "High"
data.in$ele2[data.in$local_site == "GREEN LAKE 5 ROCK GLACIER"]   <- "Snow"
data.in$ele2[data.in$local_site == "NAVAJO"]   <- "Snow"
data.in$ele2[data.in$local_site == "ARIKAREE"] <- "Snow"
data.in$ele  <- as.factor(data.in$ele)
data.in$ele2 <- as.factor(data.in$ele2)

#to reorder, and claculate low elevation sites first
data.in$ELE <- rep(NA,length(data.in$local_site))
data.in$ELE[data.in$ele == "Low"]  <- "Bio"
data.in$ELE[data.in$ele == "High"] <- "Geo"
data.in$ELE <- as.factor(data.in$ELE)

levels(data.in$ele2)
levels(data.in$type)

#convert date to day of year and total days
Diff <- function(x, start) as.numeric(x - as.Date(cut(start, "year")))
DAY1 <- as.Date("1994-01-01")
data.in <- transform(data.in, DOY = Diff(date, date), TotalDays = Diff(date, DAY1))
data.in$decY <- 1994+(data.in$TotalDays/365)
allsteps     <- length(data.in$DOY)

#remove duplicate data points
data.in$temp <- paste(data.in$local_site, data.in$date)
data.in      <- data.in[!duplicated(data.in$temp), ]
dim(data.in)
#-------------Subset data, remove NA------------------------------
AllpH  <- data.frame(pH =data.in$pH,  DOY=data.in$DOY, decY=data.in$decY, year=data.in$year,
                     local_site=data.in$local_site, ele=data.in$ele, ele2=data.in$ele2)
AllSO4.. <- data.frame(SO4..=data.in$SO4.., DOY=data.in$DOY, decY=data.in$decY, year=data.in$year,
                     local_site=data.in$local_site, ele=data.in$ele, ele2=data.in$ele2)
AllDOC <- data.frame(DOC=data.in$DOC, DOY=data.in$DOY, decY=data.in$decY, year=data.in$year,
                     local_site=data.in$local_site, ele=data.in$ele, ele2=data.in$ele2)
AllNH4. <- data.frame(NH4.=data.in$NH4., DOY=data.in$DOY, decY=data.in$decY, year=data.in$year,
                     local_site=data.in$local_site, ele=data.in$ele, ele2=data.in$ele2)
AllNO3. <- data.frame(NO3.=data.in$NO3., DOY=data.in$DOY, decY=data.in$decY, year=data.in$year,
                     local_site=data.in$local_site, ele=data.in$ele, ele2=data.in$ele2)
AllPO4... <- data.frame(PO4...=data.in$PO4..., DOY=data.in$DOY, decY=data.in$decY, year=data.in$year,
                     local_site=data.in$local_site, ele=data.in$ele, ele2=data.in$ele2)
AllCa..  <- data.frame(Ca..=data.in$Ca.., DOY=data.in$DOY, decY=data.in$decY, year=data.in$year,
                     local_site=data.in$local_site, ele=data.in$ele, ele2=data.in$ele2)
AllDON <- data.frame(DON=data.in$DON, DOY=data.in$DOY, decY=data.in$decY, year=data.in$year,
                     local_site=data.in$local_site, ele=data.in$ele, ele2=data.in$ele2)
AllpH    <- na.exclude(AllpH)
AllSO4.. <- na.exclude(AllSO4..)
AllDOC   <- na.exclude(AllDOC)
AllNH4.  <- na.exclude(AllNH4.)
AllNO3.  <- na.exclude(AllNO3.)
AllDON   <- na.exclude(AllDON)
AllPO4...<- na.exclude(AllPO4...)
AllCa..  <- na.exclude(AllCa.. )

#------------------------------------------------------
# DOC-NO3. relationships
#-------------------------------------------------

data.in$NO3._conc <- data.in$NO3. * 14 / 1000 # convert to mgN/L
data.in$DOC_molar <- data.in$DOC * 1e3 / 12 # umol/L
fout   <- paste('DOC_NO3.pdf', sep='')
pdf(fout, width=4, height=4)
par(mfrow=c(1,1), mar=c(4,4,1,1), oma=c(0,0,0,0),mai=c(0.9,0.9,0.2,0.2))

col <- as.numeric(data.in$ele2)
col[col == 3] <- 4
plot(NO3.~DOC_molar, data=data.in, col=col,
     ylim=c(0,150), #xlim=c(0,500), #omits a few points
     pch=16,
#    pch=as.numeric(droplevels(data.in$local_site)),
     ylab=expression(paste('[NO'['3']*'] (',mu, "mol N/L)")),
     xlab=expression(paste("[DOC] (",mu, "mol C/L)")),
     cex=0.6)
labs <- c('glacier', 'high', 'low')
legend('topright',labs, bty='n', cex=1.2,
       pch=16, col=c(4,1,2))
text(100,150, paste('n = ', length(data.in$local_site)))
abline(0,(1/4), lty=2)
text(620*0.65, 140*0.6, "4:1", cex=1.0)

dev.off()

par(mfrow=c(1,1), mar=c(4,4,1,1), oma=c(0,0,0,0),mai=c(0.9,0.9,0.2,0.2))
plot(NO3.~DOC_molar, data=data.in, col=col,
     ylim=c(0,50), xlim=c(0,200), #omits a few points
     pch=16,
     #    pch=as.numeric(droplevels(data.in$local_site)),
     ylab=expression(paste('[NO'['3']*'] (',mu, "mol N/L)")),
     xlab=expression(paste("[DOC] (",mu, "mol C/L)")),
     cex=0.6)
labs <- c('glacier', 'high', 'low')
legend('topright',labs, bty='n', cex=1.2,
       pch=16, col=c(4,1,2))
text(100,150, paste('n = ', length(data.in$local_site)))
abline(0,(1/4), lty=2)
text(620*0.65, 140*0.6, "4:1", cex=1.0)

data.in$DOC_NO3. <- data.in$DOC_molar/data.in$NO3.
data.in$intyear <- as.integer(data.in$decY)

plot(DOC_NO3.~DOY, data=data.in, col=col, ylim=c(0,50))
abline(h=4)

plot(DOC_NO3.~ ele2, data=data.in, c, log='y')

# --quick look at trends --
fout   <- paste('DOC_NO3.ann_trends_box.pdf', sep='')
pdf(fout, width=6, height=4) 
par(mfrow=c(3,1), mar=c(0,5,3,2))
boxplot(DOC_NO3.[data.in$ele2=='Snow']~intyear[data.in$ele2=='Snow'], 
        data=data.in, ylim=c(0,10), col=4, ylab=NA, xaxt='n')
abline(h=4, lty=2)

par(mar=c(1.5,5,1.5,2))
boxplot(DOC_NO3.[data.in$ele2=='High']~intyear[data.in$ele2=='High'], xaxt='n', 
        data=data.in, ylim=c(0,10), ylab=expression("DOC:NO"['3']*' molar'))
abline(h=4, lty=2)

par(mar=c(3,5,0,2))
boxplot(DOC_NO3.[data.in$ele2=='Low']~intyear[data.in$ele2=='Low'], 
        data=data.in, ylim=c(0,75), col=2)
abline(h=4, lty=2)

dev.off()

# !! USER BEWARE!!  
# following analyses are naive attempts to decompose seasonal cycle and annual trens 
#------------------------------------------------------
#  All data, DOC_NO3.
#------------------------------------------------------
# smoothing functiong by DOY and year
# increasing complexity of covariates
# unsure of gamm structure, 
# should we try and predict fixed effects by local_site or elevation bin
# how get random effects of local_site, here in correlation
# significant positive trend low elevation, 
# surprising not a decreasing trend in higher elevation?
# with 3 elevtions classes, no seasonal trends, ony positive trend w/ low elevation
AllDOC_NO3. <- data.frame(DOC_NO3.=data.in$DOC_NO3., DOY=data.in$DOY, decY=data.in$decY, 
                     local_site=data.in$local_site, ele=data.in$ele, ele2=data.in$ele2)
AllDOC_NO3. <- na.exclude(AllDOC_NO3.)


MOD4_DOC_NO3. <- gamm(DOC_NO3. ~ ele2 + s(DOY,bs="cc",by=ele2) + s(decY,by=ele2), data = AllDOC_NO3.,
                 correlation = corCAR1(form = ~ decY|local_site))

summary(MOD4_DOC_NO3.$gam)
plot(MOD4_DOC_NO3.$gam , residuals=TRUE,pch=16, cex=0.75,  pages=1, ylim=c(-50,50))


#------------------------------------------------------
#  All data, NO3.
#------------------------------------------------------
# model selection is really slow here
# worth doing correctly, but I'm unsure of the best way to proceed?
# models 0-4 here all fit to local_site, this is likely not necessary?

#MOD0_NO3. <- gamm(NO3. ~ local_site + s(DOY,bs="cc") + s(decY), data = AllNO3.,
#                 correlation = corCAR1(form = ~ decY|local_site))
#MOD1_NO3. <- gamm(NO3. ~ local_site + s(DOY,bs="cc", by=ele) + s(decY), data = AllNO3.,
#                 correlation = corCAR1(form = ~ decY|local_site))
#MOD2_NO3. <- gamm(NO3. ~ local_site + s(DOY,bs="cc") + s(decY,by=ele), data = AllNO3.,
#                 correlation = corCAR1(form = ~ decY|local_site))
#MOD3_NO3. <- gamm(NO3. ~ local_site + s(DOY,bs="cc",by=ele) + s(decY,by=ele), data = AllNO3.,
#                 correlation = corCAR1(form = ~ decY|local_site))
#MOD4_NO3. <- gamm(NO3. ~ local_site + s(DOY,bs="cc",by=ele) + s(decY,bs="ts",by=ele), 
#                 data = AllNO3., correlation = corCAR1(form = ~ decY|local_site))
#anova(MOD0_NO3.$lme, MOD1_NO3.$lme, MOD2_NO3.$lme, MOD3_NO3.$lme, MOD4_NO3.$lme)

# -- more justified to fit data based on elevation classes?
MOD4_NO3. <- gamm(NO3. ~ ele + s(DOY,bs="cc",by=ele) + s(decY,by=ele), data = AllNO3.,
                 correlation = corCAR1(form = ~ decY|local_site))
# Seasonal dynamics seem different in GL5_RG, Navajo, & Arikaree
# with highest NO3. concentrations at end of season?
MOD4b_NO3. <- gamm(NO3. ~ ele2 + s(DOY,bs="cc",by=ele2) + s(decY,by=ele2), data = AllNO3.,
                 correlation = corCAR1(form = ~ decY|local_site))

anova(MOD4_NO3.$lme, MOD4b_NO3.$lme)
# seasonal cycle of 'snow' different than 'high' and 'low' sites
# added complexity of breaking out 'snow' systems seems justified

summary(MOD4_NO3.$gam)
summary(MOD4_NO3.$gam)[[8]][[1]]
plot(MOD4_NO3.$gam , residuals=TRUE,pch=16, cex=0.75,  pages=1)
plot(MOD4b_NO3.$gam, residuals=TRUE,pch=16, cex=0.75,  pages=1)
with(AllNO3., tsDiagGamm(MOD4_NO3., decY, observed = NO3.))

fout   <- paste('NO3._trends.pdf', sep='')
pdf(fout, width=6, height=4)

par(mar=c(7,4,1,1), mfrow=c(1,2))
boxplot(AllNO3.$NO3.~droplevels(AllNO3.$ele),ylab='[NO3.]', las=2)
boxplot(AllNO3.$NO3.~droplevels(AllNO3.$ele2),ylab='[NO3.]',las=2)
dev.next()

par(mfrow=c(3,1), mar=c(0,5,3,2))
boxplot(AllNO3.$NO3.[AllNO3.$ele2=='Snow']~AllNO3.$year[AllNO3.$ele2=='Snow'], 
        data=AllNO3., col=4, ylab=NA, xaxt='n', ylab=NA, log='y')
abline(h=mean(AllNO3.$NO3.[AllNO3.$ele2=='Snow'], na.rm=T), lty=2)

par(mar=c(1.5,5,1.5,2))
boxplot(AllNO3.$NO3.[AllNO3.$ele2=='High']~AllNO3.$year[AllNO3.$ele2=='High'], 
        data=AllNO3., col=2, ylab=NA, xaxt='n', ylab='NO3_high', log='y')
abline(h=mean(AllNO3.$NO3.[AllNO3.$ele2=='High'], na.rm=T), lty=2)

par(mar=c(3,5,0,2))
boxplot(AllNO3.$NO3.[AllNO3.$ele2=='Low']~AllNO3.$year[AllNO3.$ele2=='Low'], 
        ylab='NO3_low',ylab=NA, log='y')
abline(h=mean(AllNO3.$NO3.[AllNO3.$ele2=='Low'], na.rm=T), lty=2)
dev.next()

par(mfrow=c(1,1), mar=c(4,5,1,1))
col <- as.numeric(AllNO3.$ele)
plot(AllNO3.$DOY,AllNO3.$NO3., 
     col=col, pch=as.numeric(droplevels(AllNO3.$local_site)),
     ylab="[NO3.] (uEQ/L)", xlab='Julian Day', cex=0.6)
labs <- levels(droplevels(AllNO3.$local_site)) 
legend('topleft',labs, bty='n', cex=0.8,
       pch=seq(1,length(text),1),
       col=c(2,2,1,2,1,1,1,2,1,2,2,2))
text(330,220, paste('n = ', length(AllNO3.$local_site)))
dev.next()

col <- as.numeric(AllNO3.$ele2)
col[col == 3] <- 4
plot(AllNO3.$DOY,AllNO3.$NO3., 
     col=col, pch=as.numeric(droplevels(AllNO3.$local_site)),
     ylab="[NO3.] (uEQ/L)", xlab='Julian Day', cex=0.6)
labs <- levels(droplevels(AllNO3.$local_site)) 
legend('topleft',labs, bty='n', cex=0.8,
       pch=seq(1,length(text),1),
       col=c(2,2,4,2,1,1,4,2,4,2,2,2))
text(330,220, paste('n = ', length(AllNO3.$local_site)))
dev.next()

ylim <- c(-15,20.)
par(mfrow=c(3,2), mar=c(1,5,5,1))
plot(MOD4b_NO3.$gam , residuals=TRUE,pch=16, cex=0.75,main="NO3. (intra-annual, snow)",  
     ylim=ylim, select=3, ylab='DOY, snow ', xlab="", xaxt='n')       
abline(0,0, lty=2, col=2)
plot(MOD4b_NO3.$gam , residuals=TRUE,pch=16, cex=0.75,main="NO3. (inter-annual, snow)",
     ylim=ylim, select=6, ylab='Year, snow',xlab="", xaxt='n')       
abline(0,0, lty=2, col=2)
p <- summary(MOD4b_NO3.$gam)[[8]][[6]]
text(2012, 12,paste('p = ',signif(p, 2)))

par(mar=c(2.5,5,2.5,1))
plot(MOD4b_NO3.$gam , residuals=TRUE,pch=16, cex=0.75, main="NO3. (intra-annual, high)",  
     ylim=ylim, select=1, ylab='DOY, High elevation', xlab="", xaxt='n')       
abline(0,0, lty=2, col=2)
plot(MOD4b_NO3.$gam , residuals=TRUE,pch=16, cex=0.75,main="NO3. (inter-annual, high)",  
     ylim=ylim, select=4, ylab='Year, High elevation', xlab="", xaxt='n')       
abline(0,0, lty=2, col=2, xlab="")
p <- summary(MOD4b_NO3.$gam)[[8]][[4]]
text(2012, 12,paste('p = ',signif(p, 2)))

par(mar=c(5,5,1,1))
plot(MOD4b_NO3.$gam , residuals=TRUE,pch=16, cex=0.75,main="NO3. (intra-annual, low)",  
     ylim=ylim, select=2, ylab='DOY, Low elevation', xlab="Day of Year")       
abline(0,0, lty=2, col=2)
plot(MOD4b_NO3.$gam , residuals=TRUE,pch=16, cex=0.75,main="NO3. (inter-annual, low)",
     ylim=ylim, select=5, ylab='Year, Low elevation',xlab="Year")       
abline(0,0, lty=2, col=2)
p <- summary(MOD4b_NO3.$gam)[[8]][[5]]
text(2012, 12,paste('p = ',signif(p, 2)))

dev.off()
#------------------------------------------------------


#------------------------------------------------------
#  All data, DOC
#------------------------------------------------------
# increasing complexity of covariates
# As with NO3., fitting to local_site results seems like overkill, 
#but it does have a lower AIC, BIC.
MOD_DOC <- gamm(DOC ~ ele + s(DOY,bs="cc") + s(decY), data = AllDOC,
                 correlation = corCAR1(form = ~ decY|local_site))
MOD0_DOC <- gamm(DOC ~ local_site + s(DOY,bs="cc") + s(decY), data = AllDOC,
                 correlation = corCAR1(form = ~ decY|local_site))
anova(MOD_DOC$lme, MOD0_DOC$lme)
#MOD1_DOC <- gamm(DOC ~ local_site + s(DOY,bs="cc", by=ele) + s(decY), data = AllDOC,
#                 correlation = corCAR1(form = ~ decY|local_site))
#MOD2_DOC <- gamm(DOC ~ local_site + s(DOY,bs="cc") + s(decY,by=ele), data = AllDOC,
#                 correlation = corCAR1(form = ~ decY|local_site))
#MOD2b_DOC <- gamm(DOC ~ local_site + s(DOY,bs="cc") + s(decY,bs="ts",by=ele), 
#                  data = AllDOC, correlation = corCAR1(form = ~ decY|local_site))
#MOD3_DOC <- gamm(DOC ~ local_site + s(DOY,bs="cc",by=ele) + s(decY,by=ele), data = AllDOC,
#                 correlation = corCAR1(form = ~ decY|local_site))
#MOD3b_DOC <- gamm(DOC ~ local_site + s(DOY,bs="cc",by=ele) + s(decY,bs="ts",by=ele), 
#                  data = AllDOC, correlation = corCAR1(form = ~ decY|local_site))
#anova(MOD0_DOC$lme, MOD1_DOC$lme, MOD2_DOC$lme,MOD2b_DOC$lme, MOD3_DOC$lme, MOD3b_DOC$lme)
#anova(MOD3_DOC$lme, MOD3b_DOC$lme, MOD4_DOC$lme)

MOD4_DOC <- gamm(DOC ~ ele + s(DOY,bs="cc",by=ele) + s(decY,by=ele), data = AllDOC,
                 correlation = corCAR1(form = ~ decY|local_site))

MOD4b_DOC <- gamm(DOC ~ ele2 + s(DOY,bs="cc",by=ele2) + s(decY,by=ele2), data = AllDOC,
                 correlation = corCAR1(form = ~ decY|local_site))

anova(MOD4_DOC$lme, MOD4b_DOC$lme)
# seasonal cycle of 'snow' NOT different than 'high' and 'low' sites
# added complexity of breaking out 'snow' systems does NOT seem justified
# Nearly a significant high elevation trend... maybe more data will help? 
summary(MOD4_DOC$gam)
plot(MOD4_DOC$gam , residuals=TRUE,pch=16, cex=0.75,  pages=1)
with(AllDOC, tsDiagGamm(MOD4_DOC, decY, observed = DOC))

par(mfrow=c(2,2), mar=c(4,5,1,1))
summary(MOD4_DOC$gam)
plot(MOD4_DOC$gam , residuals=TRUE,pch=16, cex=0.75,  pages=1)

fout   <- paste('DOC_trends.pdf', sep='')
pdf(fout, width=6, height=4)
par(mar=c(7,4,1,1), mfrow=c(1,2))
boxplot(AllDOC$DOC~droplevels(AllDOC$ele),ylab='[DOC]', las=2)
boxplot(AllDOC$DOC~droplevels(AllDOC$ele2),ylab='[DOC]',las=2)
dev.next()

par(mfrow=c(2,1), mar=c(0,5,3,2))
boxplot(AllDOC$DOC[AllDOC$ele=='High']~AllDOC$year[AllDOC$ele=='High'], 
        data=AllDOC, col=2, ylab=NA, xaxt='n', ylab='DOC', log='y')
abline(h=mean(AllDOC$DOC[AllDOC$ele=='High'], na.rm=T), lty=2)

par(mar=c(3,5,0,2))
boxplot(AllDOC$DOC[AllDOC$ele=='Low']~AllDOC$year[AllDOC$ele=='Low'], 
        ylab='DOC_low',ylab=NA, log='y')
abline(h=mean(AllDOC$DOC[AllDOC$ele=='Low'], na.rm=T), lty=2)
dev.next()


par(mfrow=c(1,1), mar=c(4,5,1,1))
col <- as.numeric(AllDOC$ele)
plot(AllDOC$DOY,AllDOC$DOC, 
     col=col, pch=as.numeric(droplevels(AllDOC$local_site)),
     ylab="[DOC] (mg C/L)", xlab='Julian Day', cex=0.6)
lab  <- levels(droplevels(AllDOC$local_site)) 
legend('topleft',lab, bty='n', cex=0.8,
       pch=seq(1,length(text),1),
       col=c(2,2,1,2,1,1,1,2,1,2,2,2))
text(330,10, paste('n = ', length(AllDOC$local_site)))
dev.next()


col <- as.numeric(AllDOC$ele2)
col[col == 3] <- 4
plot(AllDOC$DOY,AllDOC$DOC, 
     col=col, pch=as.numeric(droplevels(AllDOC$local_site)),
     ylab="[DOC] (mg C/L)", xlab='Julian Day', cex=0.6)
lab  <- levels(droplevels(AllDOC$local_site)) 
legend('topleft',lab, bty='n', cex=0.8,
       pch=seq(1,length(text),1),
       col=c(2,2,4,2,1,1,4,2,4,2,2,2))
text(330,10, paste('n = ', length(AllDOC$local_site)))
dev.next()

ylim <- c(-1,2.)
par(mfrow=c(2,2), mar=c(4,5,1,1))
plot(MOD4_DOC$gam , residuals=TRUE,pch=16, cex=0.75, main="DOC (intra-annual, high)",  
     ylim=ylim, select=1, ylab='DOY, High elevation', xlab="")       
abline(0,0, lty=2, col=2)
plot(MOD4_DOC$gam , residuals=TRUE,pch=16, cex=0.75,main="DOC (inter-annual, high)",  
     ylim=ylim, select=3, ylab='Year, High elevation', xlab="")       
abline(0,0, lty=2, col=2, xlab="")
p <- summary(MOD4_DOC$gam)[[8]][[3]]
text(2012, 1,paste('p = ',signif(p, 2)))

plot(MOD4_DOC$gam , residuals=TRUE,pch=16, cex=0.75,main="DOC (intra-annual, low)",  
     ylim=ylim, select=2, ylab='DOY, Low elevation', xlab="Day of Year")       
abline(0,0, lty=2, col=2)
plot(MOD4_DOC$gam , residuals=TRUE,pch=16, cex=0.75,main="DOC (inter-annual, low)",
     ylim=ylim, select=4, ylab='Year, Low elevation',xlab="Year")       
abline(0,0, lty=2, col=2)
p <- summary(MOD4_DOC$gam)[[8]][[4]]
text(2012, 1,paste('p = ',signif(p, 2)))
dev.off()

#------------------------------------------------------
#  All data, DON
#------------------------------------------------------
# Trends of increasing DON are more obvious in low elevation sites,
# high elevation data also shows surprising spike in DON, is this real?

#MOD0_DON <- gamm(DON ~ ele + s(DOY,bs="cc") + s(decY), data = AllDON,
#                 correlation = corCAR1(form = ~ decY|local_site))
MOD0_DON <- gamm(DON ~ local_site + s(DOY,bs="cc") + s(decY), data = AllDON,
                 correlation = corCAR1(form = ~ decY|local_site))
MOD1_DON <- gamm(DON ~ local_site + s(DOY,bs="cc", by=ele) + s(decY), data = AllDON,
                 correlation = corCAR1(form = ~ decY|local_site))
MOD2_DON <- gamm(DON ~ local_site + s(DOY,bs="cc") + s(decY,by=ele), data = AllDON,
                 correlation = corCAR1(form = ~ decY|local_site))
MOD2b_DON <- gamm(DON ~ local_site + s(DOY,bs="cc") + s(decY,bs="ts",by=ele), 
                  data = AllDON, correlation = corCAR1(form = ~ decY|local_site))
MOD3_DON <- gamm(DON ~ local_site + s(DOY,bs="cc",by=ele) + s(decY,by=ele), data = AllDON,
                 correlation = corCAR1(form = ~ decY|local_site))
MOD3b_DON <- gamm(DON ~ local_site + s(DOY,bs="cc",by=ele) + s(decY,bs="ts",by=ele), 
                  data = AllDON, correlation = corCAR1(form = ~ decY|local_site))
MOD4_DON <- gamm(DON ~ ele + s(DOY,bs="cc",by=ele) + s(decY,by=ele), data = AllDON,
                 correlation = corCAR1(form = ~ decY|local_site))

anova(MOD0_DON$lme, MOD1_DON$lme, MOD2_DON$lme,MOD2b_DON$lme, MOD3_DON$lme, MOD3b_DON$lme)
summary(MOD2b_DON$gam)

summary(MOD4_DON$gam)
plot(MOD4_DON$gam , residuals=TRUE,pch=16, cex=0.75,  pages=1)
with(AllDON, tsDiagGamm(MOD4_DON, decY, observed = DON))

fout   <- paste('DON_trends.pdf', sep='')
pdf(fout, width=6, height=4)
par(mar=c(7,4,1,1), mfrow=c(1,2))
boxplot(AllDON$DON~droplevels(AllDON$ele),ylab='[DON]', las=2)
boxplot(AllDON$DON~droplevels(AllDON$ele2),ylab='[DON]',las=2)
dev.next()

par(mfrow=c(2,1), mar=c(0,5,3,2))
boxplot(AllDON$DON[AllDON$ele=='High']~AllDON$year[AllDON$ele=='High'], 
        data=AllDON, col=2, ylab=NA, xaxt='n', ylab='DON',log='y')
abline(h=mean(AllDON$DON[AllDON$ele=='High'], na.rm=T), lty=2)

par(mar=c(3,5,0,2))
boxplot(AllDON$DON[AllDON$ele=='Low']~AllDON$year[AllDON$ele=='Low'], 
        ylab='DON_low',ylab=NA, log='y')
abline(h=mean(AllDON$DON[AllDON$ele=='Low'], na.rm=T), lty=2)
dev.next()

par(mfrow=c(1,1), mar=c(4,5,1,1))
col <- as.numeric(AllDON$ele)
plot(AllDON$DOY,AllDON$DON, 
     col=col, pch=as.numeric(droplevels(AllDON$local_site)),
     ylab="[DON] (mg C/L)", xlab='Julian Day', cex=0.6)
lab  <- levels(droplevels(AllDON$local_site)) 
legend('topleft',lab, bty='n', cex=0.8,
       pch=seq(1,length(text),1),
       col=c(2,2,1,2,1,1,1,2,1,2,2,2))
text(330,55, paste('n = ', length(AllDON$local_site)))
dev.next()

col <- as.numeric(AllDON$ele2)
col[col == 3] <- 4
plot(AllDON$DOY,AllDON$DON, 
     col=col, pch=as.numeric(droplevels(AllDON$local_site)),
     ylab="[DON] (mg C/L)", xlab='Julian Day', cex=0.6)
lab  <- levels(droplevels(AllDON$local_site)) 
legend('topleft',lab, bty='n', cex=0.8,
       pch=seq(1,length(text),1),
       col=c(2,2,4,2,1,1,4,2,4,2,2,2))
text(330,55, paste('n = ', length(AllDON$local_site)))
dev.next()

ylim <- c(-2,5.)
par(mfrow=c(2,2), mar=c(4,5,1,1))
plot(MOD4_DON$gam , residuals=TRUE,pch=16, cex=0.75, main="DON (intra-annual, high)",  
     ylim=ylim, select=1, ylab='DOY, High elevation', xlab="")       
abline(0,0, lty=2, col=2)
plot(MOD4_DON$gam , residuals=TRUE,pch=16, cex=0.75,main="DON (inter-annual, high)",  
     ylim=ylim, select=3, ylab='Year, High elevation', xlab="")       
abline(0,0, lty=2, col=2, xlab="")
p <- summary(MOD4_DON$gam)[[8]][[3]]
text(2012, 4,paste('p = ',signif(p, 2)))

plot(MOD4_DON$gam , residuals=TRUE,pch=16, cex=0.75,main="DON (intra-annual, low)",  
     ylim=ylim, select=2, ylab='DOY, Low elevation', xlab="Day of Year")       
abline(0,0, lty=2, col=2)
plot(MOD4_DON$gam , residuals=TRUE,pch=16, cex=0.75,main="DON (inter-annual, low)",
     ylim=ylim, select=4, ylab='Year, Low elevation',xlab="Year")       
abline(0,0, lty=2, col=2)
p <- summary(MOD4_DON$gam)[[8]][[4]]
text(2012, 4,paste('p = ',signif(p, 2)))
dev.off()


#------------------------------------------------------
# Legacy code, should be checked
#------------------------------------------------------

plot(pH~ANC, data=data.in, col=as.numeric(CLASS), pch=16)
legend(500, 5.2, legend=levels(data.in$CLASS), bty="n", col=seq(1,4,1), pch=16, cex=1.5 )
plot(DOC~ANC, data=data.in, col=as.numeric(class), pch=16 )
plot(DOC~NO3., data=data.in, col=as.numeric(class), pch=16 )
legend(150, 8, legend=levels(data.in$CLASS), bty="n", col=seq(1,4,1), pch=16, cex=1.5 )

par(mfrow=c(4,1), mar=c(2,5,1,2))
plot(pH~decY, data=data.in, col=as.numeric(ele), pch=16)
plot(ANC~decY, data=data.in, col=as.numeric(ele), pch=16)
plot(SO4..~decY, data=data.in, col=as.numeric(ele), pch=16, log="y")
#OMIT GL5_RG from SO4.. analysis, they go through the roof! 
#____subset=c(local_site!="GL5_RG")_____
plot(DOC~decY, data=data.in, col=as.numeric(ele), pch=16)
legend(1994, 9, legend=levels(data.in$ele), bty="n", col=seq(1,4,1), pch=16, cex=1.5 )

plot(pH~DOY, data=data.in, col=as.numeric(ele), pch=16)
plot(ANC~DOY, data=data.in, col=as.numeric(ele), pch=16)
plot(SO4..~DOY, data=data.in, col=as.numeric(ele), pch=16, log="y")
#OMIT GL5_RG from SO4.. analysis, they go through the roof! 

#____subset=c(local_site!="GL5_RG")_____
plot(DOC~DOY, data=data.in, col=as.numeric(ele), pch=16)
legend(25, 9, legend=levels(data.in$ele), bty="n", col=seq(1,4,1), pch=16, cex=1.5 )


par(mfrow=c(4,1), mar=c(2,5,1,2))
plot(DOC~decY, data=data.in, col=as.numeric(ele), pch=16)
legend(1994, 9, legend=levels(data.in$ele), bty="n", col=seq(1,4,1), pch=16, cex=1.5 )
plot(DON~decY, data=data.in, col=as.numeric(ele), pch=16, ylim=c(0,50))
plot(PO4...~decY, data=data.in, col=as.numeric(ele), pch=16)#, ylim=c(0,50))
plot(Ca.. ~decY, data=data.in, col=as.numeric(ele), pch=16, log='y')#, ylim=c(0,50))
#plot(NO3.~decY, data=data.in, col=as.numeric(ele), pch=16)
#plot(NH4.~decY, data=data.in, col=as.numeric(ele), pch=16, ylim=c(0,50))

par(mfrow=c(4,1), mar=c(2,5,1,2))
plot(DOC~DOY, data=data.in, col=as.numeric(ele), pch=16)
legend(10, 9, legend=levels(data.in$ele), bty="n", col=seq(1,4,1), pch=16, cex=1.5 )
plot(DON~DOY, data=data.in, col=as.numeric(ele), pch=16, ylim=c(0,50))
plot(PO4...~DOY, data=data.in, col=as.numeric(ele), pch=16, ylim=c(0,50))
plot(Ca.. ~DOY, data=data.in, col=as.numeric(ele), pch=16, ylim=c(0,50))
plot(NO3.~DOY, data=data.in, col=as.numeric(ele), pch=16)
plot(NH4.~DOY, data=data.in, col=as.numeric(ele), pch=16, ylim=c(0,50))


par(mfrow=c(4,1), mar=c(2,5,2,2))
plot(NO3.~decY, data=data.in, col=as.numeric(ele), pch=16, subset=c(ele=="High"),
     main="high")
legend(1994, 250, legend=levels(data.in$ele), bty="n", col=seq(1,4,1), pch=16, cex=1.5 )
plot(DON~decY, data=data.in, col=as.numeric(ele), pch=16, subset=c(ele=="High"))
plot(NO3.~decY, data=data.in, col=as.numeric(ele), pch=16, subset=c(ele=="Low"),
     main="low")
plot(DON~decY, data=data.in, col=as.numeric(ele), pch=16, subset=c(ele=="Low"))

plot(NO3.~DOY, data=data.in, col=as.numeric(ele), pch=16, subset=c(ele=="High"),
     main="high")
legend(1994, 250, legend=levels(data.in$ele), bty="n", col=seq(1,4,1), pch=16, cex=1.5 )
plot(DON~DOY, data=data.in, col=as.numeric(ele), pch=16, subset=c(ele=="High"))
plot(NO3.~DOY, data=data.in, col=as.numeric(ele), pch=16, subset=c(ele=="Low"),
     main="low")
plot(DON~DOY, data=data.in, col=as.numeric(ele), pch=16, subset=c(ele=="Low"))

#-------------------------------------------

## nested random effects and within group correlation 
# annual trend covarries w/ elevation
MOD0_pH <- gamm(pH ~ s(DOY, bs = "cc") + s(decY, by=ele), data = data.in,
                random=list(local_site=~1), correlation = corCAR1(form = ~ decY))
#AR errors defined by local_site
#produced identical results as MOD0
MOD0b_pH <- gamm(pH ~ s(DOY, bs = "cc") + s(decY, by=ele), data = data.in,
                random=list(local_site=~1), correlation = corCAR1(form = ~ decY|local_site))
#adds local_site as a fixed effect
MOD1_pH <- gamm(pH ~ local_site + s(DOY, bs = "cc") + s(decY, by=ele), data = data.in,
              correlation = corCAR1(form = ~ decY|local_site))
#add random erros ~ local_site
MOD1a_pH <- gamm(pH ~ local_site + s(DOY, bs = "cc") + s(decY, by=ele), data = data.in,
                random=list(local_site=~1), correlation = corCAR1(form = ~ decY|local_site))
#add random erros ~ year
MOD1b_pH <- gamm(pH ~ local_site + s(DOY, bs = "cc") + s(decY, by=ele), data = data.in,
                 random=list(decY=~1), correlation = corCAR1(form = ~ decY|local_site))

MOD2_pH <- gamm(pH ~ local_site + s(DOY, bs = "cc") + s(decY, by=ele), data = data.in,
                correlation = corCAR1(form = ~ 1|local_site))
# annual cycle covaries by elevation too
MOD3_pH <- gamm(pH ~ local_site + s(DOY, bs = "cc", by=ele) + s(decY, by=ele), data = data.in,
                correlation = corCAR1(form = ~ decY|local_site))
# set smoothing term to thin plate regression spline (bs="tp", "ts" not better than before)
MOD4_pH <- gamm(pH ~ local_site + s(DOY, bs = "cc", by=ele) + s(decY, bs="ts", by=ele), 
                 data = data.in, correlation = corCAR1(form = ~ decY|local_site))
MOD4b_pH <- gamm(pH ~ local_site + s(DOY, bs = "cc", by=ELE) + s(decY, bs="ts", by=ELE), 
                data = data.in, correlation = corCAR1(form = ~ decY|local_site))
# set smoothing term a factor smooth (bs="fs")
MOD5_pH <- gamm(pH ~ local_site + s(DOY, bs = "cc", by=ele) + s(decY, bs="cc", by=ele), 
                data = data.in, correlation = corCAR1(form = ~ decY|local_site))
anova(MOD0_pH$lme, MOD0b_pH$lme, MOD1_pH$lme, MOD1a_pH$lme, MOD1b_pH$lme,
      MOD2_pH$lme, MOD3_pH$lme,  MOD4_pH$lme, MOD5_pH$lme) 
anova(MOD3_pH$lme,  MOD4_pH$lme, MOD4b_pH$lme, MOD5_pH$lme) 
#lowest AIC/BIC best (highest logLik)


#best model (MOD3) with local_site fixed effects, no random effects, AR errors ~ local_site
summary(MOD5_pH$gam)       

par(mfrow=c(2,1), mar=c(6,5,1,1))
acf(resid(MOD_pH$lme, type='normalized'))
pacf(resid(MOD_pH$lme, type='normalized'))

plot(MOD5_pH$gam , residuals=TRUE,pch=16, cex=0.75,  pages=1)

par(mfrow=c(2,2), mar=c(6,5,1,1))
plot(MOD3_pH$gam , residuals=TRUE,pch=16, cex=0.75,  
     ylim=c(-0.5,0.5), select=1, ylab='DOY, High elevation')       
abline(0,0, lty=2, col=2)
plot(MOD3_pH$gam , residuals=TRUE,pch=16, cex=0.75,  
     ylim=c(-0.5,0.5), select=2, ylab='DOY, Low elevation')       
abline(0,0, lty=2, col=2)
plot(MOD3_pH$gam , residuals=TRUE,pch=16, cex=0.75,  
     ylim=c(-0.5,0.5), select=3, ylab='Year, High elevation')       
abline(0,0, lty=2, col=2)
plot(MOD3_pH$gam , residuals=TRUE,pch=16, cex=0.75,  
     ylim=c(-0.5,0.5), select=4, ylab='Year, Low elevation')       
abline(0,0, lty=2, col=2)

print(MOD3_pH$gam)

#------------------------------------------------------
# Ca..n we account for site and annual differences, 
#       and look for trends in residuals?
#------------------------------------------------------
MODx_pH <- gamm(pH ~ local_site + s(DOY, bs = "cc", by=ele), data = data.in,
                correlation = corCAR1(form = ~ decY|local_site))
MODy_pH <- gamm(pH ~ local_site + s(DOY, bs = "cc", by=ele), data = data.in,
                correlation = corCAR1(form = ~ 1|local_site))
MODz_pH <- gamm(pH ~ local_site + s(DOY, bs = "cc"), data = data.in,
                correlation = corCAR1(form = ~ 1|local_site))
anova(MODx_pH$lme, MODy_pH$lme, MODz_pH$lme)
plot(MODy_pH$gam , residuals=TRUE,pch=16, cex=0.75,  pages=1)
Pred_pH   <- predict(MODx_pH$gam)
AllpH$Resid_pH  <- residuals(MODx_pH$gam)
plot(Resid_pH ~ decY, pch=16, col=as.numeric(ele), data=AllpH)


MODx_NO3. <- gamm(NO3. ~ local_site + s(DOY, bs = "cc", by=ele), data = data.in,
                correlation = corCAR1(form = ~ decY|local_site))
MODy_NO3. <- gamm(NO3. ~ local_site + s(DOY, bs = "cc", by=ele), data = data.in,
                correlation = corCAR1(form = ~ 1|local_site))
MODz_NO3. <- gamm(NO3. ~ local_site + s(DOY, bs = "cc"), data = data.in,
                correlation = corCAR1(form = ~ 1|local_site))
anova(MODx_NO3.$lme, MODy_NO3.$lme, MODz_NO3.$lme)
summary(MODy_NO3.$gam)
par(mfrow=c(2,1), mar=c(3,5,1,1))
plot(MODy_NO3.$gam, residual=TRUE)
AllNO3.$Pred_NO3.   <- predict(MODy_NO3.$gam)
AllNO3.$Resid_NO3.  <- residuals(MODy_NO3.$gam)
par(mfrow=c(1,1))
plot(Resid_NO3. ~ decY, pch=16, col=as.numeric(ele), data=AllNO3., 
     ylim=c(-50,75))

lm_NO3._high <- lm(Resid_NO3. ~ decY, data=AllNO3., subset=c(ele=="High"))
lm_NO3._low  <- lm(Resid_NO3. ~ decY, data=AllNO3., subset=c(ele=="Low"))
summary(lm_NO3._high)
summary(lm_NO3._low)
lines(predict(lm_NO3._high)~AllNO3.$decY[AllNO3.$ele=='High'], col=1)
lines(predict(lm_NO3._low)~AllNO3.$decY[AllNO3.$ele=='Low'], col=2)
#------------------------------------------------------


#------------------------------------------------------
#  All data, SO4..
#------------------------------------------------------
# annual cycle covaries by elevation too
MOD2_SO4.. <- gamm(SO4.. ~ ele + s(DOY, bs = "cc", by=ele) + s(decY, by=ele), data = data.in,
                 correlation = corCAR1(form = ~ decY|local_site), 
                 subset=c(local_site!="GL5_RG"))
MOD3_SO4.. <- gamm(SO4.. ~ local_site + s(DOY, bs = "cc", by=ele) + s(decY, by=ele), data = data.in,
                correlation = corCAR1(form = ~ decY|local_site), 
                subset=c(local_site!="GL5_RG"))
anova(MOD2_SO4..$lme, MOD3_SO4..$lme)
summary(MOD3_SO4..$gam)
plot(MOD3_SO4..$gam , residuals=TRUE,pch=16, cex=0.75,  pages=1, 
     ylim=c(-50,100))

par(mfrow=c(2,2), mar=c(6,5,1,1))
plot(MOD3_SO4..$gam , residuals=TRUE,pch=16, cex=0.75,  
     ylim=c(-50,75), select=1, ylab='DOY, High elevation')       
abline(0,0, lty=2, col=2)
plot(MOD3_SO4..$gam , residuals=TRUE,pch=16, cex=0.75,  
     ylim=c(-50,75), select=2, ylab='DOY, Low elevation')       
abline(0,0, lty=2, col=2)
plot(MOD3_SO4..$gam , residuals=TRUE,pch=16, cex=0.75,  
     ylim=c(-50,75), select=3, ylab='Year, High elevation')       
abline(0,0, lty=2, col=2)
plot(MOD3_SO4..$gam , residuals=TRUE,pch=16, cex=0.75,  
     ylim=c(-50,75), select=4, ylab='Year, Low elevation')       
abline(0,0, lty=2, col=2)
#------------------------------------------------------

#------------------------------------------------------
#  All data, Ca.. 
#------------------------------------------------------
# annual cycle covaries by elevation too
MOD2_Ca..  <- gamm(Ca..  ~ ele + s(DOY, bs = "cc", by=ele) + s(decY, by=ele), data = data.in,
                 correlation = corCAR1(form = ~ decY|local_site), 
                 subset=c(local_site!="GL5_RG"))
MOD3_Ca..  <- gamm(Ca..  ~ local_site + s(DOY, bs = "cc", by=ele) + s(decY, by=ele), data = data.in,
                 correlation = corCAR1(form = ~ decY|local_site), 
                 subset=c(local_site!="GL5_RG"))
anova(MOD2_Ca.. $lme, MOD3_Ca.. $lme)
summary(MOD3_Ca.. $gam)
plot(MOD3_Ca.. $gam , residuals=TRUE,pch=16, cex=0.75,  pages=1, 
     ylim=c(-50,100))

par(mfrow=c(2,2), mar=c(6,5,1,1))
plot(MOD3_Ca.. $gam , residuals=TRUE,pch=16, cex=0.75,  
     ylim=c(-50,75), select=1, ylab='DOY, High elevation')       
abline(0,0, lty=2, col=2)
plot(MOD3_Ca.. $gam , residuals=TRUE,pch=16, cex=0.75,  
     ylim=c(-50,75), select=2, ylab='DOY, Low elevation')       
abline(0,0, lty=2, col=2)
plot(MOD3_Ca.. $gam , residuals=TRUE,pch=16, cex=0.75,  
     ylim=c(-50,75), select=3, ylab='Year, High elevation')       
abline(0,0, lty=2, col=2)
plot(MOD3_Ca.. $gam , residuals=TRUE,pch=16, cex=0.75,  
     ylim=c(-50,75), select=4, ylab='Year, Low elevation')       
abline(0,0, lty=2, col=2)
#------------------------------------------------------

#------------------------------------------------------
#  All data, PO4...
#------------------------------------------------------
# annual cycle covaries by elevation too
MOD2_PO4... <- gamm(PO4... ~ ele + s(DOY, bs = "cc", by=ele) + s(decY, by=ele), data = data.in,
                 correlation = corCAR1(form = ~ decY|local_site), 
                 subset=c(local_site!="GL5_RG"))
MOD3_PO4... <- gamm(PO4... ~ local_site + s(DOY, bs = "cc", by=ele) + s(decY, by=ele), data = data.in,
                 correlation = corCAR1(form = ~ decY|local_site), 
                 subset=c(local_site!="GL5_RG"))
anova(MOD2_PO4...$lme, MOD3_PO4...$lme)
summary(MOD3_PO4...$gam)
plot(MOD3_PO4...$gam , residuals=TRUE,pch=16, cex=0.75,  pages=1)

par(mfrow=c(2,2), mar=c(6,5,1,1))
plot(MOD3_PO4...$gam , residuals=TRUE,pch=16, cex=0.75,  
     ylim=c(-50,75), select=1, ylab='DOY, High elevation')       
abline(0,0, lty=2, col=2)
plot(MOD3_PO4...$gam , residuals=TRUE,pch=16, cex=0.75,  
     ylim=c(-50,75), select=2, ylab='DOY, Low elevation')       
abline(0,0, lty=2, col=2)
plot(MOD3_PO4...$gam , residuals=TRUE,pch=16, cex=0.75,  
     ylim=c(-50,75), select=3, ylab='Year, High elevation')       
abline(0,0, lty=2, col=2)
plot(MOD3_PO4...$gam , residuals=TRUE,pch=16, cex=0.75,  
     ylim=c(-50,75), select=4, ylab='Year, Low elevation')       
abline(0,0, lty=2, col=2)
#------------------------------------------------------

#for (n in 1:nsites) {
n <- 1
c <- 1
  #omit sparse sites
  #if(sites[n]!="NIWOT"){

#  temp   <- subset(data.in, local_site==sites[n], select=c(DOY,decY,pH,DOC,SO4..,NO3.,NH4.,ANC))
  temp   <- subset(data.in, ele == "High", select=c(DOY,decY,pH,DOC,SO4..,NO3.,NH4.,DON,ANC,local_site))
#  temp   <- subset(data.in, local_site=="ALBION", select=c(DOY,decY,pH,DOC,SO4..,NO3.,NH4.,DON,ANC,local_site))
  nsteps <- length(temp$DOY)

	#-------------Subset data, remove NA------------------------------
  TEMPpH  <- data.frame(pH =temp$pH,  DOY=temp$DOY, decY=temp$decY, local_site=temp$local_site)
  TEMPSO4.. <- data.frame(SO4..=temp$SO4.., DOY=temp$DOY, decY=temp$decY, local_site=temp$local_site)
  TEMPDOC <- data.frame(DOC=temp$DOC, DOY=temp$DOY, decY=temp$decY, local_site=temp$local_site)
  TEMPDON <- data.frame(DON=temp$DON, DOY=temp$DOY, decY=temp$decY, local_site=temp$local_site)
  TEMPNH4. <- data.frame(NH4.=temp$NH4., DOY=temp$DOY, decY=temp$decY, local_site=temp$local_site)
  TEMPNO3. <- data.frame(NO3.=temp$NO3., DOY=temp$DOY, decY=temp$decY, local_site=temp$local_site)
  TEMPANC <- data.frame(ANC=temp$ANC, DOY=temp$DOY, decY=temp$decY, local_site=temp$local_site)
  TEMPpH  <- na.exclude(TEMPpH)
	TEMPSO4.. <- na.exclude(TEMPSO4..)
	TEMPDOC <- na.exclude(TEMPDOC)
  TEMPDON <- na.exclude(TEMPDON)
  TEMPNH4. <- na.exclude(TEMPNH4.)
	TEMPNO3. <- na.exclude(TEMPNO3.)
  TEMPANC <- na.exclude(TEMPANC)
  
  TEMPpH$local_site  <- factor(TEMPpH$local_site)
  TEMPSO4..$local_site <- factor(TEMPSO4..$local_site)
  TEMPDOC$local_site <- factor(TEMPDOC$local_site)
  TEMPDON$local_site <- factor(TEMPDON$local_site)
  TEMPNH4.$local_site <- factor(TEMPNH4.$local_site)
  TEMPNO3.$local_site <- factor(TEMPNO3.$local_site)
  TEMPANC$local_site <- factor(TEMPANC$local_site)

#---------------look at seasonal cycle and annual trends------------
	#-----------------------pH------------------------------------------
p <- ggplot(TEMPpH, aes(x = decY, y = pH)) + geom_point()
print(p)
p + stat_smooth(method = "gam", formula = y ~ s(DOY, bs = "cc") + s(decY), size = 1)

m_ucm <- ucm(formula = pH~0, data = TEMPpH, level = TRUE)
m_ucm
levels(TEMPpH$local_site)



#-----------------------DOC------------------------------------------
p <- ggplot(TEMPDOC, aes(x = decY, y = DOC)) + geom_point()
print(p)
p + stat_smooth(method = "gam", formula = y ~ s(DOY, bs = "cc") + s(decY), size = 1)

m_ucm <- ucm(formula = pH~0, data = TEMPpH, level = TRUE)
m_ucm

xlim <- c(1994,2015)
main <- levels(data.in$CLASS)[c]
main <- levels(data.in$ele)[c]
par(mfrow=c(4,1), mar=c(1,4,3,1))
plot(pH ~ decY,  data=TEMPpH,  type="p", ylab="pH", col=as.numeric(local_site), 
     xlim=xlim, main=main)
par(mar=c(2,4,2,1))
plot(ANC ~ decY, data=TEMPANC, type="p", ylab="ANC", col=as.numeric(local_site),
     xlim=xlim)
plot(SO4.. ~ decY, data=TEMPSO4.., type="p", ylab="SO4..", col=as.numeric(local_site), 
     xlim=xlim, log="y")
plot(DOC ~ decY, data=TEMPDOC, type="p", ylab="DOC", col=as.numeric(local_site), 
     xlim=xlim)

par(mfrow=c(1,1))
plot(DOC ~ decY, data=TEMPDOC, type="p", ylab="DOC", col=as.numeric(local_site), 
     xlim=xlim, main=main)
nsite <- length(levels(TEMPDOC$local_site))
legend('topleft',levels(TEMPDOC$local_site), text.col=seq(1,nsite,1))
TEMPNO3.$NO3.[TEMPNO3.$NO3. > 200] <- NA
plot(NO3. ~ decY, data=TEMPNO3., type="p", ylab="NO3.", col=as.numeric(local_site), 
     xlim=xlim, main=main)
nsite <- length(levels(TEMPDOC$local_site))
legend('topleft',levels(TEMPDOC$local_site), text.col=seq(1,nsite,1))
TEMPDON$DON[TEMPDON$DON > 50] <- NA
plot(DON ~ decY, data=TEMPDON, type="p", ylab="DON", col=as.numeric(local_site), 
     xlim=xlim, main=main)
nsite <- length(levels(TEMPDOC$local_site))
legend('topleft',levels(TEMPDOC$local_site), text.col=seq(1,nsite,1))

par(mfrow=c(3,1))
plot(DOC ~ decY, data=TEMPDOC, type="p", ylab="DOC", col=as.numeric(local_site), 
     xlim=xlim, main=main)
nsite <- length(levels(TEMPDOC$local_site))
legend('topleft',levels(TEMPDOC$local_site), text.col=seq(1,nsite,1))
plot(DOC ~ decY, data=TEMPDOC, type="p", ylab="DOC", col=as.numeric(local_site), 
     xlim=xlim, main=main)
plot(DON ~ decY, data=TEMPDON, type="p", ylab="DON", col=as.numeric(local_site),
     xlim=xlim, log="y")
plot(NO3. ~ decY, data=TEMPNO3., type="p", ylab="NO3.", col=as.numeric(local_site),
     xlim=xlim, log="y")

plot(NH4. ~ decY, data=TEMPNH4., type="p", ylab="NH4.", col=as.numeric(local_site),
     xlim=xlim, log="y")

m0_pH <- gamm(pH ~  s(DOY, bs = "cc") + s(decY), data = TEMPpH,
              correlation = corCAR1(form = ~ decY|local_site))
m1_pH <- gamm(pH ~  s(DOY, bs = "cc", by=local_site) + s(decY, by=local_site), data = TEMPpH,
              correlation = corCAR1(form = ~ decY|local_site))
m2_pH <- gamm(pH ~ local_site + s(DOY, bs = "cc", by=local_site) + s(decY, by=local_site), data = TEMPpH,
              correlation = corCAR1(form = ~ decY|local_site))
m3_pH <- gamm(pH ~ local_site + s(DOY, bs = "cc") + s(decY), data = TEMPpH,
              correlation = corCAR1(form = ~ decY|local_site))
anova(m0_pH$lme, m1_pH$lme, m2_pH$lme, m3_pH$lme) #lowest AIC best (highest logLik)
summary(m3_pH$gam)       

par(mfrow=c(2,1), mar=c(6,5,1,1))
	acf(resid(m3_pH$lme, type='normalized'))
	pacf(resid(m3_pH$lme, type='normalized'))

plot(m3_pH$gam , residuals=TRUE,pch=16, cex=0.75, main=main, pages=1)       
  
plot(m3_pH$gam ,pch=16, cex=0.75)       
abline(0,0, lty=2, col=2)
print(m1_pH$gam)

with(TEMPpH, tsDiagGamm(m3_pH, decY, observed = pH))

pH1  <- predict(m3_pH$gam)

#-----------------------SO4..------------------------------------------
	quartz()
	plot(SO4.. ~ decY, data=TEMPSO4.., type="o", ylab="SO4..")
	m0_SO4.. <- gamm(SO4.. ~ s(DOY, bs = "cc") , data = TEMPSO4..,
            correlation = corCAR1(form = ~ decY))
	m1_SO4.. <- gamm(SO4.. ~ s(DOY, bs = "cc") + s(decY), data = TEMPSO4..,
            correlation = corCAR1(form = ~ decY))
	m2_SO4.. <- gamm(SO4.. ~ s(decY, k=30), data = TEMPSO4..,
            correlation = corCAR1(form = ~ decY))

    anova(m0_SO4..$lme, m1_SO4..$lme, m2_SO4..$lme)
    summary(m1_SO4..$gam)       

    par(mfrow=c(2,1), mar=c(6,5,1,1))
	acf(resid(m1_SO4..$lme, type='normalized'))
	pacf(resid(m1_SO4..$lme, type='normalized'))

    par(mfrow=c(2,1), mar=c(6,5,1,1))
    plot(m1_SO4..$gam, residuals=TRUE,pch=16, cex=0.75, main=sites[n])       
    abline(0,0, lty=2)

	with(TEMPSO4.., tsDiagGamm(m1_SO4.., decY, observed = SO4..))

	SO4..1  <- predict(m1_SO4..$gam)

	#-----------------------DOC------------------------------------------
	quartz()
	plot(DOC ~ decY, data=TEMPDOC, type="o", ylab="DOC")
	m0_DOC <- gamm(DOC ~ s(DOY, bs = "cc") , data = TEMPDOC,
            correlation = corCAR1(form = ~ decY))
	m1_DOC <- gamm(DOC ~ s(DOY, bs = "cc") + s(decY), data = TEMPDOC,
            correlation = corCAR1(form = ~ decY))
	m2_DOC <- gamm(DOC ~ s(decY, k=30), data = TEMPDOC,
            correlation = corCAR1(form = ~ decY))

    DOC_anova <- anova(m0_DOC$lme, m1_DOC$lme, m2_DOC$lme)
    DOC_anova 
    aAIC 	  <- DOC_anova[[4]] + 2*DOC_anova[[3]]
    (min(aAIC))
    summary <- data.frame(model = DOC_anova[[2]]-1, aAIC = aAIC)
	best <- summary$model[summary$aAIC==min(summary$aAIC)]
	
    summary(m0_DOC$gam)       
    summary(m1_DOC$gam)       

    par(mfrow=c(2,1), mar=c(6,5,1,1))
	acf(resid(m1_DOC$lme, type='normalized'))
	pacf(resid(m1_DOC$lme, type='normalized'))

    par(mfrow=c(2,1), mar=c(6,5,1,1))
    plot(m1_DOC$gam, residuals=TRUE,pch=16, cex=0.75,main="")       
    abline(0,0, lty=2, col=2)

	quartz()  #diagnostic plots
	with(TEMPDOC, tsDiagGamm(m1_DOC, decY, observed = DOC))

	DOC1  <- predict(m1_DOC$gam)

#-----------------------NO3.------------------------------------------
quartz()
plot(NO3. ~ decY, data=TEMPNO3., type="o", ylab="NO3.")
m0_NO3. <- gamm(NO3. ~ s(DOY, bs = "cc") , data = TEMPNO3.,
               correlation = corCAR1(form = ~ decY))
m1_NO3. <- gamm(NO3. ~ s(DOY, bs = "cc") + s(decY), data = TEMPNO3.,
               correlation = corCAR1(form = ~ decY))
m2_NO3. <- gamm(NO3. ~ s(decY, k=30), data = TEMPNO3.,
               correlation = corCAR1(form = ~ decY))

anova(m0_NO3.$lme, m1_NO3.$lme, m2_NO3.$lme)
summary(m1_NO3.$gam)       

par(mfrow=c(2,1), mar=c(6,5,1,1))
acf(resid(m1_NO3.$lme, type='normalized'))
pacf(resid(m1_NO3.$lme, type='normalized'))

par(mfrow=c(2,1), mar=c(6,5,1,1))
plot(m1_NO3.$gam, residuals=TRUE,pch=16, cex=0.75, main=sites[n])       
abline(0,0, lty=2)

with(TEMPNO3., tsDiagGamm(m1_NO3., decY, observed = NO3.))

NO3.1  <- predict(m1_NO3.$gam)


    fout <- paste("trends/",sites[n],".pdf", sep="")
    pdf(fout, width=7, height=7)
    range <- (c(1994,2013))
    
	par(mfrow=c(3,1), omi=c(0,0,0,0), mar=c(0,5,4,1),cex=0.9)
		plot(pH ~ decY, data=TEMPpH, main=sites[n], ylab="pH", xlab="",xaxt="n",xlim=range, 
			pch=16)
		lines(pH1~decY, data=TEMPpH,col="blue",lwd=2)
#		lines(pH2~decY, data=TEMPpH,col="red",lwd=2)
	par(mar=c(2,5,2,1))
		plot(SO4.. ~ decY, data=TEMPSO4.., ylab="SO4.. (uEQ/L)", xlab="", xaxt="n",
			xlim=range,pch=16)
		lines(SO4..1~ decY, data=TEMPSO4..,col="blue",lwd=2)
	par(mar=c(4,5,0,1))
		plot(DOC~ decY, data=TEMPDOC, ylab="DOC (mg C/L)", xlab="year",xlim=range, pch=16)
		lines(DOC1~ decY, data=TEMPDOC,col="blue",lwd=2)
	dev.off()

	remove(temp, fout)

  }
  }
  }
}

  NR1_date <- as.Date(Data.obs$DateTime)
  NR1_date[1:50]	
  
  NR1.year   <- tapply(Data.obs$Year, NR1_date, mean)
  NR1.mo     <- tapply(Data.obs$MO,   NR1_date, mean)
  NR1.dd     <- tapply(Data.obs$DD,   NR1_date, mean)
  NR1.ppt    <- tapply(Data.obs$precip_mm,  NR1_date, sum)
  NR1.tsa    <- tapply(Data.obs$T_2m, NR1_date, mean)

  remove(Data.obs)

  length(Data.sno$ST_swe)

#-------------Look at monthly mean and max SWE & snow Depth--------------------------
  nmonths  <- 12 * nyear

  ST_swe_max     <- rep(NA, nmonths)
  ST_swe_mean    <- rep(NA, nmonths)
  ST_depth_max   <- rep(NA, nmonths)
  ST_depth_mean  <- rep(NA, nmonths)
  ST_ppt_sum     <- rep(NA, nmonths)
  CLM_swe_max    <- rep(NA, nmonths)
  CLM_swe_mean   <- rep(NA, nmonths)
  CLM_depth_max  <- rep(NA, nmonths)
  CLM_depth_mean <- rep(NA, nmonths)
  CLM_snowdp_max  <- rep(NA, nmonths)
  CLM_snowdp_mean <- rep(NA, nmonths)
  CLM_ppt_sum     <- rep(NA, nmonths)
  CLM_snow_sum    <- rep(NA, nmonths)

  i <- 1
  
  for (y in 1:nyear) {
  	for (m in 1:12)   {
  		ST_swe_max[i]   <- max( Data.sno$ST_swe[Data.sno$Year==year[y]&Data.sno$mo==m], na.rm=T)
  		ST_swe_mean[i]  <- mean(Data.sno$ST_swe[Data.sno$Year==year[y]&Data.sno$mo==m], na.rm=T)
  		ST_depth_max[i] <- max( Data.sno$ST_depth[Data.sno$Year==year[y]&Data.sno$mo==m], na.rm=T)
  		ST_depth_mean[i]<- mean(Data.sno$ST_depth[Data.sno$Year==year[y]&Data.sno$mo==m], na.rm=T)
  		ST_ppt_sum[i]   <- sum(Data.sno$ST_ppt[Data.sno$Year==year[y]&Data.sno$mo==m], na.rm=T)
  		  
  		CLM_swe_max[i]     <- max( SWE_mean[Data.sno$Year==year[y] & Data.sno$mo==m], na.rm=T)
  		CLM_swe_mean[i]    <- mean(SWE_mean[Data.sno$Year==year[y] & Data.sno$mo==m], na.rm=T)
  		CLM_depth_max[i]   <- max( SNOW_DEPTH_mean[Data.sno$Year==year[y]&Data.sno$mo==m], na.rm=T)
  		CLM_depth_mean[i]  <- mean(SNOW_DEPTH_mean[Data.sno$Year==year[y]&Data.sno$mo==m],na.rm=T)
  		CLM_snowdp_max[i]  <- max( SNOWDP_mean[Data.sno$Year==year[y]&Data.sno$mo==m], na.rm=T)
  		CLM_snowdp_mean[i] <- mean(SNOWDP_mean[Data.sno$Year==year[y]&Data.sno$mo==m], na.rm=T)
  		CLM_ppt_sum[i]     <- sum(     PPT_sum[Data.sno$Year==year[y]&Data.sno$mo==m], na.rm=T)
  		CLM_snow_sum[i]    <- sum(    SNOW_sum[Data.sno$Year==year[y]&Data.sno$mo==m], na.rm=T)

  		  i <- i + 1
  	}
  }
  
  mon   <- seq(8,19,1)
  month <- c('A','S','O','N','D','J','F','M','A','M','J','J')
  
# modify HERE for CRNS
#   fout <- paste("NR1_figs/monthly_SNOWTEL_CLM.pdf", sep="")		#tower precip
    fout <- paste("NR1_figs_v2/monthly_SNOWTEL_CLM_v2.pdf", sep="") 	#CRNS  precip
    pdf(fout, width=9, height=7)
    par(mfrow=c(4,3), mar=c(0,0,0,0), oma=c(5,5,3,2))

	start <- 8 + 12
	for (y in 2:(nyear-1)) {
		print(year[y])
		end <- start + 11

		#calculate cumulative precipitation
		ann_CLM_SNOW<- CLM_snow_sum[start:end]
		ann_CLM_PPT <- CLM_ppt_sum[ start:end]
		ann_ST_PPT  <-  ST_ppt_sum[ start:end]
        cum_CLM_SNOW<- rep(NA,12)
        cum_CLM_PPT <- rep(NA,12)
        cum_ST_PPT  <- rep(NA,12)
		for (t in 1:12) {
			cum_CLM_SNOW[t] <- sum(ann_CLM_SNOW[1:t])
			cum_CLM_PPT[t]  <- sum(ann_CLM_PPT[ 1:t])
			cum_ST_PPT[t]   <- sum( ann_ST_PPT[ 1:t])
		}

		plot( mon,ST_swe_max[start:end], type="l", lwd=2, ylim=c(0,max(ST_swe_max, na.rm=T)),
		  xaxt = "n", yaxt = "n", xlab=NA, ylab="SWE (mm)")
   	    lines(mon,CLM_swe_max[start:end], col=2, lwd=2) #convert kg/m2*1m3/1000kg*1000 mm/m
		lines(mon, cum_CLM_SNOW,          col=2, lwd=2, lty=2)
		text((mon[1]*1.2),(max(ST_swe_max, na.rm=T)*0.95), paste(year[y],"-",year[(y+1)]))
		if (y == 2 || y == 5 || y == 8  || y == 11) {
	  		axis(2, tck = 0.03)				
		} else  {
	  		axis(2, labels=NA, tck = 0.03)	
		}		
		if (y == 2) {
			legend('topright',c("Snowtel","CLM4.5","cum CLM snow"), col=c(1,2,2), lty=c(1,1,2), lwd=2, bty="n")
		}
		if (y > 10) {
	  		axis(1, labels=month, at=mon,tck = 0.03)				
		} else {		
			axis(1, labels=NA, at=mon, tck = 0.03)	
		}
		if (y == 5) {mtext("SWE (mm)", side = 2, line = 3, adj = 0)}
		start <- end + 1
	}
dev.next()

	start <- 8 + 12
	for (y in 2:(nyear-1)) {
		print(year[y])
		end <- start + 11
		plot( mon,ST_depth_max[start:end], type="l", lwd=2, ylim=c(0,max(ST_depth_max, na.rm=T)),
		  xaxt = "n", yaxt = "n", xlab=NA, ylab="Max Snow depth (mm)")
   	    lines(mon,1e3*CLM_depth_max[start:end], col=2, lwd=2) #convert from m to mm
   	    text((mon[1]*1.2),(0.95*max(ST_depth_max, na.rm=T)), paste(year[y],"-",year[(y+1)]))
		if (y == 2 || y == 5 || y == 8  || y == 11) {
	  		axis(2, tck = 0.03)				
		} else  {
	  		axis(2, labels=NA, tck = 0.03)	
		}		
		if (y == 2) {
			legend('topright',c("Snowtel","CLM4.5"), col=c(1,2), lty=1, lwd=2, bty="n")
		}
		if (y > 10) {
	  		axis(1, labels=month, at=mon,tck = 0.03)				
		} else {		
			axis(1, labels=NA, at=mon, tck = 0.03)	
		}
		if (y == 5) {mtext("Max Snow Depth (mm)", side = 2, line = 3, adj = 0)}
		start <- end + 1
	}
dev.off()


#---------------------------------------------------------
#----------------  END  ----------------------------------
#---------------------------------------------------------


  # dates do not seem to align, likely off by one timestep
  cor.test(NR1.tsa, Data.sno$ST_Tavg)
  cor.test(NR1.ppt, Data.sno$ST_ppt)

  test.c <- cbind(NR.ppt, Data.sno$ST_ppt)
  test.c[1:100,] 
  ST.ppt            <- rep(NA, nday)
  ST.Tavg           <- rep(NA, nday)
  ST.swe            <- rep(NA, nday)
  ST.depth          <- rep(NA, nday)
  ST.col            <- rep(NA, nday)
  ST.ppt[1:nday-1]  <- Data.sno$ST_ppt[2:nday]
  ST.swe[1:nday-1]  <- Data.sno$ST_swe[2:nday]
  ST.depth[1:nday-1]<- Data.sno$ST_depth[2:nday]
  ST.Tavg[1:nday-1] <- Data.sno$ST_Tavg[2:nday]
  ST.col <- ifelse(ST.Tavg <= 0, 4,2)

  par(mfrow=c(2,2),mar=c(5,5,2,2))
  test <- cor.test(NR1.tsa, TSA_mean)
  r1 <- paste("r =",signif (test$estimate,digits=2))
  plot(NR1.tsa,TSA_mean, xlim=c(-30,20),ylim=c(-30,20), 
	xlab="NR1 mean daily 2m air temp (˚C)", ylab="CLM mean daily air temp (˚C)")
	abline(0,1, lty=2)
	text(15,-28,r1, adj = c(0.5,0))

  test <- cor.test(NR1.tsa, Data.sno$ST_Tavg)
  r1 <- paste("r =",signif (test$estimate,digits=2))
  plot(NR1.tsa,Data.sno$ST_Tavg,xlim=c(-30,20),ylim=c(-30,20),  
	xlab="NR1 mean daily 2m air temp (˚C)", ylab="SNOTEL air temp (˚C)")
	abline(0,1, lty=2)
	text(15,-28,r1, adj = c(0.5,0))


  test <- cor.test(NR1.tsa, ST.Tavg)
  r1 <- paste("r =",signif (test$estimate,digits=2))
  plot(NR1.tsa,ST.Tavg,xlim=c(-30,20),ylim=c(-30,20),  
	xlab="NR1 mean daily 2m air temp (˚C)", ylab="SNOTEL* air temp (˚C)")
	abline(0,1, lty=2)
	text(15,-28,r1, adj = c(0.5,0))

  test <- cor.test(TSA_mean, ST.Tavg)
  r1 <- paste("r =",signif (test$estimate,digits=2))
  plot(TSA_mean,ST.Tavg,xlim=c(-30,20),ylim=c(-30,20),  
	xlab="CLM mean daily 2m air temp (˚C)", ylab="SNOTEL* air temp (˚C)")
	abline(0,1, lty=2)
	text(15,-28,r1, adj = c(0.5,0))




  plot(NR.ppt,Data.sno$ST_ppt, ylim=c(0,60), xlim=c(0,60),
	xlab="NR1 precip (mm/d)", ylab="SNOTEL precip (mm/d)")
	abline(0,1, lty=2)


  plot(NR.Tavg_2m,Data.sno$ST_Tavg, 
	xlab="NR1 mean daily 2m air temp (˚C)", ylab="SNOTEL air temp (˚C)")
	abline(0,1, lty=2)
  plot(NR.ppt,Data.sno$ST_ppt, ylim=c(0,60), xlim=c(0,60),
	xlab="NR1 precip (mm/d)", ylab="SNOTEL precip (mm/d)")
	abline(0,1, lty=2)





#-----------------------------------------------------
# calcuate annual totals for tower data (NR) & SNOTEL (ST)
#-----------------------------------------------------
  ST.year     <- Data.sno$Year
  NR_ann_snow <- rep(NA, nyear)
  NR_ann_rain <- rep(NA, nyear)
  NR_ann_ppt  <- rep(NA, nyear)
  ST_ann_snow <- rep(NA, nyear)
  ST_ann_rain <- rep(NA, nyear)
  ST_ann_ppt  <- rep(NA, nyear)
  for (y in 1:nyear) {
  	NR_ann_snow[y] <- sum(NR.ppt[NR.year==year[y] & NR.Tavg_2m <= 0], na.rm=T)
  	NR_ann_rain[y] <- sum(NR.ppt[NR.year==year[y] & NR.Tavg_2m >  0], na.rm=T)
  	NR_ann_ppt[y]  <- sum(NR.ppt[NR.year==year[y]], na.rm=T )

  	ST_ann_snow[y] <- sum(ST.ppt[ST.year==year[y] & ST.Tavg <= 0], na.rm=T)
  	ST_ann_rain[y] <- sum(ST.ppt[ST.year==year[y] & ST.Tavg >  0], na.rm=T)
  	ST_ann_ppt[y]  <- sum(ST.ppt[ST.year==year[y]], na.rm=T )
  }
ST_ann_snow + ST_ann_rain
ST_ann_ppt
NR_ann_snow + NR_ann_rain
NR_ann_ppt

NR_ppt <-paste("  Tower  ppt =",signif (mean(NR_ann_ppt),digits=3),"+/-",signif (sd(NR_ann_ppt),digits=3),"mm/y")
ST_ppt <-paste("SNOTEL ppt =", signif (mean(ST_ann_ppt),digits=3),"+/-",signif (sd(ST_ann_ppt),digits=3),"mm/y")
t.test(NR_ann_ppt, ST_ann_ppt, paired=T)
NR_ppt 
ST_ppt

NR_snow <-paste("NR snow =",signif (mean(NR_ann_snow),digits=3),"+/-",signif(sd(NR_ann_snow),digits=3),"mm/y") 
ST_snow <-paste("ST snow =",signif (mean(ST_ann_snow),digits=3),"+/-",signif(sd(ST_ann_snow),digits=3),"mm/y")
t.test(NR_ann_snow, ST_ann_snow, paired=T)
NR_snow 
ST_snow

NR_rain<-paste("Tower  rain =",signif(mean(NR_ann_rain),digits=3),"+/-",signif(sd(NR_ann_rain),digits=3),"mm/y") 
ST_rain<-paste("SNOTEL rain =",signif(mean(ST_ann_rain),digits=3),"+/-",signif(sd(ST_ann_rain),digits=3),"mm/y")
t.test(NR_ann_rain, ST_ann_rain, paired=T)
NR_rain 
ST_rain

plot(NR_ann_ppt,ST_ann_ppt, ylim=c(450,1000), xlim=c(450,1000), pch=16,
    xlab="NR1 precip (mm/y)", ylab="SNOTEL precip (mm/y)", cex=1.2,
    main="US-NR1 1999-2007")
	abline(1,1, lty=2)
text(820,465, NR_ppt, cex=0.8)
text(820,440, ST_ppt, cex=0.8)

lab    <- ("     Rain   Snow")
tower  <- paste("US-NR1   ",signif(mean(NR_ann_rain),digits=3),"  ",signif(mean(NR_ann_snow),digits=3))
snotel <- paste("SNOTEL  ", signif(mean(ST_ann_rain),digits=3),"  ",signif(mean(ST_ann_snow),digits=3))

plot(NR_ann_rain,ST_ann_rain, ylim=c(150,600), xlim=c(150,600), cex=1.2,
    xlab="NR1 precip (mm/y)", ylab="SNOTEL precip (mm/y)",pch=16,col=2,
    main="US-NR1 1999-2007")
points(NR_ann_snow,ST_ann_snow,pch=16,col=4, cex=1.2)
legend('topleft',c("rain","snow"), col=c(2,4), pch=c(16,16), bty="n", cex=1.2)
    abline(1,1, lty=2)
text(520,210, lab, cex=0.8)
text(480,185, tower, cex=0.8)
text(480,160, snotel, cex=0.8)

ST_ann_snow + ST_ann_rain
	regu_days <- c(31,28,31,30,31,30,31,31,30,31,30,31)
	leap_days <- c(31,29,31,30,31,30,31,31,30,31,30,31)
    regu_steps <- regu_days * 48
    leap_steps <- leap_days * 48
    
	nyear <- length(year)
  Data.ann <- Data.F[Data.F$Year==2008, ]  
 
  names(Data.temp)

LAT <- 
LON <-   
#-------------------------------------------------------------
#---------------write out .nc file----------------------------
#-------------------------------------------------------------
# define the netcdf coordinate variables (name, units, type)
lat     <- dim.def.ncdf("lat","degrees_north", as.double(LAT))
lon     <- dim.def.ncdf("lon","degrees_east", as.double(LON))

# Make varables (name, units, dims, missing_value)
mv <- -9999             # missing value to use

for (y in 1:nyear) {
  Data.ann <- Data.F[Data.F$Year==year[y], ]  
  sstep <- 1	

  if(year[y]==2008 || year[y]==2012) {
  	nsteps <- leap_steps
  } else {
  	nsteps <- regu_steps
  }

  for (m in 1:12) {
	estep <- sstep + nsteps[m] - 1
     Data.mon <- Data.ann[sstep:estep, ]
     names(Data.mon)
	time  <- Data.mon$DateTime
	Rh    <- var.def.ncdf("Rh", Data.F$Rh, list(time), mv)

	print(paste(year[y],m,sstep,estep))
	
	sstep <- estep + 1

  }                        # close monthly loop
  remove(Data.ann, nsteps, sstep, estep)

}                          # close annual loop 

  
