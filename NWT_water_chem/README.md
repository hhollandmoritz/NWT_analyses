Will Wieder
Created April 2015
Modified May & Oct 2018
Plots DOC-NO3 relationship (See Taylor and Townsend, Nature 2010)
Looks for trends in pH, SO4.., and DOC from Nel Ca..ine's Niwot water chemistry data 

uses gamm function [generalized additive mixed modelling] 
in mgcv (& nlme) package

http://stackoverflow.com/questions/12623027/how-to-analyse-irregular-time-series-in-r
using an additive model to "decompose" the seasonal and trend components. 
As this is a regression-based approach you need to model the residuals as a time series process to
account for lack of independence in the residuals.
http://www.fromthebottomoftheheap.net/2011/06/12/additive-modelling-and-the-hadcrut3v-global-mean-temperature-series/

## code here is really messy 

It could like be made more efficient & user friendly
switching to tidyR and ggplot may be helpful...

