Niwot Modeling README

This is a directory to keep track of  progress for using NWT data in CLM simulations.
The broad goal of the work is to develop scripts and workflow to  more easily update CLM simulations with new NWT data.

Here's a space wehre we can keep track of what has been done and needs to be addressed.

__Get familiar with CLM__
- [ ] Set up single point case at Harvard forest
- [ ] Run with CLM5 forcings, 
- [ ] Run CLM5 with NEON data
- [ ] Begin looking at output from SP and BGC simulations

__Forcing data__ 
- [ ] Run CLM5 w/ default (GSWP3  forcing for NWT grid)
- [ ] Run CLM5 w/ Tvan  observations from 2008-2014, already formatted for  CLM.
- [ ] Run CLM5 w/ NEON NIWO met data data (2018-2019 only)

*Additional development that  we'd  need to develop to continue running w/ ongoing Tvan measurements*   
- [ ] Revisit and update data Tvan E & W tower data
- [ ] Gapfill Tvan observations with data from other tower (when necessary)
- [ ] Develop pipeline to bring in other measurements (Saddle precip, US-NR1 solar radiation, etc.)

__Observational data for surface dataset__
*These could be NEON or LTER data products*
- [ ] set vegetation to Arctic C3 grass, remove land use timeseries
- [ ] modify soil texture / porosity (see Wieder et al. 2017 from Brooks et al. 1980)
- [ ] modify soil depth for different vegetation communities (see Wieder et al. 2017 from Brooks et al. 1980)

__More ideas__ 
- [x] Try simulations with the hillslope model
- [ ] modify LAI for SP simulations?

__Observational data for parameter file__
*These could be NEON or LTER data products*
- [ ] modify allocation parameters (see Wieder et al. 2017)
- [ ] modify plant traits (leafCN, SLA, Medlyn slope).

__Observational data for model evaluation__
- [ ] Flux tower observations
- [ ] Saddle snow grid
- [ ] Saddle productivity
- [ ] Sensor network
 

