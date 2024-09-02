
import sys
import os
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import iris
import iris.plot as iplt
import iris.quickplot as qplt
import iris.coord_categorisation
import iris.analysis.cartography
import matplotlib.dates as mdates
import cartopy.crs as ccrs
import pylab as pl
import matplotlib.pyplot as plt
from iris.analysis.cartography import area_weights
import iris.analysis.cartography
import ascend
from ascend import shape
import cartopy.feature as cfeature
from collections import OrderedDict 
from iris.coord_systems import GeogCS
import h5py
import cf_units
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from iris.coords import DimCoord
from matplotlib.colors import ListedColormap, BoundaryNorm
import cftime


### First convert GFED HDF files to netCDF
'''
def burned_fraction(filename, year):
    with h5py.File(filename, 'r') as obj:
         cs = iris.coord_systems.GeogCS(6371229)
         lat_coord = iris.coords.DimCoord(obj['lat'][:, 0],
                                     standard_name='latitude',
                                     units='degrees',
                                     coord_system=cs)

         lon_coord = iris.coords.DimCoord(obj['lon'][0, :],
                                     standard_name='longitude',
                                     units='degrees',
                                     coord_system=cs)

        epoch = 'days since 2001-01-01'
        months = obj['burned_area'].keys()
        time_coord = iris.coords.DimCoord([cftime.date2num(cftime.datetime(year, int(m), 1), 
                                                           epoch) 
                                           for m in months], 
                                           standard_name='time',
                                           units=epoch)
        
        data = np.stack([obj['burned_area'][m]['burned_fraction'][:] for m in months])

        return iris.cube.Cube(
            data,
            long_name='burned_fraction',
            dim_coords_and_dims=[(time_coord, 0), 
                                 (lat_coord, 1),
                                 (lon_coord, 2)])



### load data for 1997-2016 and convert to netcdf & save
years = np.arange(1997, 2017)
for y in years:
    cube = burned_fraction(f'/data/cr1/cburton/GFED/GFED4.1s_{y}.hdf5', y)
    iris.save(cube.collapsed('time', iris.analysis.SUM), f'/scratch/cburton/GFED/bf_{y}.nc')
    print(y)

### load & save the data for 1997-2016 as a single cube
bf_years = iris.load_cube([f'/scratch/cburton/GFED/bf_{y}.nc' for y in years])
bf_years.coord('latitude').guess_bounds()
bf_years.coord('longitude').guess_bounds()
iris.save(bf_years, '/data/cr1/cburton/GFED/GFED4s_AnnualTotalBA_1997-2016.nc')
'''    


### Second read in data and plot maps
date = iris.Constraint(time=lambda cell: 2001 <= cell.point.year <= 2014)

def PrepareData(model, obs):
    cs_new = iris.coord_systems.GeogCS(6371229.)
    obs.coord('latitude').coord_system = cs_new
    obs.coord('longitude').coord_system = cs_new
    model.coord('latitude').coord_system = cs_new
    model.coord('longitude').coord_system = cs_new
    obs = obs.regrid(model, iris.analysis.Linear())
    obs.units = cf_units.Unit('%')
    model.units = cf_units.Unit('%')
    obs.rename('burntFractionAll')
    model.rename('burntFractionAll')
    return model, obs

def CollapseToTimeseries(cube):
    coords = ('longitude', 'latitude')
    for coord in coords:
        if not cube.coord(coord).has_bounds():
            cube.coord(coord).guess_bounds()
    area = iris.analysis.cartography.area_weights(cube, normalize=False)
    cube = cube.collapsed(coords, iris.analysis.SUM, weights=area)/1E10
    return cube



### UKESM ###
folder = '/scratch/cburton/scratch/OptimESM/'
var_constraint = iris.Constraint(name="m01s08i700")

#With evolving ice sheets
ens1 = iris.load_cube(folder+'cy623a.py*.pp', var_constraint)
ens2 = iris.load_cube(folder+'da914a.py*.pp', var_constraint)
ens3 = iris.load_cube(folder+'da916a.py*.pp', var_constraint)
ens4 = iris.load_cube(folder+'da917a.py*.pp', var_constraint)
IceMean = (ens1+ens2+ens3+ens4)/4

UKESM = IceMean.extract(date)*86400*365 


var_constraint = iris.Constraint(name="burntFractionAll")
### CNRM ###
CNRM = iris.load_cube(folder+'burntFractionAll*CNRM-ESM2-1*.nc', var_constraint)
CNRM = CNRM.extract(date)   
iris.coord_categorisation.add_year(CNRM, 'time', name='year')
CNRM = CNRM.aggregated_by(['year'],iris.analysis.SUM)/100


### ECEarth ###
ECEarth2 = iris.load(folder+'burntFractionAll*EC-Earth3*r3i1p1f1*.nc', var_constraint)
for cube in ECEarth2:
    cube.attributes = None
ECEarth2 = ECEarth2.concatenate_cube()
ECEarth3 = iris.load(folder+'burntFractionAll*EC-Earth3*r5i1p1f1*.nc', var_constraint)
for cube in ECEarth3:
    cube.attributes = None
ECEarth3 = ECEarth3.concatenate_cube()

ECEarth = (ECEarth2 + ECEarth3)/2
ECEarth = ECEarth.extract(date)  
iris.coord_categorisation.add_year(ECEarth, 'time', name='year')
ECEarth = ECEarth.aggregated_by(['year'],iris.analysis.SUM)/100


### GFED4 ###
GFED4 = iris.load_cube('/data/cr1/cburton/GFED/GFED4s_AnnualTotalBA_1997-2016.nc')
GFED4 = GFED4.extract(date)   
UKESM,GFED4 = PrepareData(UKESM,GFED4) 

### GFED5 ###
GFED5 = iris.load_cube('/data/cr1/cburton/GFED/GFED5/GFED5_Burned_Percentage.nc')
GFED5 = GFED5.extract(date)   
iris.coord_categorisation.add_year(GFED5, 'time', name='year')
GFED5 = GFED5.aggregated_by(['year'],iris.analysis.SUM)/100
UKESM,GFED5 = PrepareData(UKESM,GFED5) 

### FireCCI ###
CCI = iris.load_cube('/data/cr1/cburton/FireCCI/2001-2016_Annual_BAFrac_CB.nc')
CCI = CCI.extract(date)
UKESM,CCI = PrepareData(UKESM,CCI) 

CNRM,GFED4 = PrepareData(CNRM,GFED4)
CNRM,GFED5 = PrepareData(CNRM,GFED5)
ECEarth,GFED4 = PrepareData(ECEarth,GFED4)
ECEarth,GFED5 = PrepareData(ECEarth,GFED5)


UKESM = CollapseToTimeseries(UKESM)
CNRM = CollapseToTimeseries(CNRM)
ECEarth = CollapseToTimeseries(ECEarth)
GFED4 = CollapseToTimeseries(GFED4)
GFED5 = CollapseToTimeseries(GFED5)
CCI = CollapseToTimeseries(CCI)
    

############### Make plot ###################

x = np.arange(2001,2015)
plt.plot(x, GFED4.data, color='black', label='GFED4.1s')
plt.plot(x, GFED5.data, color='darkblue', label='GFED5')
plt.plot(x, CCI.data, color='blue', label='FireCCI')
plt.plot(x, UKESM.data, color='Red', label='UKESM')
plt.plot(x, CNRM.data, color='Green', label='CNRM')
plt.plot(x, ECEarth.data, color='Orange', label='EC-Earth')


plt.ylabel('Mha yr$^{-1}$')
plt.xlabel('Years')
plt.legend(loc='best', ncol=2)
plt.title('Annual global total burned area')
plt.show()
 


