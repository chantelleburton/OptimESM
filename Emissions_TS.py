
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
def emitted_fraction(filename, year):
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
         months = obj['emissions'].keys()

         time_coord = iris.coords.DimCoord([cftime.date2num(cftime.datetime(year, int(m), 1), 
                                                           epoch) 
                                           for m in months], 
                                           standard_name='time',
                                           units=epoch)


         data = np.stack([obj['emissions'][m]['C'] for m in months])

         return iris.cube.Cube(
            data,
            long_name='emissions',
            dim_coords_and_dims=[(time_coord, 0), 
                                 (lat_coord, 1),
                                 (lon_coord, 2)])



### load data for 1997-2016 and convert to netcdf & save
years = np.arange(1997, 2017)
for y in years:
    cube = emitted_fraction(f'/data/cr1/cburton/GFED/GFED4.1s_{y}.hdf5', y)
    iris.save(cube.collapsed('time', iris.analysis.SUM), f'/scratch/cburton/GFED/ef_{y}.nc')
    print(y)

### load & save the data for 1997-2016 as a single cube
bf_years = iris.load_cube([f'/scratch/cburton/GFED/ef_{y}.nc' for y in years])
bf_years.coord('latitude').guess_bounds()
bf_years.coord('longitude').guess_bounds()
iris.save(bf_years, '/data/cr1/cburton/GFED/GFED4s_AnnualTotalEmissions_1997-2016.nc')
exit()    
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
    obs.units = cf_units.Unit('kg/m2/s')
    model.units = cf_units.Unit('kg/m2/s')
    obs.rename('fFire')
    model.rename('fFire')
    return model, obs

def CollapseToTimeseries(cube):
    coords = ('longitude', 'latitude')
    for coord in coords:
        if not cube.coord(coord).has_bounds():
            cube.coord(coord).guess_bounds()
    area = iris.analysis.cartography.area_weights(cube, normalize=False)
    cube = cube.collapsed(coords, iris.analysis.SUM, weights=area)/1E12
    return cube

### UKESM ###
folder = '/scratch/cburton/scratch/OptimESM/'
var_constraint = iris.Constraint(name="m01s08i701")

#With evolving ice sheets
ens1 = iris.load_cube(folder+'cy623a.py*.pp', var_constraint)
ens2 = iris.load_cube(folder+'da914a.py*.pp', var_constraint)
ens3 = iris.load_cube(folder+'da916a.py*.pp', var_constraint)
ens4 = iris.load_cube(folder+'da917a.py*.pp', var_constraint)
UKESM = ((ens1+ens2+ens3+ens4)/4)*86400*365

var_constraint = iris.Constraint(name="fFire")
### CNRM ###
CNRM = iris.load_cube(folder+'fFire*CNRM-ESM2-1*.nc', var_constraint)
CNRM = CNRM.extract(date)   
iris.coord_categorisation.add_year(CNRM, 'time', name='year')
CNRM = CNRM.aggregated_by(['year'],iris.analysis.SUM)*1E6
LandFrac = iris.load_cube('/scratch/cburton/scratch/OptimESM/sftlf_fx_CNRM-ESM2-1_esm-hist_r1i1p2f2_gr.nc', 'sftlf')/100
CNRM = CNRM*LandFrac
CNRM,UKESM = PrepareData(CNRM,UKESM)
CNRM = CollapseToTimeseries(CNRM)


### ECEarth ###
ECEarth1 = iris.load(folder+'fFire*EC-Earth3*r1i1p1f1*.nc', var_constraint)
for cube in ECEarth1:
    cube.attributes = None
ECEarth1 = ECEarth1.concatenate_cube()
ECEarth1 = ECEarth1.extract(date)  

ECEarth2 = iris.load(folder+'fFire*EC-Earth3*r3i1p1f1*.nc', var_constraint)
for cube in ECEarth2:
    cube.attributes = None
ECEarth2 = ECEarth2.concatenate_cube()
ECEarth2 = ECEarth2.extract(date)  

ECEarth3 = iris.load(folder+'fFire*EC-Earth3*r5i1p1f1*.nc', var_constraint)
for cube in ECEarth3:
    cube.attributes = None
ECEarth3 = ECEarth3.concatenate_cube()
ECEarth3 = ECEarth3.extract(date)  

ECEarth = (ECEarth1 + ECEarth2 + ECEarth3)/3
iris.coord_categorisation.add_year(ECEarth, 'time', name='year')
ECEarth = ECEarth.aggregated_by(['year'],iris.analysis.SUM)*1E6
LandFrac = iris.load_cube('/scratch/cburton/scratch/OptimESM/sftlf_fx_EC-Earth3-ESM-1_esm-hist_r5i1p1f1_gr.nc', 'sftlf')/100
ECEarth.data = ECEarth.data * LandFrac.data
ECEarth,UKESM = PrepareData(ECEarth,UKESM)
ECEarth = CollapseToTimeseries(ECEarth)

### GFED4 ###
GFED4 = iris.load_cube('/data/cr1/cburton/GFED/GFED4s_AnnualTotalEmissions_1997-2016.nc')
GFED4 = GFED4.extract(date)   
GFED4 = GFED4.regrid(UKESM, iris.analysis.Linear())
GFED4 =  CollapseToTimeseries(GFED4)/1000

### Get GFAS data
var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'cfire')
GFAS = iris.load_cube('/data/cr1/cburton/GFAS/GFAS.nc', var_constraint)*86400*365 
GFAS,UKESM = PrepareData(GFAS,UKESM)
years = range(2000, 2014)
F = []
for year in years:
    F.append(GFAS[year-2000:year-1999, :, :])

coords = ('longitude', 'latitude')
for f in F:
    for coord in coords:
        if not f.coord(coord).has_bounds():
            f.coord(coord).guess_bounds()

weights = iris.analysis.cartography.area_weights(F[0])
F_sum = []
for f in F:
    F_sum.append(f.collapsed(coords, iris.analysis.SUM, weights=weights) / 1E12)

xa = tuple(years)
GFAS = tuple(f.data for f in F_sum)

UKESM =  CollapseToTimeseries(UKESM)

###Make plot
x = np.arange(2001,2015)
plt.plot(x, GFED4.data, color='black', label='GFED4.1s')
plt.plot(xa, GFAS, color='blue', label='GFAS')
plt.plot(x, UKESM.data, color='Red', label='UKESM')
plt.plot(x, CNRM.data, color='Green', label='CNRM')
plt.plot(x, ECEarth.data, color='Orange', label='EC-Earth')

plt.ylabel('GtC')
plt.xlabel('Years')
plt.legend(loc='best')
plt.title('Annual global total fire emissions')
plt.show()
 





