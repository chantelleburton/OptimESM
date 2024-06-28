
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

### UKESM ###
folder = '/scratch/cburton/scratch/OptimESM/'
var_constraint = iris.Constraint(name="m01s08i701")

#With evolving ice sheets
ens1 = iris.load_cube(folder+'cy623a.py*.pp', var_constraint)
ens2 = iris.load_cube(folder+'da914a.py*.pp', var_constraint)
ens3 = iris.load_cube(folder+'da916a.py*.pp', var_constraint)
ens4 = iris.load_cube(folder+'da917a.py*.pp', var_constraint)
IceMean = (ens1+ens2+ens3+ens4)/4

#Without evolving ice sheets
ensA = iris.load_cube(folder+'cy690a.py*.pp', var_constraint)
ensB = iris.load_cube(folder+'cy691a.py*.pp', var_constraint)
ensC = iris.load_cube(folder+'cy692a.py*.pp', var_constraint)
ensD = iris.load_cube(folder+'cy693a.py*.pp', var_constraint)
NoIceMean = (ensA+ensB+ensC+ensD)/4

coords = ('longitude', 'latitude')
for coord in coords:
    if not IceMean.coord(coord).has_bounds():
        IceMean.coord(coord).guess_bounds()
area = iris.analysis.cartography.area_weights(IceMean, normalize=False)

date = iris.Constraint(time=lambda cell: 2001 <= cell.point.year <= 2014)

### GFED4 ###
GFED4 = iris.load_cube('/data/cr1/cburton/GFED/GFED4s_AnnualTotalEmissions_1997-2016.nc')
GFED4 = GFED4.extract(date)   
GFED4 = GFED4.regrid(IceMean, iris.analysis.Linear())
GFED4 = GFED4.collapsed(coords, iris.analysis.SUM, weights=area)/1E15


### Get GFAS data
var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'cfire')
GFAS = iris.load_cube('/data/cr1/cburton/GFAS/GFAS.nc', var_constraint)*86400*360 
cs_new = iris.coord_systems.GeogCS(6371229.)
GFAS.coord('latitude').coord_system = cs_new
GFAS.coord('longitude').coord_system = cs_new
IceMean.coord('latitude').coord_system = cs_new
IceMean.coord('longitude').coord_system = cs_new
GFAS = GFAS.regrid(IceMean, iris.analysis.Linear())
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


IceMean = IceMean.collapsed(coords, iris.analysis.SUM, weights=area)*86400*365/1E12
NoIceMean = NoIceMean.collapsed(coords, iris.analysis.SUM, weights=area)*86400*365/1E12


###Make plot
x = np.arange(2001,2015)
plt.plot(x, GFED4.data, color='black', label='GFED4.1s')
plt.plot(xa, GFAS, color='blue', label='GFAS')
plt.plot(x, IceMean.data, color='Red', label='UKESM')
plt.ylabel('GtC')
plt.xlabel('Years')
plt.legend(loc='best')
plt.title('Annual global total fire emissions')
plt.show()
 





