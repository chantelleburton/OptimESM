
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

### UKESM ###
folder = '/scratch/cburton/scratch/OptimESM/'
var_constraint = iris.Constraint(name="m01s08i700")

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

IceMean = IceMean.collapsed('time', iris.analysis.MEAN)*86400*365*100
NoIceMean = NoIceMean.collapsed('time', iris.analysis.MEAN)*86400*365*100

date = iris.Constraint(time=lambda cell: 2001 <= cell.point.year <= 2014)

### GFED4 ###
GFED4 = iris.load_cube('/data/cr1/cburton/GFED/GFED4s_AnnualTotalBA_1997-2016.nc')
GFED4 = GFED4.extract(date)   
GFED4 = GFED4.collapsed('time', iris.analysis.MEAN)*100
GFED4 = GFED4.regrid(IceMean, iris.analysis.Linear())

### GFED5 ###
GFED5 = iris.load_cube('/data/cr1/cburton/GFED/GFED5/GFED5_Burned_Percentage.nc')
GFED5 = GFED5.extract(date)   
iris.coord_categorisation.add_year(GFED5, 'time', name='year')
GFED5 = GFED5.aggregated_by(['year'],iris.analysis.SUM)
GFED5 = GFED5.collapsed('time', iris.analysis.MEAN)
GFED5 = GFED5.regrid(IceMean, iris.analysis.Linear())

### FireCCI ###
CCI = iris.load_cube('/data/cr1/cburton/FireCCI/2001-2016_Annual_BAFrac_CB.nc')
CCI = CCI.extract(date)   
CCI = CCI.collapsed('time', iris.analysis.MEAN)*100
cs_new = iris.coord_systems.GeogCS(6371229.)
CCI.coord('latitude').coord_system = cs_new
CCI.coord('longitude').coord_system = cs_new
IceMean.coord('latitude').coord_system = cs_new
IceMean.coord('longitude').coord_system = cs_new
CCI = CCI.regrid(IceMean, iris.analysis.Linear())

#Bias
diff4 = IceMean - GFED4
diff5 = GFED5.copy()
diff5.data = IceMean.data - GFED5.data

crs_latlon = ccrs.PlateCarree()#
def make_plot(projection_name, projection_crs):

###############
    fig = plt.figure()

    jet = plt.get_cmap('jet')
    colors = jet(np.linspace(0.1, 1, 256))
    colors = np.vstack(([1, 1, 1, 1], colors))
    cmap = ListedColormap(colors)
    tick_values = [0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
    levels = MaxNLocator(nbins=10).tick_values=tick_values
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
  
    ax = plt.subplot(2,3,1)
    iplt.pcolormesh(CCI, cmap=cmap, norm=norm )
    plt.gca().coastlines()
    plt.title("FireCCI", fontsize=12)
    plt.colorbar(orientation='horizontal')

    ax = plt.subplot(2,3,2)
    iplt.pcolormesh(GFED4, cmap=cmap, norm=norm )
    plt.gca().coastlines()
    plt.title("GFED4.1s", fontsize=12)
    plt.colorbar(orientation='horizontal')

    ax = plt.subplot(2,3,3)
    iplt.pcolormesh(GFED5, cmap=cmap, norm=norm )
    plt.gca().coastlines()
    plt.title("GFED5", fontsize=12)
    plt.colorbar(orientation='horizontal')
    
    ax = plt.subplot(2,3,4)
    iplt.pcolormesh(IceMean, cmap=cmap, norm=norm )
    plt.gca().coastlines()
    plt.title("UKESM", fontsize=12)
    plt.colorbar(orientation='horizontal', label='Burned Fraction (% yr$^{-1}$)')

    ax = plt.subplot(2,3,5)
    mesh = iplt.pcolormesh(diff4, cmap='RdBu_r', vmin=-20, vmax=20 )
    plt.gca().coastlines()
    plt.title("Difference UKESM-GFED4.1s", fontsize=12)
    plt.colorbar(orientation='horizontal')

    ax = plt.subplot(2,3,6)
    mesh = iplt.pcolormesh(diff5, cmap='RdBu_r', vmin=-20, vmax=20 )
    plt.gca().coastlines()
    plt.title("Difference UKESM-GFED5", fontsize=12)
    plt.colorbar(orientation='horizontal')

    plt.suptitle('2001-2014 Burned Area')
    plt.show()
            
def main():
    make_plot('Plot', ccrs.PlateCarree())
    
if __name__ == '__main__':
    main()


