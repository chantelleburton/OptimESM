
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

### UKESM ###
folder = '/scratch/cburton/scratch/OptimESM/'
var_constraint = iris.Constraint(name="m01s08i700")

#With evolving ice sheets
ens1 = iris.load_cube(folder+'cy623a.py*.pp', var_constraint)
ens2 = iris.load_cube(folder+'da914a.py*.pp', var_constraint)
ens3 = iris.load_cube(folder+'da916a.py*.pp', var_constraint)
ens4 = iris.load_cube(folder+'da917a.py*.pp', var_constraint)
IceMean = (ens1+ens2+ens3+ens4)/4

#Without evolving ice sheets (not currently used, just here for reference)
ensA = iris.load_cube(folder+'cy690a.py*.pp', var_constraint)
ensB = iris.load_cube(folder+'cy691a.py*.pp', var_constraint)
ensC = iris.load_cube(folder+'cy692a.py*.pp', var_constraint)
ensD = iris.load_cube(folder+'cy693a.py*.pp', var_constraint)
NoIceMean = (ensA+ensB+ensC+ensD)/4

#Multiply by sec and days in a year, and 100 to convert frac to %. Find mean over 2001-2014
UKESM = IceMean.extract(date)*86400*365*100 
UKESM = UKESM.collapsed('time', iris.analysis.MEAN)  
NoIceMean = NoIceMean.extract(date)*86400*365*100   


var_constraint = iris.Constraint(name="burntFractionAll")
### CNRM ###
#Convert daily % to monthly % (as advised by Roland), and aggregate to get annual total. Find mean over 2001-2014
CNRM = iris.load_cube(folder+'burntFractionAll*CNRM-ESM2-1*.nc', var_constraint)
CNRM = CNRM.extract(date)   
iris.coord_categorisation.add_year(CNRM, 'time', name='year')
CNRM = CNRM.aggregated_by(['year'],iris.analysis.SUM)*365/12
CNRM = CNRM.collapsed('time', iris.analysis.MEAN)  

### ECEarth ###
# Aggregate monthly data to get annual total. Find mean over 2001-2014
# Only using 2 ensemble members for now. r1i1p1f1 is missing 2010-2013, and data not yet available from Lund
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
ECEarth = ECEarth.aggregated_by(['year'],iris.analysis.SUM)
ECEarth = ECEarth.collapsed('time', iris.analysis.MEAN)  

### GFED4 ###
# In fraction so multiply by 100 to convert to %
GFED4 = iris.load_cube('/data/cr1/cburton/GFED/GFED4s_AnnualTotalBA_1997-2016.nc')
GFED4 = GFED4.extract(date)*100   
GFED4 = GFED4.collapsed('time', iris.analysis.MEAN)  
UKESM,GFED4 = PrepareData(UKESM,GFED4) 

### GFED5 ###
# Already in %
GFED5 = iris.load_cube('/data/cr1/cburton/GFED/GFED5/GFED5_Burned_Percentage.nc')
GFED5 = GFED5.extract(date)   
iris.coord_categorisation.add_year(GFED5, 'time', name='year')
GFED5 = GFED5.aggregated_by(['year'],iris.analysis.SUM)
GFED5 = GFED5.collapsed('time', iris.analysis.MEAN)  
UKESM,GFED5 = PrepareData(UKESM,GFED5) 

### FireCCI ###
# In fraction so multiply by 100 to convert to %
CCI = iris.load_cube('/data/cr1/cburton/FireCCI/2001-2016_Annual_BAFrac_CB.nc')
CCI = CCI.extract(date)*100
CCI = CCI.collapsed('time', iris.analysis.MEAN)  
UKESM,CCI = PrepareData(UKESM,CCI) 

#Bias
diff1 = UKESM - GFED4
diff2 = GFED5.copy()
diff2.data = UKESM.data - GFED5.data

CNRM,GFED4 = PrepareData(CNRM,GFED4)
diff3 = CNRM - GFED4
CNRM,GFED5 = PrepareData(CNRM,GFED5)
diff4 = GFED5.copy()
diff4.data = CNRM.data - GFED5.data

ECEarth,GFED4 = PrepareData(ECEarth,GFED4)
diff5 = GFED4.copy()
diff5.data = ECEarth.data - GFED4.data
ECEarth,GFED5 = PrepareData(ECEarth,GFED5)
diff6 = GFED5.copy()
diff6.data = ECEarth.data - GFED5.data


############### Make plots ###################
crs_latlon = ccrs.PlateCarree()
def make_plot(projection_name, projection_crs):

    fig = plt.figure()

    jet = plt.get_cmap('jet')
    colors = jet(np.linspace(0.1, 1, 256))
    colors = np.vstack(([1, 1, 1, 1], colors))
    cmap = ListedColormap(colors)
    tick_values = [0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
    levels = MaxNLocator(nbins=10).tick_values=tick_values
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
  
    ax = plt.subplot(4,3,1)
    iplt.pcolormesh(CCI, cmap=cmap, norm=norm )
    plt.gca().coastlines()
    plt.title("FireCCI", fontsize=12)
    plt.colorbar(orientation='horizontal')

    ax = plt.subplot(4,3,2)
    iplt.pcolormesh(GFED4, cmap=cmap, norm=norm )
    plt.gca().coastlines()
    plt.title("GFED4.1s", fontsize=12)
    plt.colorbar(orientation='horizontal')

    ax = plt.subplot(4,3,3)
    iplt.pcolormesh(GFED5, cmap=cmap, norm=norm )
    plt.gca().coastlines()
    plt.title("GFED5", fontsize=12)
    plt.colorbar(orientation='horizontal')
    
    ax = plt.subplot(4,3,4)
    iplt.pcolormesh(UKESM, cmap=cmap, norm=norm )
    plt.gca().coastlines()
    plt.title("UKESM", fontsize=12)
    plt.colorbar(orientation='horizontal')

    ax = plt.subplot(4,3,5)
    mesh = iplt.pcolormesh(diff1, cmap='RdBu_r', vmin=-20, vmax=20 )
    plt.gca().coastlines()
    plt.title("Difference UKESM-GFED4.1s", fontsize=12)
    plt.colorbar(orientation='horizontal')

    ax = plt.subplot(4,3,6)
    mesh = iplt.pcolormesh(diff2, cmap='RdBu_r', vmin=-20, vmax=20 )
    plt.gca().coastlines()
    plt.title("Difference UKESM-GFED5", fontsize=12)
    plt.colorbar(orientation='horizontal')

    ax = plt.subplot(4,3,7)
    iplt.pcolormesh(CNRM, cmap=cmap, norm=norm )
    plt.gca().coastlines()
    plt.title("CNRM", fontsize=12)
    plt.colorbar(orientation='horizontal')

    ax = plt.subplot(4,3,8)
    mesh = iplt.pcolormesh(diff3, cmap='RdBu_r', vmin=-20, vmax=20 )
    plt.gca().coastlines()
    plt.title("Difference CNRM-GFED4.1s", fontsize=12)
    plt.colorbar(orientation='horizontal')

    ax = plt.subplot(4,3,9)
    mesh = iplt.pcolormesh(diff4, cmap='RdBu_r', vmin=-20, vmax=20 )
    plt.gca().coastlines()
    plt.title("Difference CNRM-GFED5", fontsize=12)
    plt.colorbar(orientation='horizontal')

    ax = plt.subplot(4,3,10)
    iplt.pcolormesh(ECEarth, cmap=cmap, norm=norm )
    plt.gca().coastlines()
    plt.title("EC-Earth", fontsize=12)
    plt.colorbar(orientation='horizontal', label='Burned Fraction (% yr$^{-1}$)')

    ax = plt.subplot(4,3,11)
    mesh = iplt.pcolormesh(diff5, cmap='RdBu_r', vmin=-20, vmax=20 )
    plt.gca().coastlines()
    plt.title("Difference EC-Earth-GFED4.1s", fontsize=12)
    plt.colorbar(orientation='horizontal', label='Difference in Burned Fraction (% yr$^{-1}$)')

    ax = plt.subplot(4,3,12)
    mesh = iplt.pcolormesh(diff6, cmap='RdBu_r', vmin=-20, vmax=20 )
    plt.gca().coastlines()
    plt.title("Difference EC-Earth-GFED5", fontsize=12)
    plt.colorbar(orientation='horizontal', label='Difference in Burned Fraction (% yr$^{-1}$)')

    plt.suptitle('2001-2014 Burned Area')
    plt.show()
            
def main():
    make_plot('Plot', ccrs.PlateCarree())
    
if __name__ == '__main__':
    main()


