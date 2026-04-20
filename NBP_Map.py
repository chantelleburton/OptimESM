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
from scipy import stats

date = iris.Constraint(time=lambda cell: 2000 <= cell.point.year <= 2014)

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
    cube = cube.collapsed(coords, iris.analysis.SUM, weights=area)/1E12
    return cube
    
def TimeMean(cube):
    cube = cube.collapsed('time', iris.analysis.MEAN)
    return cube     

### OBS ###
pwdin = "/data/users/eleanor.burke/obs_data/nbp"

cams_cube = iris.load_cube(pwdin + "/cams/cams73_latest_co2_flux_surface_mm.nc")
cams_cube = cams_cube.extract(date)
cams_cube = cams_cube.aggregated_by(['year'],iris.analysis.SUM)
cams = TimeMean(cams_cube) *86400*365/12*100  #(kg/m2/sec to g/m2/month))

carboscope_cube = iris.load_cube(pwdin + "/carboscope/r76nbetEXToc_v2025.flux_land.mon.nc")
carboscope = carboscope_cube.extract(date)
carboscope = carboscope.aggregated_by(['year'],iris.analysis.SUM)
carboscope = TimeMean(carboscope) *86400*365/12*100  #(kg/m2/sec to g/m2/month))

gcp2024 = iris.load_cube(pwdin + "/ct2022/GCP2024.flux1x1-monthly.processed.nc")
gcp2024 = gcp2024.extract(date)
gcp2024 = gcp2024.aggregated_by(['year'],iris.analysis.SUM)
gcp2024 = TimeMean(gcp2024)   *86400*365/12*100  #(kg/m2/sec to g/m2/month))


### MODEL ###
DataFolder = '/data/scratch/chantelle.burton/OptimESM/NBP/'
var_constraint = iris.Constraint(name="nbp")

### UKESM ###
UKESM = iris.load_cube(DataFolder+'nbp_Lmon_UKESM1-2-LL_esm-hist_r1i1p1f1_gn_195001-201412.nc', var_constraint)
UKESM = UKESM.extract(date)

iris.coord_categorisation.add_season_year(UKESM, 'time', name='year')
UKESM = UKESM.aggregated_by(['year'],iris.analysis.SUM)
UKESM = UKESM

UKESM = TimeMean(UKESM)
UKESM = UKESM*86400*365/12*100  #(kg/m2/sec to g/m2/month))

### CNRM ### 
CNRM = iris.load_cube(DataFolder+'nbp_Lmon_CNRM-ESM2-1_esm-hist_r1i1p2f2_gr_185001-201412.nc', var_constraint)
CNRM = CNRM.extract(date)
LandFrac = iris.load_cube(DataFolder+'sftlf_fx_CNRM-ESM2-1_esm-hist_r1i1p2f2_gr.nc', 'sftlf')/100
CNRM = CNRM*LandFrac

iris.coord_categorisation.add_season_year(CNRM, 'time', name='year')
CNRM = CNRM.aggregated_by(['year'],iris.analysis.SUM)
CNRM = CNRM

CNRM = TimeMean(CNRM)
CNRM = CNRM*86400*365/12*100  #(kg/m2/sec to g/m2/month))

### IPSL ### 
IPSL = iris.load_cube(DataFolder+'nbp_Lmon_IPSL-CM6-ESMCO2_esm-hist_r1i1p3f1_gr_195001-201412.nc', var_constraint)
IPSL = IPSL.extract(date)
LandFrac = iris.load_cube(DataFolder+'sftlf_fx_CNRM-ESM2-1_esm-hist_r1i1p2f2_gr.nc', 'sftlf')/100
CNRM = CNRM*LandFrac
iris.coord_categorisation.add_season_year(IPSL, 'time', name='year')
IPSL = IPSL.aggregated_by(['year'],iris.analysis.SUM)
IPSL = IPSL

IPSL = TimeMean(IPSL)
IPSL = IPSL*86400*365/12*100  #(kg/m2/sec to g/m2/month))

### ECEARTH ### 
ECEARTH = iris.load(DataFolder+'nbp_LPmon_EC-Earth3-ESM-1_esm-hist_r1i1p1f1_gr_20*.nc', var_constraint)
iris.util.equalise_attributes(ECEARTH)# new version
ECEARTH = ECEARTH.concatenate_cube()
ECEARTH = ECEARTH.extract(date)
#LandFrac = iris.load_cube(DataFolder+'sftlf_fx_CNRM-ESM2-1_esm-hist_r1i1p2f2_gr.nc', 'sftlf')/100
#CNRM = CNRM*LandFrac
iris.coord_categorisation.add_season_year(ECEARTH, 'time', name='year')
ECEARTH = ECEARTH.aggregated_by(['year'],iris.analysis.SUM)
ECEARTH = ECEARTH
ECEARTH = TimeMean(ECEARTH)
ECEARTH = ECEARTH*86400*365/12*100  #(kg/m2/sec to g/m2/month))
LandFrac = LandFrac.regrid(ECEARTH, iris.analysis.Linear())
ECEARTH.data = np.ma.masked_where(LandFrac.data == 0, ECEARTH.data)
LandFrac = LandFrac.regrid(gcp2024, iris.analysis.Linear())
gcp2024.data = np.ma.masked_where(LandFrac.data == 0, gcp2024.data)


CNRM=CNRM.extract(iris.Constraint(latitude=lambda cell: (-20.0) < cell < (20.0))) 
UKESM=UKESM.extract(iris.Constraint(latitude=lambda cell: (-20.0) < cell < (20.0))) 
IPSL=IPSL.extract(iris.Constraint(latitude=lambda cell: (-20.0) < cell < (20.0))) 
ECEARTH=ECEARTH.extract(iris.Constraint(latitude=lambda cell: (-20.0) < cell < (20.0))) 
gcp2024=gcp2024.extract(iris.Constraint(latitude=lambda cell: (-20.0) < cell < (20.0))) 
cams=cams.extract(iris.Constraint(latitude=lambda cell: (-20.0) < cell < (20.0))) 
carboscope=carboscope.extract(iris.Constraint(latitude=lambda cell: (-20.0) < cell < (20.0))) 


max_val = 10
plt.subplot(2,4,1)
mesh = iplt.pcolormesh(ECEARTH, vmin=-max_val, vmax=max_val)
plt.gca().coastlines(linewidth=0.2)
plt.title('EC-Earth')
plt.subplot(2,4,2)
iplt.pcolormesh(IPSL, vmin=-max_val, vmax=max_val)
plt.gca().coastlines(linewidth=0.2)
plt.title('IPSL')
plt.subplot(2,4,3)
iplt.pcolormesh(UKESM, vmin=-max_val, vmax=max_val)
plt.gca().coastlines(linewidth=0.2)
plt.title('UKESM')
plt.subplot(2,4,4)
iplt.pcolormesh(CNRM, vmin=-max_val, vmax=max_val)
plt.gca().coastlines(linewidth=0.2)
plt.title('CNRM')

plt.subplot(2,4,5)
iplt.pcolormesh(cams, vmin=-max_val, vmax=max_val)
plt.gca().coastlines(linewidth=0.2)
plt.title('CAMS')
plt.subplot(2,4,6)
iplt.pcolormesh(carboscope, vmin=-max_val, vmax=max_val)
plt.gca().coastlines(linewidth=0.2)
plt.title('CarboScope')
plt.subplot(2,4,7)
iplt.pcolormesh(gcp2024, vmin=-max_val, vmax=max_val)
plt.gca().coastlines(linewidth=0.2)
plt.title('GCP2024')




                         #Position: Left, Bottom, width, height
colorbar_axes = plt.gcf().add_axes([0.3, 0.14, 0.45, 0.04])
colorbar = plt.colorbar(mesh, colorbar_axes, orientation='horizontal', label='g/m2/month') 
plt.suptitle('NBP')
plt.show()


exit()



xs = np.arange(0,16)
plt.plot(xs,CNRM.data, 'blue', label='CNRM')
plt.plot(xs,UKESM.data, 'green',label='UKESM')
plt.plot(xs,IPSL.data, 'yellow', label='CNRM')
plt.plot(xs,ECEARTH.data, 'red',label='EC-EARTH')
years = ('2000', '2001', '2002', '2003', '2004', '2005','2006','2007','2008','2009','2010','2011','2012','2013','2014')
x_pos = (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14)
plt.xticks(x_pos, years)
plt.legend()
plt.title('NBP')
plt.ylabel('PgC')
plt.legend()
plt.show()
exit()


xs = np.arange(0,180)

ys = CNRM.data
slope, intercept, r_value, p_value, std_err = stats.linregress(xs,ys)
print (slope,r_value,p_value)
plt.plot(xs, intercept + slope*xs, 'b')
print ("r-squared:", r_value**2)
plt.plot(xs, CNRM.data, 'blue', label='CNRM')

ys = UKESM.data
slope, intercept, r_value, p_value, std_err = stats.linregress(xs,ys)
print (slope,r_value,p_value)
plt.plot(xs, intercept + slope*xs, 'g')
print ("r-squared:", r_value**2)
plt.plot(xs, UKESM.data, 'green', label='UKESM')

ys = IPSL.data
slope, intercept, r_value, p_value, std_err = stats.linregress(xs,ys)
print (slope,r_value,p_value)
plt.plot(xs, intercept + slope*xs, 'orange')
print ("r-squared:", r_value**2)
plt.plot(xs, IPSL.data, 'orange', label='IPSL')

ys = ECEARTH.data
slope, intercept, r_value, p_value, std_err = stats.linregress(xs,ys)
print (slope,r_value,p_value)
plt.plot(xs, intercept + slope*xs, 'r')
print ("r-squared:", r_value**2)
plt.plot(xs, ECEARTH.data, 'red', label='EC-EARTH')

years = ('2000', '2001', '2002', '2003', '2004', '2005','2006','2007','2008','2009','2010','2011','2012','2013','2014')
x_pos = (0,12,24,36,48,60,72,84,96,108,120,132,144,156,168)
plt.xticks(x_pos, years)
plt.legend()
plt.title('NBP')
plt.ylabel('PgC')
plt.show()






