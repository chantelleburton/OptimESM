#Same as 3 panels script but makes does Boreal, Tropics and Global

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
from iris.analysis.cartography import area_weights
import iris.analysis.cartography
from iris.coord_systems import GeogCS
import cf_units
from matplotlib.colors import BoundaryNorm, ListedColormap
from iris.coords import DimCoord
import cftime
from scipy import stats

# ── Global settings ───────────────────────────────────────────────────────────
date           = iris.Constraint(time=lambda cell: 2000 <= cell.point.year <= 2014)
ChoosePeriod   = 'Annual'
DataFolder     = '/data/scratch/chantelle.burton/OptimESM/NBP/'
var_constraint = iris.Constraint(name="nbp")

# ── Region definitions ────────────────────────────────────────────────────────
REGIONS = [
    {
        'name': 'Boreal',
        'lat' : lambda cell: 60.0 < cell < 90.0,
        'lon' : None,
    },
    {
        'name': 'Tropics',
        'lat' : lambda cell: -20.0 < cell < 20.0,
        'lon' : None,
    },
    {
        'name': 'Global',
        'lat' : lambda cell: -90.0 < cell < 90.0,
        'lon' : None,
    },
]

# ── Model definitions ─────────────────────────────────────────────────────────
MODELS = [
    {
        'label'   : 'CNRM',
        'color'   : 'red',
        'sign'    : 1,
        'r1_file' : DataFolder + 'nbp_Lmon_CNRM-ESM2-1_esm-hist_r1i1p2f2_gr_185001-201412.nc',
        'landfrac': DataFolder + 'sftlf_fx_CNRM-ESM2-1_esm-hist_r1i1p2f2_gr.nc',
        'lf_name' : 'sftlf',
        'ensemble': [],
    },
    {
        'label'   : 'UKESM',
        'color'   : 'black',
        'sign'    : 1,
        'r1_file' : DataFolder + 'nbp_Lmon_UKESM1-2-LL_esm-hist_r1i1p1f1_gn_195001-201412.nc',
        'landfrac': DataFolder + 'qrparm.landfrac.nc',
        'lf_name' : None,
        'ensemble': [
            DataFolder + 'nbp_Lmon_UKESM1-2-LL_esm-hist_r2i1p1f1_gn_195001-201412.nc',
            DataFolder + 'nbp_Lmon_UKESM1-2-LL_esm-hist_r3i1p1f1_gn_195001-201412.nc',
        ],
    },
    {
        'label'   : 'IPSL',
        'color'   : 'goldenrod',
        'sign'    : -1,   # wrong sign — pers comms Lars & Paul @ Lund
        'r1_file' : DataFolder + 'nbp_Lmon_IPSL-CM6-ESMCO2_esm-hist_r1i1p3f1_gr_195001-201412.nc',
        'landfrac': DataFolder + 'sftlf_fx_CNRM-ESM2-1_esm-hist_r1i1p2f2_gr.nc',
        'lf_name' : 'sftlf',
        'ensemble': [
            DataFolder + 'nbp_Lmon_IPSL-CM6-ESMCO2_esm-hist_r2i1p3f1_gr_195001-201412.nc',
            DataFolder + 'nbp_Lmon_IPSL-CM6-ESMCO2_esm-hist_r3i1p3f1_gr_195001-201412.nc',
            DataFolder + 'nbp_Lmon_IPSL-CM6-ESMCO2_esm-hist_r4i1p3f1_gr_195001-201412.nc',
        ],
    },
    {
        'label'   : 'EC-EARTH',
        'color'   : 'purple',
        'sign'    : 1,
        'r1_file' : DataFolder + 'nbp_LPmon_EC-Earth3-ESM-1_esm-hist_r1i1p1f1_gr_20*.nc',
        'landfrac': DataFolder + 'sftlf_fx_EC-Earth3-ESM-1_esm-hist_r5i1p1f1_gr.nc',
        'lf_name' : 'sftlf',
        'ensemble': [
            DataFolder + "nbp_LPmon_EC-Earth3-ESM-1_esm-hist_r2i1p1f1_gr_*.nc",
            DataFolder + "nbp_LPmon_EC-Earth3-ESM-1_esm-hist_r3i1p1f1_gr_*.nc",
            DataFolder + "nbp_LPmon_EC-Earth3-ESM-1_esm-hist_r5i1p1f1_gr_*.nc"
        ],
    },
]


# ── Observations ──────────────────────────────────────────────────────────────
def read_obs_nbp():
    pwdin = "/data/users/eleanor.burke/obs_data/nbp"
    cams_cube       = iris.load_cube(pwdin + "/cams/cams73_latest_co2_flux_surface_mm.nc")
    carboscope_cube = iris.load_cube(pwdin + "/carboscope/r76nbetEXToc_v2025.flux_land.mon.nc")
    gcp2024_cube    = iris.load_cube(pwdin + "/ct2022/GCP2024.flux1x1-monthly.processed.nc")
    return cams_cube, carboscope_cube, gcp2024_cube

OBS_CUBES  = read_obs_nbp()
OBS_LABELS = ["CAMS", "CarboScope", "GCP2024"]
OBS_STYLES = [dict(color=f"C{i}", ls='dashed') for i in range(len(OBS_LABELS))]

# ── Helper functions ──────────────────────────────────────────────────────────
def remove_duplicate_times(cube):
    """Remove duplicate time points based on year coordinate."""
    try:
        year_points = cube.coord('year').points
    except iris.exceptions.CoordinateNotFoundError:
        year_points = cube.coord('time').points
    _, unique_idx = np.unique(year_points, return_index=True)
    return cube[np.sort(unique_idx)]


def get_years(cube):
    """Extract year values using the year coordinate if present."""
    try:
        return cube.coord('year').points.astype(int)
    except iris.exceptions.CoordinateNotFoundError:
        time_coord = cube.coord('time')
        dates = time_coord.units.num2date(time_coord.points)
        return np.array([d.year for d in dates])

def CollapseToTimeseries(cube):
    """Area-weighted spatial sum → PgC/yr (Annual) or PgC/month."""
    coords = ('longitude', 'latitude')
    for coord in coords:
        if not cube.coord(coord).has_bounds():
            cube.coord(coord).guess_bounds()
    area = iris.analysis.cartography.area_weights(cube, normalize=False)
    cube = cube.collapsed(coords, iris.analysis.SUM, weights=area) / 1e12
    cube = cube * 86400 * 365 / 12
    return cube


def apply_region(cube, region):
    """Subset a cube to a region dict. Handles optional longitude cut."""
    if region['lon'] is not None:
        cube = cube.intersection(longitude=(-180, 180))
        cube = cube.extract(
            iris.Constraint(latitude=region['lat'], longitude=region['lon'])
        )
    else:
        cube = cube.extract(iris.Constraint(latitude=region['lat']))
    return cube


def load_and_preprocess(filepath, landfrac_file, lf_name, sign, is_glob=False):
    """Load one NBP file (or glob), apply land-fraction, deduplicate, return cube."""
    if is_glob:
        cubelist = iris.load(filepath, var_constraint)
        iris.util.equalise_attributes(cubelist)
        cube = cubelist.concatenate_cube()
    else:
        cube = iris.load_cube(filepath, var_constraint)

    cube = cube.extract(date)

    if landfrac_file is not None:
        if lf_name is not None:
            lf = iris.load_cube(landfrac_file, lf_name) / 100.0
        else:
            lf = iris.load_cube(landfrac_file)
        try:
            cube = cube * lf
        except Exception:
            pass

    iris.coord_categorisation.add_season_year(cube, 'time', name='year')
    if ChoosePeriod == 'Annual':
        cube = cube.aggregated_by(['year'], iris.analysis.SUM)

    cube = remove_duplicate_times(cube)
    return cube * sign


def timeseries_for_region(filepath, landfrac_file, lf_name, sign, region, is_glob=False):
    """Full pipeline: load → preprocess → region extract → collapse."""
    cube = load_and_preprocess(filepath, landfrac_file, lf_name, sign, is_glob)
    cube = apply_region(cube, region)
    return CollapseToTimeseries(cube)


# ── Main plotting loop ────────────────────────────────────────────────────────
fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
fig.subplots_adjust(hspace=0.12)

for ax, region in zip(axes, REGIONS):

    # ---- Observations --------------------------------------------------------
    for ob_cube, ob_label, ob_style in zip(OBS_CUBES, OBS_LABELS, OBS_STYLES):
        ob = ob_cube.copy()
        ob = ob.extract(date)
        if ChoosePeriod == 'Annual':
            ob = ob.aggregated_by(['year'], iris.analysis.SUM)
        ob = remove_duplicate_times(ob)
        ob = apply_region(ob, region)
        ob_ts = CollapseToTimeseries(ob)
        years = get_years(ob_ts)
        ax.plot(years, ob_ts.data, label=ob_label, **ob_style)

    # ---- Models --------------------------------------------------------------
    for mdef in MODELS:
        is_glob = '*' in mdef['r1_file']

        r1_ts    = timeseries_for_region(
            mdef['r1_file'], mdef['landfrac'], mdef['lf_name'],
            mdef['sign'], region, is_glob=is_glob
        )
        r1_years = get_years(r1_ts)

        if mdef['ensemble']:
            ens_years_list = []
            ens_data       = []
            for ens_file in mdef['ensemble']:
                ens_is_glob = '*' in ens_file
                try:
                    ens_ts = timeseries_for_region(
                        ens_file, mdef['landfrac'], mdef['lf_name'],
                        mdef['sign'], region, is_glob=ens_is_glob
                    )
                    ens_years_list.append(get_years(ens_ts))
                    ens_data.append(ens_ts.data)
                except Exception as e:
                    print(f"  WARNING: {ens_file}: {e}")

            if ens_data:
                all_year_sets = [set(r1_years.tolist())] + [set(y.tolist()) for y in ens_years_list]
                common_years  = np.array(sorted(set.intersection(*all_year_sets)))

                def subset(data_array, years_array):
                    idx = np.array([i for i, y in enumerate(years_array.tolist())
                                    if y in common_years.tolist()], dtype=int)
                    return np.array(data_array)[idx]

                r1_common  = subset(r1_ts.data, r1_years)
                ens_common = [subset(d, y) for d, y in zip(ens_data, ens_years_list)]

                # Stack all members including r1 and compute mean, min, max
                all_members  = np.vstack([r1_common] + ens_common)
                ens_mean     = np.mean(all_members, axis=0)
                ens_min      = np.min(all_members, axis=0)
                ens_max      = np.max(all_members, axis=0)

                ax.fill_between(common_years, ens_min, ens_max,
                                color=mdef['color'], alpha=0.15, linewidth=0)
                # Solid line is now the ensemble mean
                ax.plot(common_years, ens_mean,
                        color=mdef['color'], label=mdef['label'])
        else:
            ax.plot(r1_years, r1_ts.data, color=mdef['color'], label=mdef['label'])
            

    ax.set_title(region['name'], fontsize=11, loc='left', fontweight='bold')
    ax.set_ylabel('NBP (PgC yr$^{-1}$)')
    ax.axhline(0, color='grey', linewidth=0.6, linestyle=':')
    ax.legend(loc='best', fontsize=8, ncol=2)

axes[-1].set_xlabel('Year')
axes[-1].xaxis.set_major_locator(plt.MultipleLocator(2))
fig.suptitle('Net Biome Productivity by Region', fontsize=13, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig('NBP_regions_ensemble.png', dpi=150, bbox_inches='tight')
plt.show()
