import os

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd

from datetime import datetime
from calendar import monthrange

import seaborn as sns

from warnings import filterwarnings
filterwarnings('ignore',category=FutureWarning)

varname = 'i10fg'
year = 2020
basepath = "/g/data/rt52/era5/single-levels/reanalysis/"
varname = "i10fg"

dataPath = os.path.join(basepath, varname)

outputPath = "/scratch/w85/cxa547/reanalysis/timeseries"
# Load stations from a file:
# Just needs to be three lists - a station id, the longitude and latitude coordinates
stndf = gpd.read_file("stationlist.shp")
stnid = stndf.stnNum.values
stnx = stndf.stnLon.values
stny = stndf.stnLat.values

x = xr.DataArray(stnx, coords=[stnid], dims=['stnid'])
y = xr.DataArray(stny, coords=[stnid], dims=['stnid'])

def to_dataset(filename, lons, lats, timeperiod='1D', func=np.max):
    """
    Extract data from `filename` at the given locations, resampling to the time period
    given by `timeperiod`

    :param str filename: path to the netcdf file
    :param lons: `xarray.DataArray` of x-coordinates
    :param lats: `xarray.DataArray` of y-coordinates
    :param str timeperiod: time period string to resample the time series data. Default '1D'
    :param func: `callable()` function which can be called in the form _func(x, axis=axis, **kwargs)_
                  to return the result of collapsing an np.ndarray over an integer valued axis

    """
    try:
        ds = xr.open_dataset(filename)
    except:
        log.exception(f"Cannot open {filename}")

    ts = ds.interp(
            longitude=lons,
            latitude=lats
            ).resample(time=timeperiod).reduce(func)

    return ts


fullts = []
for month in range(1, 4):
    startdate = datetime(year, month, 1)
    enddate = datetime(year, month, monthrange(year, month)[1])
    filedatestr = f"{startdate.strftime('%Y%m%d')}-{enddate.strftime('%Y%m%d')}"
    tfile = os.path.join(dataPath, f'{year}', f'{varname}_era5_oper_sfc_{filedatestr}.nc')
    # Get the daily maximum value at each point:
    print(f"Extracting time series data from {tfile}")
    ts = to_dataset(tfile, x, y, '1D', np.max)
    fullts.append(ts)


ds = xr.concat(fullts, dim='time')
ds.to_netcdf(os.path.join(outputPath, f"ts.{varname}.{year}.nc"))

