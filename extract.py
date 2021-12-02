import os
import sys
from os.path import join as pjoin
import logging
import argparse
import getpass
import datetime
from calendar import monthrange

import numpy as np
import xarray as xr
import pandas as pd
import geopandas as gpd

from git import Repo

from warnings import filterwarnings
filterwarnings('ignore',category=FutureWarning)

LOGGER = logging.getLogger()

r = Repo('', search_parent_directories=True)
COMMIT = str(r.commit('HEAD'))

MODULE_DESCRIPTION="""
Extract a time series of point data from ERA5 reanalysis data for a single year"""

def to_dataframe(da, lons, lats, timeperiod='1D', func=np.max, method='nearest'):
    """
    Extract data from an `xarray.DataArray` at the given locations,
    resampling to the time period given by `timeperiod`
    """
    outdf = pd.DataFrame([])
    for idx, row in stnlist.iterrows():
        df = (
            da
            .sel(latitude=lats, longitude=lons, method=method)
            .resample(time=timeperiod)
            .reduce(func)
            .to_dataframe()
        )
        outdf = outdf.append(df, sort=True)
    return outdf

def to_dataset(filename, lons, lats, timeperiod='1D', func=np.max):
    """
    Extract data from `filename` at the given locations, resampling to the time period
    given by `timeperiod`. This will extract all variables in the given file.

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
        LOGGER.exception(f"Cannot open {filename}")
        return

    ts = (
        ds
        .interp(longitude=lons, latitude=lats)
        .resample(time=timeperiod)
        .reduce(func)
    )

    return ts


def main():
    startime = datetime.datetime.ctime(datetime.datetime.now())
    p = argparse.ArgumentParser(description=MODULE_DESCRIPTION)

    p.add_argument('-y', '--year', help="Year to process", default='2020')
    p.add_argument('-p', '--path', help="Base path to data")
    p.add_argument('-o', '--output', help="Output path")
    p.add_argument('-v', '--variable', help="Variable name", default="i10fg")
    p.add_argument('-s', '--stations', help="Station listing")
    args = p.parse_args()

    year = int(args.year)
    datapath = args.path
    variable = args.variable
    outputpath = args.output
    stationfile = args.stations

    logFile = f"extract.{variable}.{year}.log"
    logLevel = "INFO"
    logging.basicConfig(filename=logFile, level=logLevel,
                        format="%(asctime)s: %(funcName)s: %(message)s",
                        filemode='w', datefmt="%Y-%m-%d %H:%M:%S")
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(getattr(logging, logLevel))
    formatter = logging.Formatter('%(asctime)s: %(funcName)s:  %(message)s',
                                  datefmt='%H:%M:%S', )
    console.setFormatter(formatter)
    LOGGER.addHandler(console)

    if not os.path.isdir(outputpath):
        try:
            os.makedirs(outputpath)
        except OSError:
            LOGGER.exception(f"Cannot create output directory {outputpath}")
            raise

    stnlist = gpd.read_file(stationfile)
    stnid = stnlist.stnNum.values
    stnx = stnlist.stnLon.values
    stny = stnlist.stnLat.values
    LOGGER.info(f"There are {len(stnid)} stations")
    lons = xr.DataArray(stnx, coords=[stnid], dims=['stnid'])
    lats = xr.DataArray(stny, coords=[stnid], dims=['stnid'])

    fullts = []
    for month in range(1, 13):
        startdate = datetime.datetime(year, month, 1)
        enddate = datetime.datetime(year, month, monthrange(year, month)[1])
        filedatestr = f"{startdate.strftime('%Y%m%d')}-{enddate.strftime('%Y%m%d')}"
        tfile = os.path.join(datapath, f'{year}', f'{variable}_era5_oper_sfc_{filedatestr}.nc')
        # Get the daily maximum value at each point:
        LOGGER.info(f"Extracting time series data from {tfile}")
        ts = to_dataset(tfile, lons, lats, '1D', np.max)
        fullts.append(ts)

    ds = xr.concat(fullts, dim='time')
    ds.attrs['title'] = "ERA5 single-level reanalysis instantaneous_10m_wind_gust"
    ds.attrs['license'] = "Licence to use Copernicus Products: https://apps.ecmwf.int/datasets/licences/copernicus/"
    ds.attrs['source'] = "Fifth generation ECMWF atmospheric reanalysis of the global climate"
    ds.attrs['local_source'] = datapath
    ds.attrs['code_version'] = COMMIT
    ds.attrs['created_by'] = getpass.getuser()
    ds.attrs['history'] = f"{startime}: {' '.join(sys.argv)}"
    ds.to_netcdf(os.path.join(outputpath, f"ts.{variable}.{year}.nc"))

if __name__ == "__main__":
    main()
