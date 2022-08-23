import os
import io
import sys
from os.path import join as pjoin
import numpy as np
import matplotlib.pyplot as plt

import pandas as pd
import geopandas as gpd
from datetime import datetime

from extremes import returnLevels, empReturnPeriod
from distributions import fittedPDF

import seaborn as sns
sns.set_context("talk")
sns.set_style("whitegrid")

obsPath = "X:/georisk/HaRIA_B_Wind/projects/tcha/data/derived/observations/daily"
obsPath = "X:/georisk/HaRIA_B_Wind/data/derived/tcobs"
stationFilePath = "X:/georisk/HaRIA_B_Wind/data/raw/from_bom/daily/2017/"
stnfile = pjoin(stationFilePath, "DC02D_StnDet_999999999425050.txt")

STNTYPES = [('st', 'S2'), ('stnId', 'i'), ('stnDistCode', 'S4'), ('stnName', 'S'), 
            ('stnDateOpen', 'S10'), ('stnDateClosed', 'S10'), ('stnLat', 'f8'), 
            ('stnLon', 'f8'), ('method', 'S15'), ('state', 'S3'), 
            ('stnElevation', 'f8'), ('baroElev', 'i'), ('stnWMONumber', 'i'), ('stnDataStart', 'i'), 
            ('stnDataEnd', 'i'), ('blank', 'S3'), ('percentcomplete', 'f8'), ('pcqualy', 'f8'), 
            ('pcqualn', 'f8'), ('pcqualw', 'f8'), ('pcquals', 'f8'), ('pcquali', 'f8'), ('end', 'S1')]
STNCONVERT = {'stnName' : str.rstrip}
stndf = pd.read_csv(stnfile, parse_dates=[4, 5],
                        usecols=(1,2,3,4,5,6,7,9,10,12,13,14,16), 
                        names = np.dtype(STNTYPES).names,
                        skiprows=1, engine='python', index_col='stnId', 
                        converters=STNCONVERT)

def loadObservations(stnId):
    """
    Load the observations from file for a given BoM station, where the observations have
    been selected from the complete digital history of daily maximum wind speeds, where a cyclone has 
    passed within 200 km of the station, and the station was open at the time of passage.
    
    :param int stnId: Bureau of Meteorology Station identification number
    
    :returns: data frame containing the gust wind speed, direction and cyclone name
              based on passage of cyclones near the selected station. If no observation
              file is found (`Exception.FileNotFoundError`), return `None`
    """
    
    names = ['recid', 'stnId', 'datetime', 'gust',
             'direction', 'quality', 'cycName', 'cycId']
    
    filename = pjoin(obsPath, "bom_{0:06d}.csv".format(stnId))
    try:
        obsdf = pd.read_csv(filename, skiprows=1, names=names,
                            parse_dates=[2], infer_datetime_format=True)
    except IOError:
        print("No data file for stnId: {0}".format(stnId))
        return None
    return obsdf

def getStationDates(stndf, stnId):
    """
    Retrieve the length of observed record in years, based on the start and end dates
    of the observational data. 
    
    :param int stnId: Bureau of Meteorology Station identification number
    
    :returns: number of years between the start and end of the data on file.
    """
    startYear = stndf.loc[stnId]['stnDataStart']
    endYear = stndf.loc[stnId]['stnDataEnd']
    numYears = endYear - startYear + 1
    return numYears

def calcObservedEP(locId):
    obsdf = loadObservations(locId)
    obsdf.sort_values('gust', axis=0, inplace=True)
    if obsdf is None:
        return None #ax
    numYears = getStationDates(stndf, locId)
    data = np.zeros(int(numYears * 365.25))
    wspd = np.sort(np.array(obsdf['gust'])) # Include conversion to 0.2 second wind gust
    data[-len(wspd):] = wspd
    obsdf['pprp'] = empReturnPeriod(data)[-len(wspd):]
    obsdf['ppaep'] = 1. - np.exp(-1./obsdf['pprp'])
    return obsdf
    #ax.scatter(emprp[emprp > 1], data[emprp > 1], s=50,
    #            color='k', marker='x', label="Observed ARI", zorder=100)
    #return ax

locId = 200283
obsdf = calcObservedEP(locId)
obsdf.to_csv(pjoin(obsPath, f"ppari_{locId:06d}.csv"), index=False)