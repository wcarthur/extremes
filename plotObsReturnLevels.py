# coding: utf-8

# Comparison of observed and simulated ARI wind speeds
# 
# This notebook plots the average recurrence interval (ARI) wind
# speeds based on observed wind speeds corresponding to the passage of
# TCs (within 200 km of a station). It adds a plot of the fitted ARI
# wind speeds from a TCRM simulation.

import os
import io
import sys
from os.path import join as pjoin
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


import pandas as pd
import geopandas as gpd
from datetime import datetime

from extremes import returnLevels, empReturnPeriod
from distributions import fittedPDF

import seaborn as sns
sns.set_context("poster")
sns.set_style("whitegrid")


def loadObservations(stnId):
    names = ['recid', 'stnId', 'datetime', 'gust',
             'direction', 'quality', 'cycName']
    
    filename = pjoin(obsPath, "bom_{0:06d}.csv".format(stnId))
    try:
        obsdf = pd.read_csv(filename, 
                            parse_dates=[2], infer_datetime_format=True)
    except FileNotFoundError:
        print("No data file for stnId: {0}".format(stnId))
        return None
    except TypeError:
        import pdb; pdb.set_trace()
    return obsdf

def getStationDates(stnId):
    startYear = stndf.loc[stnId]['stnDataStart']
    endYear = stndf.loc[stnId]['stnDataEnd']
    numYears = endYear - startYear + 1
    return numYears


STNTYPES = [('st', 'S2'), ('stnId', 'i'), ('stnDistCode', 'S4'),
            ('stnName', 'S'), ('stnDateOpen', 'S10'),
            ('stnDateClosed', 'S10'), ('stnLat', 'f8'),
            ('stnLon', 'f8'), ('method', 'S15'), ('state', 'S3'),
            ('stnElevation', 'f8'), ('baroElev', 'i'),
            ('stnWMONumber', 'i'), ('stnDataStart', 'i'),
            ('stnDataEnd', 'i'), ('blank', 'S3'),
            ('percentcomplete', 'f8'), ('pcqualy', 'f8'),
            ('pcqualn', 'f8'), ('pcqualw', 'f8'), ('pcquals', 'f8'),
            ('pcquali', 'f8'), ('end', 'S1')]
STNCONVERT = {'stnName' : str.rstrip}


# Start with loading the observation station information. This is from
# the daily maximum wind gust dataset (Geosciene Australia eCat
# #110561), starting with the station details file.

obsPath = "C:/WorkSpace/data/derived/tcobs/daily"
outputPath = "C:/Workspace/temp/tcobs"
stationFilePath = "C:/WorkSpace/data/raw/daily_max_wind_gust/"
stnfile = pjoin(stationFilePath, "DC02D_StnDet_999999999425050.txt")

stndf = pd.read_csv(stnfile, parse_dates=[4, 5],
                    usecols=(1,2,3,4,5,6,7,9,10,12,13,14,16),
                    names = np.dtype(STNTYPES).names,
                    skiprows=1, engine='python', index_col='stnId',
                    converters=STNCONVERT)
stationNameList = list(stndf['stnName'])

# Now load a shape file that contains the observed stations joined
# with the TCRM simulation locations. Note in this dataframe, we need
# to add an index, and so we index by both the location id number
# (TCRM simulation locations) *and* the station number (observations).

locationFilePath = "C:/WorkSpace/data/derived/tcobs/merged.shp"
locdf = gpd.read_file(locationFilePath)
locdf = locdf.set_index(["locId", 'stnId'])
locationNameList = list(locdf['Station_Na'])

# Indexing using the `locId` first, then selecting the `index`
# attribute returns the `stnId`, which is used to load the observed
# data


# The parameters of the fitted distribution are contained in another
# data file, and this is indexed using the TCRM location id number.

paramFile = "C:/WorkSpace/data/derived/tc/tcha/parameters.csv"
paramNames = ['locId', "locName", "it_scale", "it_shape", "it_thresh", 
              "it_rate", "gpd_rate", "gpd_shape", "gpd_thresh", "gpd_scale"]
gpddf = pd.read_csv(paramFile, names=paramNames, skiprows=1, index_col='locId')

def plotObservedHazard(locId, ax):
    obsdf = loadObservations(locId)
    if obsdf is None:
        return ax
    numYears = getStationDates(locId)
    data = np.zeros(int(numYears * 365.25))
    # Include conversion to 0.2 second wind gust
    wspd = np.sort(np.array(obsdf['gust']))*1.114 
    if len(wspd) == 0:
        return ax
    data[-len(wspd):] = wspd
    emprp = empReturnPeriod(data)
    
    ax.scatter(emprp[emprp > 1], data[emprp > 1], s=50,
                color='k', marker='x', label="Empirical ARI")
    return ax

def plotFittedHazard(gpd_params, ax):
    """
    Plot a fitted distribution, with approximate 90% confidence interval
    and empirical return period values.

    :param data: :class:`numpy.ndarray` of observed data values.
    :param float mu: Selected threshold value.
    :param float xi: Fitted shape parameter.
    :param float sigma: Fitted scale parameter.
    :param str title: Title string for the plot.
    :param str figfile: Path to store the file (includes image format)

    """
    
    rp = np.array([1, 2, 5, 10, 20, 50, 100, 200,
                   500, 1000, 2000, 5000, 10000])
    mu, xi, sigma, rate = gpd_params
    rval = returnLevels(rp, mu, xi, sigma, rate)

    ax.semilogx(rp, rval, label="Fitted hazard curve")
    return ax


def loadParameters(locationName):
    locId = locdf.index[locationNameList.index(locationName)][0]   
    try:
        stnId = locdf.loc[locId].index[0]
    except KeyError:
        print("No index for given location id: {0}".format(locId))
    else:
        stnName = locdf.loc[locId, stnId]['Station_Na']

        stnObsFile = pjoin(obsPath, "bom_{0:06d}.csv".format(stnId))
        if os.path.exists(stnObsFile):
            print("Observation file exists for {0}".format(stnName))
        else:
            print("No observations for {0}".format(stnName))
            
    if locId in gpddf.index.values:
        gpd_rate = gpddf.loc[locId]['gpd_rate']
        gpd_shape = gpddf.loc[locId]['gpd_shape']
        gpd_scale = gpddf.loc[locId]['gpd_scale']
        gpd_thresh = gpddf.loc[locId]['gpd_thresh']
        
        fig, ax = plt.subplots(1, 1)
        plotFittedHazard((gpd_thresh, gpd_shape, gpd_scale, gpd_rate), ax)
        plotObservedHazard(stnId, ax)
        
        title_str = (stnName)  # + "\n" +
                 #r"$\mu$ = {0:.3f}, $\xi$ = {1:.5f}, $\sigma$ = {2:.4f}".
                 #format(mu, xi, sigma))
        ax.set_title(title_str)
        ax.set_ylim((0, 100))
        ax.set_yticks(np.arange(0, 101, 10))
        ax.set_xlim((1, 10000))
        ax.set_ylabel('Wind speed (m/s)')
        ax.set_xlabel('Average recurrence interval (years)')
        ax.grid(which='major', linestyle='-')
        ax.grid(which='minor', linestyle='--', linewidth=1)
        ax.axhline(45.6, c='lime', linestyle='--', linewidth=2)
        ax.axhline(62.5, c='darkorange', linestyle='--', linewidth=2)
        ax.axhline(77.8, c='darkred', linestyle='--', linewidth=2)
        ax.text(20000, 45.6, 'Cat 3', ha='center')
        ax.text(20000, 62.5, 'Cat 4', ha='center')
        ax.text(20000, 77.8, 'Cat 5', ha='center')
        ax.legend(loc=2)
        fig.tight_layout()
        plt.savefig(pjoin(outputPath, "{0:06d}.png".format(stnId)),
                    bbox_inches='tight')
        plt.clf()
    else:
        print("No index in GPD parameter file for {0}".format(locId))
        
    


for locName in locationNameList:
    loadParameters(locName)

