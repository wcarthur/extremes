from os.path import join as pjoin
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator


import pandas as pd

from extremes import empReturnPeriod

import seaborn as sns
sns.set_context("paper")
sns.set_style("whitegrid")

majloc = LogLocator()
minloc = LogLocator(base=10.0, subs=np.arange(0.1, 1, 0.1), numticks=12)

obsPath = "X:/georisk/HaRIA_B_Wind/projects/tcha/data/derived/observations/daily"
obsPath = "X:/georisk/HaRIA_B_Wind/data/derived/tcobs"
stationFilePath = "X:/georisk/HaRIA_B_Wind/data/raw/from_bom/2019/Daily"
stnfile = pjoin(stationFilePath, "DC02D_StnDet_999999999632559_updated.txt")

STNTYPES = [('st', 'S2'), ('stnId', 'i'), ('stnDistCode', 'S4'), ('stnName', 'S'), 
            ('stnDateOpen', 'S10'), ('stnDateClosed', 'S10'), ('stnLat', 'f8'), 
            ('stnLon', 'f8'), ('method', 'S15'), ('state', 'S3'), 
            ('stnElevation', 'f8'), ('baroElev', 'i'), ('stnWMONumber', 'i'), ('stnDataStart', 'i'), 
            ('stnDataEnd', 'i'), ('blank', 'S3'), ('percentcomplete', 'f8'), ('pcqualy', 'f8'), 
            ('pcqualn', 'f8'), ('pcqualw', 'f8'), ('pcquals', 'f8'), ('pcquali', 'f8'), ('end', 'S1')]
STNCONVERT = {'stnName' : str.rstrip}


def loadTCObs(obsfile: str) -> pd.DataFrame:
    """
    Load a file containing observed TC wind gusts. Ensure the date/time columns are `pd.Datetime` types.

    :param str obsfile: Path to the file containing observations

    :returns: `pd.DataFrame` containing records of observed TC wind gusts
    """
    df = pd.read_csv(obsfile)
    df['dtObs'] = pd.to_datetime(df['dtObs'])
    df['dtTC'] = pd.to_datetime(df['dtTC'])
    return df

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

def calcObservedEP(stnNum: int, df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate the exceedance probability for wind speeds based on observed data

    :param int stnNum: BoM station number
    :param df: `pd.DataFrame` of observed wind speeds at that station

    :returns: updated `pd.DataFrame` with the exceedance probability and recurrence interval for each wind speed

    """

    df.sort_values('gust', axis=0, inplace=True)
    numYears = getStationDates(stndf, stnNum)
    data = np.zeros(int(numYears * 365.25))
    wspd = np.sort(np.array(df['gust']))
    data[-len(wspd):] = wspd
    df['pprp'] = empReturnPeriod(data)[-len(wspd):]
    df['ppaep'] = 1. - np.exp(-1./df['pprp'])
    return df

def plotObsEP(stnNum: int, obsdf: pd.DataFrame):
    """
    Plot an EP curve for the set of observarions

    :param int stnNum: BoM station number
    :param obsdf: `pd.DataFrame` containing TC-related observations

    :returns: None
    """
    fig, ax = plt.subplots(1, 1,)
    ax.scatter(obsdf.gust, obsdf.ppaep, s=50, color='k', marker='x',
               label='Observed AEP')
    locName = stndf.loc[stnNum, 'stnName']
    locState = stndf.loc[stnNum, 'state']
    locStart = stndf.loc[stnNum, 'stnDataStart']
    locEnd = stndf.loc[stnNum, 'stnDataEnd']
    ax.set_title(f"{locName} ({locState}) {locStart}-{locEnd}")
    ax.set_yscale('log')
    ax.grid(which='major', linestyle='-')
    ax.grid(which='minor', linestyle='--', linewidth=1)
    ax.set_xlim((0, 100))
    ax.set_xticks(np.arange(0, 101, 10))
    ax.set_ylim((10e-3, 1))
    ax.set_xlabel('Wind speed [m/s]')
    ax.set_ylabel('Annual Exceedance Probability')
    #ax.yaxis.set_minor_formatter(NullFormatter())
    fig.tight_layout()
    plt.savefig(pjoin(obsPath, f"ppari_{stnNum:06d}.png"), bbox_inches='tight')
    plt.close()
 

def _calcObservedEP(locId):
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

"""locId = 200283
obsdf = _calcObservedEP(locId)
obsdf.to_csv(pjoin(obsPath, f"ppari_{locId:06d}.csv"), index=False)
"""
stndf = pd.read_csv(stnfile, parse_dates=[4, 5],
                    usecols=(1,2,3,4,5,6,7,9,10,12,13,14,16), 
                    names = np.dtype(STNTYPES).names,
                    skiprows=1, engine='python', index_col='stnId', 
                    converters=STNCONVERT)

tcobsdf = loadTCObs(pjoin(obsPath, "stncpa_obs.csv"))

for stnNum, obs in tcobsdf.groupby('stnNum'):
    if len(obs) > 10:
        print(f"Calculating exceedance probability for {stnNum}")
        obsdf = calcObservedEP(stnNum, obs)
        obsdf.to_csv(pjoin(obsPath, f"ppari_{stnNum:06d}.csv"), index=False)
        plotObsEP(stnNum, obsdf)
    else:
        print(f"Fewer than 10 observations for {stnNum}")
