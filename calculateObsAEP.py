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
stationFilePath = "X:/georisk/HaRIA_B_Wind/data/derived/tcobs"
stnfile = pjoin(stationFilePath, "DC02D_StnDet_999999999632559_M3.csv")

STNTYPES = [('st', 'S2'), ('stnId', 'i'), ('stnDistCode', 'S4'), ('stnName', 'S'),
            ('stnDateOpen', 'S10'), ('stnDateClosed', 'S10'), ('stnLat', 'f8'),
            ('stnLon', 'f8'), ('method', 'S15'), ('state', 'S3'),
            ('stnElevation', 'f8'), ('baroElev', 'i'), ('stnWMONumber', 'i'), ('stnDataStart', 'i'),
            ('stnDataEnd', 'i'), ('percentcomplete', 'f8'), ('pcqualy', 'f8'),
            ('pcqualn', 'f8'), ('pcqualw', 'f8'), ('pcquals', 'f8'), ('pcquali', 'f8'), ('end', 'S1'),
            ('rindex', 'i',), ('stName', 'S'), ('tname', 'S'),
            ('Mt_n', 'f8'), ('Mz_n', 'f8'), ('M3_n', 'f8'),
            ('Mt_ne', 'f8'), ('Mz_ne', 'f8'), ('M3_ne', 'f8'),
            ('Mt_e', 'f8'), ('Mz_e', 'f8'), ('M3_e', 'f8'),
            ('Mt_se', 'f8'), ('Mz_se', 'f8'), ('M3_se', 'f8'),
            ('Mt_s', 'f8'), ('Mz_s', 'f8'), ('M3_s', 'f8'),
            ('Mt_sw', 'f8'), ('Mz_sw', 'f8'), ('M3_sw', 'f8'),
            ('Mt_w', 'f8'), ('Mz_w', 'f8'), ('M3_w', 'f8'),
            ('Mt_nw', 'f8'), ('Mz_nw', 'f8'), ('M3_nw', 'f8'),]
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

def getStationM3(stndf, stnId):
    """
    Retrieve the Mt and M3 values for each direction for a given station

    :param stndf:

    """
    directions = ['n', 'ne', 'e', 'se', 's', 'sw', 'w', 'nw']
    Mtcols = [f"Mt_{d}" for d in directions]
    Mzcols = [f"Mz_{d}" for d in directions]
    try:
        mt = stndf.loc[stnId][Mtcols]
        # Sometimes the direction is missing, so set value to 1
        mt["Mt_nan"] = 1.0
        mz = stndf.loc[stnId][Mzcols]
        mz["Mz_nan"] = 1.0

    except KeyError:
        mt = np.ones((9,))
        mz = np.ones((9,))
        return mt, mz
    return mt, mz

def categoriseDirection(df: pd.DataFrame) -> pd.DataFrame:
    """
    Categorise the direction to enable selection of appropriate directional site multiplier

    :param df: `pd.Dataframe` containing the data

    :returns: updated `pd.DataFrame` with additional column "dir"
    """
    labels = ['n', 'ne', 'e', 'se', 's', 'sw', 'w', 'nw', 'n']
    bins = [0, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5, 292.5, 337.5, 360]

    df['dir'] = pd.cut(df.direction, bins=bins, right=True, ordered=False, labels=labels)
    return df

def calcObservedEP(stnNum: int, df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate the exceedance probability for wind speeds based on observed data

    :param int stnNum: BoM station number
    :param df: `pd.DataFrame` of observed wind speeds at that station

    :returns: updated `pd.DataFrame` with the exceedance probability and recurrence interval for each wind speed

    """
    df = categoriseDirection(df)
    mt, mz = getStationM3(stndf, stnNum)
    mtvals = mt.loc[[f"Mt_{x}" for x in df['dir']]]
    mzvals = mz.loc[[f"Mz_{x}" for x in df['dir']]]

    # Calculate standardised wind speed by removing local topographic and terrain effects
    df['gust_std'] = df['gust'] / mtvals.values / mzvals.values
    df.sort_values('gust_std', axis=0, inplace=True)
    numYears = getStationDates(stndf, stnNum)
    data = np.zeros(int(numYears * 365.25))
    wspd = np.sort(np.array(df['gust_std']))
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
    ax.scatter(obsdf.gust_std, obsdf.ppaep, s=50, color='k', marker='x',
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

stndf = pd.read_csv(stnfile, parse_dates=[4, 5],
                    usecols=(1,2,3,4,5,6,7,9,10,12,13,14,16,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49),
                    names = np.dtype(STNTYPES).names,
                    skiprows=1, engine='python', index_col='stnId',
                    converters=STNCONVERT)

tcobsdf = loadTCObs(pjoin(obsPath, "stncpa_obs.csv"))

for stnNum, obs in tcobsdf.groupby('stnNum'):
    if stnNum not in stndf.index:
        # Assume no multiplier data available
        # Extraction only completed for QLD as at Sept 2022
        continue
    if len(obs) > 10: # Need a decent number of obs
        print(f"Calculating exceedance probability for {stnNum}")
        obsdf = calcObservedEP(stnNum, obs)
        obsdf.to_csv(pjoin(obsPath, f"ppari_{stnNum:06d}.csv"), index=False)
        plotObsEP(stnNum, obsdf)
    else:
        print(f"Fewer than 10 observations for {stnNum}")
