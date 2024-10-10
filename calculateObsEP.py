import os
import numpy as np
import pandas as pd

from extremes import ppoints



def loadData(inputPath, stnnum):
    """
    Load data from a csv file. Remove all records where the quality flag is not 'Y', and any missing values.

    :param str inputPath: Directory where observation data files are stored
    :param int stnnum: BoM Station Number

    :returns: `np.array` of gust wind speed values.

    """
    fname = f"DC02D_Data_{stnnum:06d}_999999999632559.txt"
    df = pd.read_csv(os.path.join(inputPath, fname), skipinitialspace=True)
    gust = df[df['Quality of maximum gust speed']=='Y']['Speed of maximum wind gust in m/s'].dropna().values

    return gust

def calcEP(gust, datafile, npyr=365.25):
    n = len(gust)
    empprob = ppoints(n)
    empep = 1 - np.exp(npyr * np.log(empprob))
    idx = np.where(empep < 0.999)[0]
    np.savetxt(datafile,
               np.vstack([np.sort(gust)[idx], empep[idx]]).T,
               delimiter=",", fmt="%.5f")


# -----------------------------------------------------------------
inputPath = r"X:\georisk\HaRIA_B_Wind\data\raw\from_bom\2019\Daily"
outputPath = r"X:\georisk\HaRIA_B_Wind\data\derived\obs\daily_max_wind\wind\obsep"
stationfile = os.path.join(inputPath,'DC02D_StnDet_999999999632559_updated.txt')
stndf = pd.read_csv(stationfile)
if not os.path.isdir(outputPath):
    os.makedirs(outputPath)

for idx, stn in stndf.iterrows():
    stnnum = stn['Bureau of Meteorology Station Number']
    stnname = stn['Station Name'].strip()
    print(f"Processing {stnname} ({stnnum})")
    gust = loadData(inputPath, stnnum)
    if len(gust) < 365.25 * 5:
        print(f"< 5 years data for {stnnum}")
        continue
    epdatafile = os.path.join(outputPath, f"{stnnum:06d}.ep.csv")
    calcEP(gust, epdatafile)