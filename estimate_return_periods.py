from __future__ import division

import pandas as pd
import numpy as np
import logging
from functools import wraps
import time
import os

from os.path import join as pjoin
from datetime import datetime
from scipy.stats import genpareto


from plot import plotDiagnostics, plotFit
from return_period import returnLevels
from Utilities.files import flStartLog

CONVERTERS = {'Speed': lambda s: float(s or 0)/3.6}

def parse(yr, month, day, time):
    """
    Parse year, month and day as strings and return a datetime.
    
    Handles the case of a missing time string (Pandas returns nan
    if the field is empty).
    """
    if time is np.nan:
        time = '0000'
    timestr = '{0}-{1}-{2} {3}'.format(yr, month, day, time)

    return datetime.strptime(timestr, '%Y-%m-%d %H%M')

def timer(func):
    """
    Decorator to report execution time of a function/script.
    """
    @wraps(func)
    def wrap(*args, **kwargs):
        t1 = time.time()
        res = func(*args, **kwargs)

        tottime = time.time() - t1
        msg = "%02d:%02d:%02d " % \
          reduce(lambda ll, b: divmod(ll[0], b) + ll[1:],
                 [(tottime,), 60, 60])

        log.info("Time for {0}: {1}".format(func.func_name, msg))
        return res

    return wrap

def find_nearest_index(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx

def calculateShape(mu, data):
    """
    :param float mu: threshold parameter for the GPD distribution.
    :param data: :class:`numpy.ndarray` of data values to fit.
    """
    nobs = len(data)
    nexc = len(data[data > mu])
    rate = float(nexc)/float(nobs)
    gpd = genpareto.fit(data[data > mu] - mu)

    return gpd

@timer
def selectThreshold(data, minexc=10):
    """
    Select an appropriate threshold for fitting a generalised pareto
    distribution.
    The only constraint placed on the selection is that the shape
    parameter is negative (such that the distribution is bounded).
    :param data: :class:`numpy.ndarray` containing the observed values (with
                 missing values removed).
    :param int minexc: Minimum number of exceedances required.
    :returns: tuple of the shape, scale and threshold.
    """

    sh = []
    sc = []
    t = []
    q1000list = []
    q10000list = []

    eps = -0.01
    nobs = len(data)
    mu = np.median(data)
    while mu < data.max():
        #    for mu in np.arange(np.median(data), data.max(), 0.002):
        nexc = len(data[data > mu])
        rate = nexc / nobs
        if nexc < minexc:
            break

        pp = calculateShape(mu, data)
        q1000, q10000 = returnLevels(np.array([1000, 10000]),
                                     mu, pp[0], pp[2], rate)
        if np.isnan(q1000) or np.isnan(q10000):
            continue

        qdiff = np.abs(q10000 - q1000)
        if pp[0] < eps: # and qdiff < 0.2*q10000:# and qdiff > -eps:
            t.append(mu)
            sh.append(pp[0])
            sc.append(pp[2])
            q1000list.append(q1000)
            q10000list.append(q10000)
        mu += 0.002

    if len(t) == 0:
        log.warn("No suitable shape parameters identified")
        return 0, 0, 0
    Av1000 = np.mean(np.array(q1000list))
    Av10000 = np.mean(np.array(q10000list))
    Av1000 = np.ceil(Av1000 + 0.05*Av1000)
    Av10000 = np.ceil(Av10000 + 0.05*Av10000)

    idx1000 = find_nearest_index(np.array(q1000list), Av1000)
    idx10000 = find_nearest_index(np.array(q10000list), Av10000)

    u1000 = t[idx1000]
    u10000 = t[idx10000]

    if u1000 > u10000:
        shmax = sh[idx1000]
        scmax = sc[idx1000]
    else:
        shmax = sh[idx10000]
        scmax = sc[idx10000]

    return shmax, scmax, u1000

def ymParser(dateStr):
    if dateStr == "":
        return np.nan
    else:
        return datetime.strptime(dateStr, '%m/%Y')

def readDataFile(filename):
    """
    Read a data file containing daily maximum wind speed records.

    """    
    # Column names of the data files:
    names = ['dc', 'StnNum', 'Year', 'Month', 'Day', 'Speed',
             'QSpeed', 'Dir', 'QDir', 'Time', 'QTime']

    df = pd.read_csv(filename, skipinitialspace=True,
                     skiprows=1, names=names,
                     parse_dates=[['Year', 'Month', 'Day', 'Time']], 
                     date_parser=parse, index_col=False)
                     #converters=CONVERTERS)

    return df

def readStationFile(stnFile):
    """
    Read a station details file to get a list of stations to process. The 
    ``dtOpen`` and ``dtClose`` fields are not converted to datetime instances,
    even though they represent dates. 
    
    :param str stnFile: Full path of the station file.
    
    :returns: :class:`pandas.DataFrame` holding details of the stations
    """
    log.info("Reading stations from {0}".format(stnFile))
    names = ["id", "stnNum", "rfdist", "stnName", "dtOpen", "dtClose", 
             "stnLat", "stnLon", "stnSource", "stnState", "stnElev", 
             "stnBaroElev", "stnWMO", "DataStartYear", "DataEndYear", 
             "PercentComplete", "PercentQuality", "PercentQualityN",
             "PercentQualityW", "PercentQualityS", "PercentQualityI"]
    
    dtype = {"stnNum":str}
    stndf = pd.read_csv(stnFile, skipinitialspace=True, names=names, 
                        index_col=False, dtype=dtype) #, skiprows=1) 
    log.info("Loaded details of {0} stations".format(len(stndf.index)))
    return stndf

#------------------------------------------------
# Main execution code:

input_path = "C:/WorkSpace/data/daily/input/"
output_path = "C:/WorkSpace/data/daily/output/2006"
logfile = pjoin(output_path, "estimate_return_period.log")
loglevel = "INFO"

log = flStartLog(logfile, loglevel, True, True)


stationFile = pjoin(input_path, "DC02D_StnDet_99999999720437.txt")
stndf = readStationFile(stationFile)
basefile = "DC02D_Data_{0}_99999999720437.txt"

gpdfile = open(pjoin(output_path, "gpdtable_ms.csv"), 'w')
gpdfile.write("Station, Name, shape, scale, threshold, rate\n")

rp = np.array([1, 2, 5, 10, 20, 50, 100, 200,
               500, 1000, 2000, 5000, 10000])

rpfile = open(pjoin(output_path, "rptable_ms.csv"), 'w')
rpfile.write("Station, Name," + ", ".join(rp.astype(str)) + "\n")

for i in stndf.index:
    stnNum = stndf['stnNum'][i]
    filename = pjoin(input_path, basefile.format(stnNum))
    dataRange = "({0} - {1})".format(stndf['DataStartYear'][i], 
                                     stndf['DataEndYear'][i])
    stnName = stndf['stnName'][i].title().strip() + " " + dataRange
    fitname = pjoin(output_path, '{0}_gpdfit.png'.format(stnNum))
    diagname = pjoin(output_path, '{0}_gpddiag.png'.format(stnNum))
    if os.path.exists(filename):
        log.info("Processing {0}".format(stnName))
        df = readDataFile(filename)
        quality = df['QSpeed'].fillna("X").map(lambda x: x in 
                                               ['Y', 'N', 'X', ' ', np.nan])
        dmax = df['Speed'][df['Speed'].notnull() & quality]
        if len(dmax) == 0:
            log.info("No valid data")
            continue
        xi, sigma, mu = selectThreshold(dmax, minexc=10)
        log.debug("Parameters: {0}, {1}, {2}".format(xi, sigma, mu))
        rate = float(len(dmax[dmax > mu])) / float(len(dmax))
        if xi == 0:
            continue
        plotFit(dmax, mu, xi, sigma, stnName, fitname)
        plotDiagnostics(dmax, mu, xi, sigma, diagname)

        gpdfile.write("{0}, {1}, {2:.6f}, {3:.6f}, {4:.3f}, {5:.4f}\n".
                      format(stnNum, stnName, xi, sigma, mu, rate))
        rpvals = returnLevels(rp, mu, xi, sigma, rate)
        rpstr = ", ".join(['{:.3f}']*len(rpvals)).format(*rpvals)
        rpfile.write("{0}, {1}, {2}\n".format(stnNum, stnName, rpstr))
    else:
        log.info("No data file for {0}".format(stnName))
gpdfile.close()
rpfile.close()
