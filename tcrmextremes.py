from __future__ import print_function, division

from os.path import join as pjoin, realpath, isdir, dirname, exists
import logging as log
import os
import sys
import traceback
from datetime import datetime


import argparse
import database
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import statsmodels.api as sm

from scipy.stats import genpareto, scoreatpercentile

from Utilities.config import ConfigParser
from return_period import returnLevels, empiricalReturnPeriod, returnPeriodUncertainty
from calculate_return_periods import selectThreshold
from distributions import fittedPDF


from Utilities.files import flStartLog, flConfigFile

from Utilities.config import ConfigParser

import seaborn as sns
sns.set_context("poster")
sns.set_style("whitegrid")

def run_fit(recs, locId, locName, outputPath):
    log.info("Processing {0}".format(locName))
    xi, sigma, mu = selectThreshold(recs['wspd'])
    thresh = np.percentile(recs['wspd'], 95)
    gpd = genpareto.fit(recs['wspd'][recs['wspd'] > thresh], floc=thresh)
    data = np.zeros(int(5000 * 365.25))
    data[-len(recs):] = recs['wspd']
    rp = np.array([1, 2, 5, 10, 20, 50, 100, 200,
                   500, 1000, 2000, 5000, 10000])

    rate = float(len(data[data > mu])) / float(len(data))
    if xi == 0:
        rval = np.zeros(len(rp))
    else:
        rval = returnLevels(rp, mu, xi, sigma, rate)

    rate2 = float(len(data[data > thresh])) / float(len(data))
    rval2 = returnLevels(rp, thresh, gpd[0], gpd[2], rate2)

    emprp = empiricalReturnPeriod(data)
    sortedmax = np.sort(data)
    fig, ax1 = plt.subplots(1, 1)
    ax1.semilogx(rp, rval, label="Fitted hazard curve")
    ax1.semilogx(rp, rval2 , label=r"$\mu$ = {0}".format(max(data)/2), color='0.5')
    ax1.scatter(emprp[emprp > 1], sortedmax[emprp > 1], s=100,
                    color='r', label="Empirical ARI")

    title_str = (locName   + "\n" +
                 r"$\mu$ = {0:.3f}, $\xi$ = {1:.5f}, $\sigma$ = {2:.4f}".
                 format(mu, xi, sigma))
    ax1.set_title(title_str)
    ax1.set_ylim((0, 100))
    ax1.set_yticks(np.arange(0, 101, 10))
    ax1.set_xlim((1, 10000))
    ax1.set_ylabel('Wind speed (m/s)')
    ax1.set_xlabel('Average recurrence interval (years)')
    ax1.grid(which='major', linestyle='-')
    ax1.grid(which='minor', linestyle='--', linewidth=1)
    ax1.axhline(45.6, c='lime', linestyle='--', linewidth=2)#, label='Cat 3')
    ax1.axhline(62.5, c='darkorange', linestyle='--', linewidth=2)#, label='Cat 4')
    ax1.axhline(77.8, c='darkred', linestyle='--', linewidth=2)#, label='Cat 5')
    ax1.text(20000, 45.6, 'Cat 3', ha='center')
    ax1.text(20000, 62.5, 'Cat 4', ha='center')
    ax1.text(20000, 77.8, 'Cat 5', ha='center')
    ax1.legend(loc=2)
    plt.savefig(pjoin(outputPath, "{0}.png".format(locId)))
    plt.close()
    return (mu, xi, sigma), (gpd)

def main(configFile):
    log.info("Running tcrmextremes")
    config = ConfigParser()
    config.read(configFile)
    numYears = config.getint('TrackGenerator', 'NumSimulations')

    db = database.HazardDatabase(configFile)
    locations = db.getLocations()
    log.info("There are {0} locations in the database".format(len(locations)))
    locNameList = list(locations['locName'])
    outputPath = ""#"/c/WorkSpace/data/Derived/tc/fitting"

    fh = open(pjoin(outputPath, "parameters.csv"), "w")
    for loc in locNameList:
        locId = locations['locId'][locNameList.index(loc)]
        recs = database.locationRecords(db, locId)
        (mu, xi, sigma), (gpd) = run_fit(recs, locId, loc, outputPath)
        fh.write("{0}, {1:.5f}, {2:.5f}, {3:.5f}, {4:.5f}, {5}\n".format(loc, xi, sigma, mu, recs['wspd'].max()/2, gpd))

    fh.close()
    
    

def startup():
    """
    Parse command line arguments and call the :func:`main` function.
    """
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_file',
                        help='Path to configuration file')
    parser.add_argument('-v', '--verbose', help='Verbose output',
                        action='store_true')
    parser.add_argument('-d', '--debug', help='Allow pdb traces',
                        action='store_true')
    args = parser.parse_args()
    print(args)
    configFile = args.config_file
    config = ConfigParser()
    config.read(configFile)
    logfile = config.get('Logging', 'LogFile')
    logdir = dirname(realpath(logfile))

    # If log file directory does not exist, create it
    if not isdir(logdir):
        try:
            os.makedirs(logdir)
        except OSError:
            logfile = flConfigFile('.log')

    logLevel = config.get('Logging', 'LogLevel')
    verbose = config.getboolean('Logging', 'Verbose')
    datestamp = config.getboolean('Logging', 'Datestamp')
    debug = False

    if args.verbose:
        verbose = True

    if args.debug:
        debug = True

    log = flStartLog(logfile, logLevel, verbose, datestamp)
    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    warnings.filterwarnings("ignore", category=UserWarning, module="pytz")
    warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
    warnings.filterwarnings("ignore", category=UserWarning,
                            module="matplotlib")

    warnings.filterwarnings("ignore", category=RuntimeWarning)

    if debug:
        main(configFile)
    else:
        try:
            main(configFile)
        except ImportError as e:
            log.critical("Missing module: {0}".format(e.strerror))
        except Exception:  # pylint: disable=W0703
            # Catch any exceptions that occur and log them (nicely):
            tblines = traceback.format_exc().splitlines()
            for line in tblines:
                log.critical(line.lstrip())
                
if __name__ == "__main__":
    startup()
