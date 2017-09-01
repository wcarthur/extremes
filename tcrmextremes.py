from __future__ import print_function, division

import os
from os.path import join as pjoin, realpath, isdir, dirname
import logging as log
import traceback

import argparse
import database
import numpy as np

import matplotlib
matplotlib.use('Agg', warn=False)  # Use matplotlib backend
import matplotlib.pyplot as plt
plt.ioff()

from scipy.stats import genpareto

from extremes import gpdSelectThreshold, returnLevels, empReturnPeriod

from Utilities.files import flStartLog, flConfigFile
from Utilities.config import ConfigParser
from Utilities.parallel import attemptParallel
from Utilities.version import version


import seaborn as sns
sns.set_context("poster")
sns.set_style("whitegrid")


def runFit(recs, locId, locName, numYears, outputPath):
    """
    Run a GPD fitting routine, using both the threshold selection algorithm
    and a simple regression fit (using `scipy.stats.rv_distribution.fit`).

    """
    log.info("Processing {0}".format(locName))

    # Run the threshold selection algorithm:
    xi, sigma, mu = gpdSelectThreshold(recs['wspd'])

    # Determine a GPD fit using a fixed threshold (99.5 percentile):
    thresh = np.percentile(recs['wspd'], 99.5)
    gpd = genpareto.fit(recs['wspd'][recs['wspd'] > thresh], floc=thresh)

    data = np.zeros(int(numYears * 365.25))
    data[-len(recs):] = recs['wspd']
    rp = np.array([1, 2, 5, 10, 20, 50, 100, 200,
                   500, 1000, 2000, 5000, 10000])

    # Rate of exceedances above the algorithmically-selected threshold:
    rate = float(len(data[data > mu])) / float(len(data))
    if xi == 0:
        rval = np.zeros(len(rp))
    else:
        rval = returnLevels(rp, mu, xi, sigma, rate)

    # Rate of exceedances above 99.5th percentile:
    rate2 = float(len(data[data > thresh])) / float(len(data))
    rval2 = returnLevels(rp, thresh, gpd[0], gpd[2], rate2)

    emprp = empReturnPeriod(data)
    sortedmax = np.sort(data)
    fig, ax1 = plt.subplots(1, 1)
    ax1.semilogx(rp, rval, label="Fitted hazard curve")
    ax1.semilogx(rp, rval2, label=r"$\mu$ = {0}".format(max(data)/2),
                 color='0.5')
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
    ax1.axhline(45.6, c='lime', linestyle='--', linewidth=2)
    ax1.axhline(62.5, c='darkorange', linestyle='--', linewidth=2)
    ax1.axhline(77.8, c='darkred', linestyle='--', linewidth=2)
    ax1.text(20000, 45.6, 'Cat 3', ha='center')
    ax1.text(20000, 62.5, 'Cat 4', ha='center')
    ax1.text(20000, 77.8, 'Cat 5', ha='center')
    ax1.legend(loc=2)
    fig.tight_layout()
    plt.savefig(pjoin(outputPath, "{0}.png".format(locId)))
    plt.close()
    return locId, mu, xi, sigma, thresh, (gpd)

def main(configFile):
    """
    Execute the fitting routine for all locations in the database.

    """
    log.info("Running tcrmextremes")
    config = ConfigParser()
    config.read(configFile)
    numYears = config.getint('TrackGenerator', 'NumSimulations')
    outputPath = config.get("Output", "Path")
    plotPath = pjoin(outputPath, "plots/")
    processPath = pjoin(outputPath, "process")
    db = database.HazardDatabase(configFile)
    locations = db.getLocations()
    log.info("There are {0} locations in the database".\
             format(len(locations)))
    locNameList = list(locations['locName'])
    pp.barrier()

    work_tag = 0
    result_tag = 1

    # On the head node:
    if (pp.rank() == 0) and (pp.size() > 1):
        #db = database.HazardDatabase(configFile)
        #locations = db.getLocations()
        #log.info("There are {0} locations in the database".\
        #         format(len(locations)))
        #locNameList = list(locations['locName'])
        fh = open(pjoin(processPath, "parameters.csv"), "w")
        w = 0
        p = pp.size() - 1

        for d in range(1, pp.size()):
            if w < len(locations):
                locName = locNameList[w]
                locId = locations['locId'][locNameList.index(locName)]
                log.info("Running calculations for {0}".format(locName))
                args = (locId, locName, numYears, plotPath)
                pp.send(args, destination=d, tag=work_tag)
                w += 1
            else:
                pp.send(None, destination=d, tag=work_tag)
                p = w

        terminated = 0

        while terminated < p:
            result, status = pp.receive(pp.any_source, tag=result_tag,
                                        return_status=True)
            log.debug("Receiving results from node {0}".\
                      format(status.source))
            locId, mu, sigma, xi, thresh, gpd = result
            locName = locations['locName'][locNameList.index(locId)]
            fh.write("{0}, {1:.5f}, {2:.5f}, {3:.5f}, {4:.5f}, {5}\n".
                     format(locName, xi, sigma, mu, thresh, gpd))

            d = status.source
            if w < len(locNameList):
                locName = locNameList[w]
                locId = locations['locId'][locNameList.index(locName)]
                log.info("Running calculations for {0}".format(locName))
                args = (locId, locName, numYears, plotPath)
                pp.send(args, destination=d, tag=work_tag)
                w += 1
            else:
                pp.send(None, destination=d, tag=work_tag)
                terminated += 1

    elif (pp.size() > 1) and (pp.rank() != 0):
        while True:
            args = pp.receive(source=0, tag=work_tag)
            if args is None:
                break

            locId, locName, numYears, plotPath = args
            log.info("Processing {0} on node {1}".\
                     format(args[1], pp.rank()))
            recs = database.locationRecords(db, locId)
            result = runFit(recs, locId, locName, numYears, plotPath)
            pp.send(result, destination=0, tag=result_tag)

    elif pp.size() == 1 and pp.rank() == 0:
        # Assume no Pypar
        db = database.HazardDatabase(configFile)
        locations = db.getLocations()
        log.info("There are {0} locations in the database".\
                 format(len(locations)))
        locNameList = list(locations['locName'])
        fh = open(pjoin(processPath, "parameters.csv"), "w")

        for locName in locNameList:
            log.info("Running calculations for {0}".format(locName))
            locId = locations['locId'][locNameList.index(locName)]
            recs = database.locationRecords(db, locId)
            args = (recs, locId, locName, numYears, plotPath)

            locId, mu, sigma, xi, thresh, gpd =\
                        runFit(recs, locId, locName, numYears, plotPath)
            fh.write("{0}, {1:.5f}, {2:.5f}, {3:.5f}, {4:.5f}, {5}\n".
                     format(locName, xi, sigma, mu, thresh, gpd))

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

    global pp
    pp = attemptParallel()
    import atexit
    atexit.register(pp.finalize)

    if pp.size() > 1 and pp.rank() > 0:
        # MPI execution:
        logfile += '-' + str(pp.rank())
        verbose = False
    else:
        pass

    log = flStartLog(logfile, logLevel, verbose, datestamp)
    log.info("Code version: {0}".format(version()))
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
            tblines = traceback.format_exc().splitlines()
            for line in tblines:
                log.critical(line.lstrip())
        except Exception:  # pylint: disable=W0703
            # Catch any exceptions that occur and log them (nicely):
            tblines = traceback.format_exc().splitlines()
            for line in tblines:
                log.critical(line.lstrip())

if __name__ == "__main__":
    startup()
