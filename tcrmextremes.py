import os
from os.path import join as pjoin, realpath, isdir, dirname
import sys
import logging as log
import traceback
import atexit
import warnings

import argparse
import numpy as np
import seaborn as sns
from scipy.stats import genpareto

import matplotlib
matplotlib.use('Agg', warn=False)  # Use matplotlib backend
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, LogLocator, NullFormatter
plt.ioff()

from extremes import gpdSelectThreshold, returnLevels, empReturnPeriod, calculateUncertainty

import database
from Utilities.files import flStartLog, flConfigFile
from Utilities.config import ConfigParser
from Utilities.parallel import attemptParallel
from Utilities.version import version

#sns.set_context("poster")
sns.set_style("whitegrid")

majloc = LogLocator()
minloc = LogLocator(base=10.0, subs=np.arange(0.1, 1, 0.1), numticks=12)

mpl_logger = log.getLogger('matplotlib')
mpl_logger.setLevel(log.WARNING)

def runFit(recs, locId, locName, numYears, outputPath):
    """
    Run a GPD fitting routine, using both the threshold selection algorithm
    and a simple regression fit (using `scipy.stats.rv_distribution.fit`).

    """
    log.info("Processing {0}".format(locName))

    # Run the threshold selection algorithm:
    xi, sigma, mu = gpdSelectThreshold(recs['wspd'], nexc=10)

    # Determine a GPD fit using a fixed threshold (99.5 percentile):
    thresh = np.percentile(recs['wspd'], 99.9)
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

    params, crp, lrp, urp = calculateUncertainty(recs['wspd'], rp, *gpd)
    if crp is None:
        return None
    sortedmax = np.sort(data)
    fig, ax1 = plt.subplots(1, 1, figsize=(10, 8))
    ax1.semilogx(rp, rval, label=r"Iterative threshold $\mu$ = {0:.4f}".format(mu))
    ax1.semilogx(rp, crp, label=r"Percentile threshold $\mu$ = {0:.4f}".format(thresh),
                 color='0.5')
    ax1.semilogx(rp, lrp, label=r"90% CI", color='0.5', linestyle='--')
    ax1.semilogx(rp, urp, color='0.5', linestyle='--')
    ax1.scatter(emprp[emprp > 1], sortedmax[emprp > 1], s=100,
                color='r', label="Empirical ARI")

    #title_str = (locName   + "\n" +
    #             r"$\mu$ = {0:.3f}, $\xi$ = {1:.5f}, $\sigma$ = {2:.4f}".
    #             format(mu, xi, sigma))
    title_str = locName
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
    ax1.text(15000, 45.6, 'Cat 3', ha='center')
    ax1.text(15000, 62.5, 'Cat 4', ha='center')
    ax1.text(15000, 77.8, 'Cat 5', ha='center')
    ax1.legend(loc=2)
    fig.tight_layout()
    plt.savefig(pjoin(outputPath, "{0}.png".format(locId)),
                bbox_inches="tight")
    plt.close()
    return locId, mu, xi, sigma, rate, rval, thresh, rate2, (params), rval2, lrp, urp

def processCI(recs, locId, locName, numYears, outputPath, pctl=99.):
    print(f"Plotting ARI curve for {locName} ({locId})")
    wspd = recs['wspd'][recs['wspd'] > 0]
    years = np.array([1, 2, 5, 10, 20, 50, 100, 200,
                      500, 1000, 2000, 5000, 10000])
    data = np.zeros(int(365.25 * numYears))
    data[-len(wspd):] = wspd
    emprp = empReturnPeriod(data)

    xi, sigma, mu = gpdSelectThreshold(wspd, nexc=10)
    rate = float(len(data[data > mu])) / float(len(data))
    if xi == 0:
        itrval = np.zeros(len(years))
    else:
        itrval = returnLevels(years, mu, xi, sigma, rate)

    threshold = np.percentile(wspd, pctl)
    rate2 = float(len(data[data > threshold])) / float(len(data))
    gpd = genpareto.fit(wspd[wspd > threshold], loc=threshold)


    params, crp, lrp, urp = calculateUncertainty(wspd[wspd > threshold], years, *gpd)
    if crp is None:
        return None
    fig, ax = plt.subplots(1, 1)
    ax.semilogx(years, crp, label=rf"Fitted ARI ($\mu = P_{{{pctl}}}$)")
    ax.semilogx(years, lrp, '0.5', label=r"90% CI", ls='--')
    ax.semilogx(years, urp, '0.5', ls='--')
    ax.axvline(500, color='k', zorder=1)
    ax.text(500, 0, "1:500", ha='right', va='bottom', rotation='vertical', fontsize='x-small')
    ax.axvline(2000, color='k', zorder=1)
    ax.text(2000, 0, "1:2000", ha='right', va='bottom', rotation='vertical', fontsize='x-small')
    ax.semilogx(years, itrval, label=rf'Iterative ARI ($\mu = {{{mu:.4f}}}$)', color='k')
    ax.scatter(emprp[emprp > 1], data[emprp > 1], s=50, alpha=0.5,
               color='r', label="Empirical ARI")

    ax.set_ylim((0, 120))
    ax.set_yticks(np.arange(0, 121, 10))
    ax.set_xlim((1, 10000))
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.set_ylabel('Wind speed [m/s]')
    ax.set_xlabel('Average recurrence interval [years]')
    ax.grid(which='major', linestyle='-')
    ax.grid(which='minor', linestyle='--', linewidth=1)
    ax.set_title("{0}".format(locName))
    ax.legend(loc=2)
    fig.tight_layout()
    plt.savefig(pjoin(outputPath, "{0}.png".format(locId)),
                bbox_inches="tight")
    #plt.close()
    aep = 1. / years
    fig2, ax2 = plt.subplots(1, 1)
    ax2.semilogy(crp, aep, label=rf'Fitted AEP ($\mu = P_{{{pctl}}}$)')
    ax2.semilogy(lrp, aep, '0.5', ls='--', label=r"90% CI")
    ax2.semilogy(urp, aep, '0.5', ls='--')
    ax2.semilogy(itrval, aep, color='k', label=rf'Iterative AEP  ($\mu = {{{mu:.4f}}}$)')
    ax2.scatter(data[emprp > 1], 1./emprp[emprp > 1], s=50, alpha=0.5,
                color='r', label='Empirical AEP')
    ax2.yaxis.set_major_locator(majloc)
    ax2.yaxis.set_minor_locator(minloc)
    ax2.yaxis.set_minor_formatter(NullFormatter())
    ax2.legend(loc=1, fontsize='small')

    ax2.set_ylabel("Annual exceedance probability [%]")
    ax2.set_xlabel("Wind speed [m/s]")
    ax2.set_xlim((0, 120))
    ax2.axhline(1./500, color='k', zorder=1)
    ax2.text(0, 1./500, "1:500", ha='left', va='bottom', fontsize='x-small')
    ax2.axhline(1./2000, color='k', zorder=1)
    ax2.text(0, 1./2000, "1:2000", ha='left', va='bottom', fontsize='x-small')
    ax2.grid(which='major', linestyle='-')
    ax2.grid(which='minor', linestyle='--', linewidth=0.5)
    fig2.tight_layout()

    plt.savefig(pjoin(outputPath, "{0}.AEP.png".format(locId)),
                bbox_inches="tight")
    plt.close()

    return locId, mu, xi, sigma, rate, crp, threshold, rate2, (params), crp, lrp, urp

def main(configFile):
    """
    Execute the fitting routine for all locations in the database.

    """
    log.info("Running tcrmextremes")
    config = ConfigParser()
    config.read(configFile)
    numYears = config.getint('TrackGenerator', 'NumSimulations')
    outputPath = config.get("Output", "Path")
    plotPath = pjoin(outputPath, "plots/fit")
    processPath = pjoin(outputPath, "process")
    db = database.HazardDatabase(configFile)
    locations = db.getLocations()
    log.info("There are {0} locations in the database".\
             format(len(locations)))
    locNameList = list(locations['locName'])
    locIdList = list(locations['locId'])
    comm.barrier()

    status = MPI.Status()

    work_tag = 0
    result_tag = 1

    rp = np.array([1, 2, 5, 10, 20, 50, 100, 200,
                   500, 1000, 2000, 5000, 10000])
    paramheader = ("locId, locName, it_shape, it_thresh, it_scale, it_rate,"
                   " gpd_rate, gpd_shape, gpd_thresh, gpd_scale\n")
    paramfmt = "{}, {}, " + ", ".join(["{:.8f}"] * 8) + "\n"
    rvalheader = "locId, locName, " + ", ".join(["{:d}"]*len(rp)).format(*rp) + "\n"
    rvalfmt = "{}, {}, " + ", ".join(["{:.2f}"] * len(rp)) + "\n"
    # On the head node:
    if (comm.rank == 0) and (comm.size > 1):
        fh = open(pjoin(processPath, "parameters.csv"), "w")
        rval1fh = open(pjoin(processPath, "iterative_rl.csv"), "w")
        rval2fh = open(pjoin(processPath, "fitted_rl.csv"), "w")
        rvalufh = open(pjoin(processPath, "fitted_rl_u.csv"), "w")
        rvallfh = open(pjoin(processPath, "fitted_rl_l.csv"), "w")
        fh.write(paramheader)
        rval1fh.write(rvalheader)
        rval2fh.write(rvalheader)
        rvalufh.write(rvalheader)
        rvallfh.write(rvalheader)

        w = 0
        p = comm.size - 1

        for d in range(1, comm.size):
            if w < len(locations):
                locName = locNameList[w]
                locId = locations['locId'][locNameList.index(locName)]
                log.info("Running calculations for {0}".format(locName))
                args = (locId, locName, numYears, plotPath)
                comm.send(args, dest=d, tag=work_tag)
                w += 1
            else:
                comm.send(None, dest=d, tag=work_tag)
                p = w

        terminated = 0

        while terminated < p:
            result = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
            log.debug("Receiving results from node {0}".\
                      format(status.source))
            if result is None:
                pass
            else:
                locId, mu, sigma, xi, rate1, rval1, thresh, rate2, gpd, rval2, lrp, urp = result
                locName = locations['locName'][locIdList.index(locId)]
                fh.write(paramfmt.format(locId, locName, sigma, mu, xi, rate1, rate2, *gpd))
                rval1fh.write(rvalfmt.format(locId, locName, *rval1))
                rval2fh.write(rvalfmt.format(locId, locName, *rval2))
                rvalufh.write(rvalfmt.format(locId, locName, *urp))
                rvallfh.write(rvalfmt.format(locId, locName, *lrp))

            d = status.source
            if w < len(locNameList):
                locName = locNameList[w]
                locId = locations['locId'][locNameList.index(locName)]
                log.info("Running calculations for {0}".format(locName))
                args = (locId, locName, numYears, plotPath)
                comm.send(args, dest=d, tag=work_tag)
                w += 1
            else:
                comm.send(None, dest=d, tag=work_tag)
                terminated += 1

        fh.close()

    # On the worker nodes:
    elif (comm.size > 1) and (comm.rank != 0):
        while True:
            args = comm.recv(source=0, tag=work_tag, status=status)
            if args is None:
                break

            locId, locName, numYears, plotPath = args
            log.info("Processing {0} on node {1}".\
                     format(args[1], comm.rank))
            recs = database.locationRecords(db, locId)
            result = processCI(recs, locId, locName, numYears, plotPath)
            comm.send(result, dest=0, tag=result_tag)

    elif comm.size == 1 and comm.rank == 0:
        # Assume no Pypar
        db = database.HazardDatabase(configFile)
        locations = db.getLocations()
        log.info("There are {0} locations in the database".\
                 format(len(locations)))
        locNameList = list(locations['locName'])
        fh = open(pjoin(processPath, "parameters.csv"), "w")
        rval1fh = open(pjoin(processPath, "iterative_rl.csv"), "w")
        rval2fh = open(pjoin(processPath, "fitted_rl.csv"), "w")
        rvalufh = open(pjoin(processPath, "fitted_rl_u.csv"), "w")
        rvallfh = open(pjoin(processPath, "fitted_rl_l.csv"), "w")

        fh.write(paramheader)
        rval1fh.write(rvalheader)
        rval2fh.write(rvalheader)
        rvalufh.write(rvalheader)
        rvallfh.write(rvalheader)

        for locName in locNameList:
            log.info("Running calculations for {0}".format(locName))
            locId = locations['locId'][locNameList.index(locName)]
            recs = database.locationRecords(db, locId)
            args = (recs, locId, locName, numYears, plotPath)
            result = processCI(recs, locId, locName, numYears, plotPath)
            if result is None:
                continue
            else:
                locId, mu, sigma, xi, rate1, rval1, thresh, rate2, gpd, rval2, lrp, urp = result
                fh.write(paramfmt.format(locId, locName, sigma, mu, xi, rate1, rate2, *gpd))
                rval1fh.write(rvalfmt.format(locId, locName, *rval1))
                rval2fh.write(rvalfmt.format(locId, locName, *rval2))
                rvalufh.write(rvalfmt.format(locId, locName, *urp))
                rvallfh.write(rvalfmt.format(locId, locName, *lrp))

        fh.close()
        rval1fh.close()
        rval2fh.close()
        rvalufh.close()
        rvallfh.close()


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

    global MPI, comm
    MPI = attemptParallel()
    atexit.register(MPI.Finalize)
    comm = MPI.COMM_WORLD
    if comm.size > 1 and comm.rank > 0:
        # MPI execution:
        logfile += '-' + str(comm.rank)
        verbose = False
    else:
        pass

    log = flStartLog(logfile, logLevel, verbose, datestamp)
    log.info("Code version: {0}".format(version()))
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
            # Gracefully handle failure on parallel execution:
            if comm.size > 1:
                for d in range(1, comm.size):
                    comm.send(None, dest=d, tag=0)

    comm.barrier()
    log.info("Finished running {0}".format(sys.argv[0]))

if __name__ == "__main__":
    startup()
