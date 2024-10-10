import os
import argparse
import logging
from datetime import time, timedelta, datetime as dt
from configparser import ConfigParser, ExtendedInterpolation
from os.path import join as pjoin, realpath, isdir, dirname, splitext

import numpy as np
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri
from rpy2.rinterface_lib.sexp import NULLType
from statsmodels.distributions.empirical_distribution import ECDF
import pandas as pd
import matplotlib
matplotlib.use("Agg") # Preven bitmap allocation error
import matplotlib.pyplot as plt
import seaborn as sns

from process import pAlreadyProcessed, pWriteProcessedFile, pArchiveFile, pInit
from files import flStartLog, flGetStat
sns.set_style('whitegrid')

g_files = {}
g_output_path = os.getcwd()
LOGGER = logging.getLogger()

base = importr('base')
utils = importr('utils')
evmix = importr('evmix')
stats = importr('stats')
numpy2ri.activate()

np.random.seed(1234)

class lognormgpd():
    def __init__(self, fit):
        self.lnmean = fit.rx2('lnmean')[0]
        self.lnsd = fit.rx2('lnsd')[0]
        self.u = fit.rx2('u')[0]
        self.sigmau = fit.rx2('sigmau')[0]
        self.xi = fit.rx2('xi')[0]
        self.phiu = fit.rx2('phiu')[0]



def ppoints(n, a=0.375):
    """
    Generates the sequence of probabilities at which to evaluate an inverse distribution.

    :param n: either int or `np.array` - either number of points, or an array of observations.
    :param float a: offset fraction, typically in (0, 1), default=0.375
    """
    if isinstance(n, int):
        m = np.arange(1, n + 1)
        if n > 10:
            a = 0.5
    else:
        m = np.arange(1, len(n) + 1)
        if len(n) > 10:
            a = 0.5
    pp = (m - a)/(m[-1] + (1 - a) - a)
    return pp

def histplot(gust, paramdf, stnname, plotfile):
    """
    Plot a histogram of the data with the median fitted gamma-GPD (median u and median xi)

    :param gust: `np.array` of observed gust wind speed values
    :param paramdf: `pd.DataFrame` of the fitted parameter values
    :param str stnname: Name of the station for plot title
    :param str plotfile: Path to save the image to

    """

    xx = np.arange(0, gust.max()+1, 0.01)
    ximed = paramdf['xi'].median()
    idxxi = np.argmin(np.abs(paramdf['xi'] - ximed))
    pmedxi = paramdf.iloc[idxxi]
    umed = paramdf['u'].median()
    idxu = np.argmin(np.abs(paramdf['u'] - umed))
    pmedu = paramdf.iloc[idxu]
    dfitxi = evmix.dlognormgpd(xx, pmedxi.lnmean, pmedxi.lnsd, pmedxi.u, pmedxi.sigmau, pmedxi.xi, pmedxi.phiu)
    dfitu = evmix.dlognormgpd(xx, pmedu.lnmean, pmedu.lnsd, pmedu.u, pmedu.sigmau, pmedu.xi, pmedu.phiu)
    fig, ax = plt.subplots(1, 1)
    ax.hist(gust, bins=100, density=True, fc='white', ec='k')
    ax.plot(xx, dfitxi, color='red',
        label=rf"Bulk tail fraction: $(u, \xi ) = {{{pmedxi.u:.1f}}}, {{{pmedxi.xi:.3f}}}$")
    ax.plot(xx, dfitu, color='blue',
        label=rf"Bulk tail fraction: $(u, \xi ) = {{{pmedu.u:.1f}}}, {{{pmedu.xi:.3f}}}$")

    ax.legend(fontsize='xx-small',)
    ax.set_xlabel("Gust wind speed [m/s]")
    ax.set_ylabel("Density")
    ax.set_title(stnname)
    plt.savefig(plotfile, bbox_inches="tight")

def epplot(gust, paramdf, stnname, plotfile):
    n = len(gust)
    empprob = ppoints(n)
    transempprob = -1/np.log(empprob)
    minemppower = -np.log10(1 - 1/n/100)
    maxemppower = np.ceil(np.log10(np.max(transempprob))) + 1
    theprob = 1 - np.power(10, -np.arange(minemppower, maxemppower + 0.05, 0.05))
    theprob = np.sort(np.concatenate((theprob, 1 - theprob)))
    transtheprob = -1/np.log(theprob)
    ep = 1 - np.exp(-1/(transtheprob/365.25))
    empep = 1 - np.exp(-1/(transempprob/365.25))

    fig, ax = plt.subplots(1, 1)
    for idx, row in paramdf.iterrows():
        thequant = evmix.qlognormgpd(theprob, row.lnmean, row.lnsd, row.u, row.sigmau, row.xi, row.phiu)
        ax.plot(thequant, ep, color='0.75', linewidth=1)

    ximed = paramdf['xi'].median()
    umed = paramdf['u'].median()

    idxxi = np.argmin(np.abs(paramdf['xi'] - ximed)) # Median shape (xi < 0 only)
    idxu = np.argmin(np.abs(paramdf['u'] - umed)) # Median threshold
    try:
        idxxise = paramdf[paramdf['xise']!=0]['xise'].idxmin() # Minimum standard error of xi
        plotse = True
    except ValueError:
        print(f"No valid standard errors for {stnname}")
        plotse = False
    pmed = paramdf.iloc[idxxi]
    medquant = evmix.qlognormgpd(theprob, pmed.lnmean, pmed.lnsd, pmed.u, pmed.sigmau, pmed.xi, pmed.phiu)
    ax.plot(medquant, ep, color='r', linewidth=2,
        label=rf"$(u, \xi ) = {{{pmed.u:0.3f}}}, {{{pmed.xi:0.3f}}}$")

    pmed = paramdf.iloc[idxu]
    medquant = evmix.qlognormgpd(theprob, pmed.lnmean, pmed.lnsd, pmed.u, pmed.sigmau, pmed.xi, pmed.phiu)
    ax.plot(medquant, ep, color='b',
        linewidth=2, label=rf"$(u, \xi ) = {{{pmed.u:0.3f}}}, {{{pmed.xi:0.3f}}}$")

    # Minimum xi SE above the initial threshold
    if plotse:
        pmin = paramdf.iloc[idxxise]
        minquant = evmix.qlognormgpd(theprob, pmin.lnmean, pmin.lnsd, pmin.u, pmin.sigmau, pmin.xi, pmin.phiu)
        ax.plot(minquant, ep, color='g',
            linewidth=2, label=rf"$(u, \xi ) = {{{pmin.u:0.3f}}}, {{{pmin.xi:0.3f}}}$")

    idx = np.where(empep < 0.999)[0]
    xmin = int(np.sort(gust)[idx][0] / 5) * 5
    xmax = (int(gust.max() / 5) + 2) * 5
    ax.scatter(np.sort(gust)[idx], empep[idx], s=20, color='k', marker='x', label='Data', zorder=100,)
    ylims = (1e-3, 1)
    xlims = (xmin, xmax)
    ax.grid(which='major', linestyle='-')
    ax.grid(which='minor', linestyle='--', linewidth=1)
    ax.set_yscale('log')
    ax.set_ylabel("AEP")
    ax.set_xlabel("Return level [m/s]")
    ax.legend(loc=1, fontsize='xx-small')
    ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    ax.set_title(stnname)
    plt.savefig(plotfile, bbox_inches="tight")
    plt.close(fig)

def rlplot(gust, paramdf, stnname,  plotfile):
    """
    Plot the fitted gamma-GPD curve on a return level plot

    :param gust: `np.array` of observed gust wind speed values
    :param paramdf: `pd.DataFrame` of the fitted parameter values
    :param str stnname: Name of the station for plot title
    :param str plotfile: Path to save the image to
    """
    n = len(gust)
    empprob = ppoints(n)
    transempprob = -1/np.log(empprob)
    minemppower = -np.log10(1 - 1/n/100)
    maxemppower = np.ceil(np.log10(np.max(transempprob))) + 1
    theprob = 1 - np.power(10, -np.arange(minemppower, maxemppower + 0.05, 0.05))
    theprob = np.sort(np.concatenate((theprob, 1 - theprob)))
    transtheprob = -1/np.log(theprob)

    fig, ax = plt.subplots(1, 1)
    for idx, row in paramdf.iterrows():
        thequant = evmix.qlognormgpd(theprob, row.lnmean, row.lnsd, row.u, row.sigmau, row.xi, row.phiu)
        ax.plot(transtheprob/365.25, thequant, color='0.75', linewidth=1)

    ximed = paramdf['xi'].median()
    umed = paramdf['u'].median()
    idxxi = np.argmin(np.abs(paramdf['xi'] - ximed))
    idxu = np.argmin(np.abs(paramdf['u'] - umed))

    pmed = paramdf.iloc[idxxi]
    medquant = evmix.qlognormgpd(theprob, pmed.lnmean, pmed.lnsd, pmed.u, pmed.sigmau, pmed.xi, pmed.phiu)
    ax.plot(transtheprob/365.25, medquant, color='r', linewidth=2,
        label=rf"$(u, \xi ) = {{{pmed.u:0.3f}}}, {{{pmed.xi:0.3f}}}$")

    pmed = paramdf.iloc[idxu]
    medquant = evmix.qlognormgpd(theprob, pmed.lnmean, pmed.lnsd, pmed.u, pmed.sigmau, pmed.xi, pmed.phiu)
    ax.plot(transtheprob/365.25, medquant, color='b',
        linewidth=2, label=rf"$(u, \xi ) = {{{pmed.u:0.3f}}}, {{{pmed.xi:0.3f}}}$")

    ax.scatter(transempprob/365.25, np.sort(gust), s=20, color='k', marker='x', label='Data', zorder=100,)

    xlims = (0.05, 1000)
    ax.grid(which='major', linestyle='-')
    ax.grid(which='minor', linestyle='--', linewidth=1)
    ax.set_xscale('log')
    ax.set_xlabel("ARI [years]")
    ax.set_ylabel("Return level [m/s]")
    ax.legend(loc=4, fontsize='xx-small')
    ax.set_xlim(xlims)
    ax.set_title(stnname)
    plt.savefig(plotfile, bbox_inches="tight")
    plt.close(fig)

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

def calcMAE(gust, t, fit):
    """
    Calculate the MAE of the distribution fit

    :param gust: `np.array` of wind speed values
    :param float t: threshold value
    :param fit: `r` object containing the fit parameters

    :returns: MAE for the fit, separating above and below the threshold
    """
    sgust = np.sort(gust)
    empcdf = ECDF(sgust)
    fgpd = lognormgpd(fit)
    dfit = evmix.dlognormgpd(sgust, fgpd.u, fgpd.sigmau, fgpd.xi, fgpd.phiu)
    delta = np.abs(empcdf(sgust) - dfit)
    idx = np.argmin(np.abs(sgust - t))
    maeupper = np.mean(delta[idx:])
    maelower = np.mean(delta[:idx])
    return maelower, maeupper

def fitGPD(stnnum, stnname, inputPath, outputPath, jitter=False):
    """
    Fit a Gamma-GPD mixture model to the gust data for a given station
    It's known that the data has an issue with rounded values. The data are stored in the files in
    units of m/s to 1 decimal place, but it's understood the data are recorded in knots, and rounded to whole values
    prior to storing in the database. Use the `jitter` kwarg to add a uniform random variation to the values.

    :param int stnnum: BoM Station Number
    :param str stnname: BoM Station Name - for plot titles
    :param str inputPath: Directory where observation data files are stored
    :param str outputPath: Directory where output should be stored
    :param bool jitter: Add random variation around the values. Default `False`
    """
    gust = loadData(inputPath, stnnum)
    if len(gust) < 365.25 * 5:
        print(f"Insufficient data for {stnnum}")
        return
    jgust = gust
    if jitter:
        jgust += np.random.uniform(-0.2572, 0.2572, len(gust))
    tlower = np.quantile(gust, 0.75)
    tupper = np.max(gust)
    paramdf = pd.DataFrame(columns=['t', 'lnmean', 'lnsd', 'u', 'sigmau', 'xi', 'phiu', 'sigmase', 'xise'])
    for t in np.arange(tlower, tupper, 0.05):
        if len(gust[gust > t]) < 10:
            print("Less than 10 records greater than the threshold")
            break
        fit = evmix.flognormgpd(jgust, useq=t, std_err=True)
        if fit.rx2('xi')[0] < 0:
            if isinstance(fit.rx2('se'), np.ndarray):
                sigmase = fit.rx2('se')[3]
                xise = fit.rx2('se')[4]
            else:
                sigmase = 0
                xise = 0
            #maelower, maeupper = calcMAE(gust, t, fit)
            paramdf = paramdf.append({"t": t, "lnmean": fit.rx2('lnmean')[0], "lnsd": fit.rx2('lnsd')[0],
                            "u": fit.rx2('u')[0], 'sigmau': fit.rx2('sigmau')[0],
                            "xi": fit.rx2('xi')[0], "phiu": fit.rx2('phiu')[0],
                            "sigmase": sigmase, "xise": xise},
                            ignore_index=True)
    rlplotfile = os.path.join(outputPath, f"{stnnum:06d}.rl.png")
    epplotfile = os.path.join(outputPath, f"{stnnum:06d}.ep.png")
    histplotfile = os.path.join(outputPath, f"{stnnum:06d}.hist.png")
    if len(paramdf > 0):
        rlplot(gust, paramdf, stnname, rlplotfile)
        epplot(gust, paramdf, stnname, epplotfile)
        histplot(jgust, paramdf, stnname, histplotfile)

        paramdf.to_csv(os.path.join(outputPath, f"{stnnum:06d}.csv"), index=False)

def start():
    p = argparse.ArgumentParser()
    p.add_argument('-c', '--config_file', help="Configuration file")
    p.add_argument('-v', '--verbose', help="Verbose output",
                   action='store_true')
    args = p.parse_args()

    global configFile

    configFile = args.config_file
    verbose = args.verbose
    config = ConfigParser(allow_no_value=True,
                          interpolation=ExtendedInterpolation())
    config.optionxform = str
    config.read(configFile)

    logFile = config.get('Logging', 'LogFile')
    logdir = dirname(realpath(logFile))
    # if log file directory does not exist, create it
    if not isdir(logdir):
        try:
            os.makedirs(logdir)
        except OSError:
            logFile = pjoin(os.getcwd(), 'evmix.log')

    verbose = config.getboolean('Logging', 'Verbose', fallback=False)

    if args.verbose:
        verbose = True

    pInit(configFile)
    main(config, verbose)

def main(config, verbose=False):
    logFile = config.get('Logging', 'LogFile')
    logLevel = config.get('Logging', 'LogLevel', fallback='INFO')
    verbose = config.getboolean('Logging', 'Verbose', fallback=verbose)
    datestamp = config.getboolean('Logging', 'Datestamp', fallback=False)
    LOGGER = flStartLog(logFile, logLevel, verbose, datestamp)

    ListAllFiles(config)
    processFiles(config)

# -----------------------------------------------------------------
inputPath = r"X:\georisk\HaRIA_B_Wind\data\raw\from_bom\2019\Daily"
outputPath = r"X:\georisk\HaRIA_B_Wind\data\derived\obs\daily_max_wind\wind\lognormgpd_iterative"
stationfile = os.path.join(inputPath,'DC02D_StnDet_999999999632559_updated.txt')
stndf = pd.read_csv(stationfile)

for idx, stn in stndf.iterrows():
    stnnum = stn['Bureau of Meteorology Station Number']
    stnname = stn['Station Name'].strip()
    print(f"Processing {stnname} ({stnnum})")
    fitGPD(stnnum, stnname, inputPath, outputPath)

