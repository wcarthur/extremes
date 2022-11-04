import os
import pandas as pd
import numpy as np
import seaborn as sns
import scipy.stats as stats
from matplotlib import pyplot as plt

from extremes import returnLevels, empReturnPeriod, returnPeriodUncertainty

sns.set_context("talk")
figsize = (16, 9)

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
    gpd = stats.genpareto.fit(data[data > mu] - mu)

    return gpd

def selectThreshold(data, start=None):
    """
    Select an appropriate threshold for fitting a generalised pareto
    distribution.

    The only constraint placed on the selection is that the shape
    parameter is negative (such that the distribution is bounded).

    :param data: :class:`numpy.ndarray` containing the observed values (with
                 missing values removed).
    :param start: starting point for the threshold value. If not given,
                  defaults to the median of the ``data`` variable.
    :returns: tuple of the shape, scale and threshold.
    """

    sh = []
    sc = []
    t = []
    q1000list = []
    q10000list = []

    eps = -0.01
    datamax = data.max()
    nobs = len(data)
    if start:
        startValue = start
    else:
        startValue = np.median(data)
    for mu in np.arange(startValue, datamax, 0.01):
        nexc = len(data[data > mu]) 
        rate = nexc / nobs
        if nexc < 5:
            break

        pp = calculateShape(mu, data)
        q1000, q10000 = returnLevels(np.array([1000, 10000]), mu, pp[0], pp[2], rate)
        if np.isnan(q1000):
            continue

        if np.isnan(q10000):
            continue

        qdiff = np.abs(q10000 - q1000)
        if pp[0] < eps and qdiff < 0.2*q10000 and qdiff > -eps:
            t.append(mu)
            sh.append(pp[0])
            sc.append(pp[2])
            q1000list.append(q1000)
            q10000list.append(q10000)

    if len(t) == 0:
        #print "No suitable shape parameters identified"
        return 0, 0, 0
    av1000 = np.mean(np.array(q1000list))
    av10000 = np.mean(np.array(q10000list))
    av1000 = np.ceil(av1000 + 0.05*av1000)
    av10000 = np.ceil(av10000 + 0.05*av10000)

    idx1000 = find_nearest_index(np.array(q1000list), av1000)
    idx10000 = find_nearest_index(np.array(q10000list), av10000)

    u1000 = t[idx1000]
    u10000 = t[idx10000]

    if u1000 > u10000:
        shmax = sh[idx1000]
        scmax = sc[idx1000]
    else:
        shmax = sh[idx10000]
        scmax = sc[idx10000]

    return shmax, scmax, u1000

def plotFit(data, mu, xi, sigma, title):
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
    rate = float(len(data[data > mu])) / float(len(data))
    rval = returnLevels(rp, mu, xi, sigma, rate)

    emprp = empReturnPeriod(data)
    err = returnPeriodUncertainty(data, mu, xi, sigma, rp)

    sortedmax = np.sort(data)
    fig, ax1 = plt.subplots(1, 1, figsize=figsize)
    ax1.semilogx(rp, rval, label="Fitted RP curve")
    ax1.semilogx(rp, rval + 1.96 * err, label="95% CI",
                 linestyle='--', color='0.5')
    ax1.semilogx(rp, rval - 1.96 * err, linestyle='--', color='0.5')

    ax1.scatter(emprp[emprp > 1], sortedmax[emprp > 1], marker='x', s=100,
                color='k', label="Empirical RP", zorder=100)

    title_str = (title + "\n" +
                 r"$\mu$ = {0:.3f}, $\xi$ = {1:.5f}, $\sigma$ = {2:.4f}".
                 format(mu, xi, sigma))
    ax1.set_title(title_str)
    ax1.legend(loc=2)
    ax1.set_ylim((0, 100))
    ax1.set_xlim((1, 10000))
    ax1.set_ylabel('Wind speed (m/s)')
    ax1.set_xlabel('Return period (years)')
    ax1.grid(which='major')
    ax1.grid(which='minor', linestyle='--', linewidth=1)



datapath = "X:/georisk/HaRIA_B_Wind/projects/tcha/data/derived/observations/daily"
station = 14015
filename = os.path.join(datapath, f"bom_{station:06d}.csv")
df = pd.read_csv(filename, index_col=0, parse_dates=[2], infer_datetime_format=True )
timerange = df.datetime.max().year - df.datetime.min().year + 1
xi, sigma, mu = selectThreshold(1.098*df.gust) #, start=np.min(df.gust))
dummydata = np.zeros(int(timerange)*365)
ndata = len(df.gust)
dummydata[-ndata:] = 1.098*df.gust.values

plotFit(dummydata, mu, xi, sigma, "Darwin Airport")
plt.savefig(os.path.join(datapath, f"{station:06d}.fit.png"), dpi=300, bbox_inches='tight')
