"""
:mod:`plot` -- basic plotting routines for extremes
===================================================

.. module:: plot
    :synopsis: Provide some basic plotting routines for extreme value
               distribution modelling.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""

from __future__ import division

import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import genpareto
import statsmodels.api as sm

import logging
LOG = logging.getLogger(__name__)

from extremes import returnLevels, empReturnPeriod, returnPeriodUncertainty
from distributions import fittedPDF

sns.set_context("poster")
sns.set_style("ticks")


def plotFit(data, mu, xi, sigma, title, figfile):
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
    LOG.info("Plotting fitted return period curve")

    rp = np.array([1, 2, 5, 10, 20, 50, 100, 200,
                   500, 1000, 2000, 5000, 10000])
    rate = float(len(data[data > mu])) / float(len(data))
    rval = returnLevels(rp, mu, xi, sigma, rate)

    emprp = empReturnPeriod(data)
    err = returnPeriodUncertainty(data, mu, xi, sigma, rp)

    sortedmax = np.sort(data)
    fig, ax1 = plt.subplots(1, 1, figsize=(12, 12))
    ax1.semilogx(rp, rval, label="Fitted RP curve")
    ax1.semilogx(rp, rval + 1.96 * err, label="95% CI",
                 linestyle='--', color='0.5')
    ax1.semilogx(rp, rval - 1.96 * err, linestyle='--', color='0.5')

    ax1.scatter(emprp[emprp > 1], sortedmax[emprp > 1], s=100,
                color='r', label="Empirical RP")

    title_str = (title + "\n" +
                 r"$\mu$ = {0:.2f}, $\xi$ = {1:.5f}, $\sigma$ = {2:.4f}, $\rho$ = {3:.4f}".
                 format(mu, xi, sigma, rate))
    ax1.set_title(title_str)
    ax1.legend(loc=2)
    ax1.set_ylim((0, 100))
    ax1.set_xlim((1, 10000))
    ax1.set_ylabel('Wind speed (m/s)')
    ax1.set_xlabel('Return period (years)')
    ax1.grid(which='major')
    ax1.grid(which='minor', linestyle='--', linewidth=1)
    ax1.axhline(45.6, c='lime', linestyle='--', linewidth=2)
    ax1.axhline(62.5, c='darkorange', linestyle='--', linewidth=2)
    ax1.axhline(77.8, c='darkred', linestyle='--', linewidth=2)
    ax1.text(20000, 45.6, 'Cat 3', ha='center')
    ax1.text(20000, 62.5, 'Cat 4', ha='center')
    ax1.text(20000, 77.8, 'Cat 5', ha='center')
    plt.savefig(figfile)
    plt.close()


def plotDiagnostics(data, mu, xi, sigma, figfile):
    """
    Create a 4-panel diagnostics plot of the fitted distribution.

    :param data: :class:`numpy.ndarray` of observed data values (in units
                 of metres/second).
    :param float mu: Selected threshold value.
    :param float xi: Fitted shape parameter.
    :param float sigma: Fitted scale parameter.
    :param str figfile: Path to store the file (includes image format)

    """
    LOG.info("Plotting diagnostics")
    fig, ax = plt.subplots(2, 2)
    axes = ax.flatten()
    # Probability plots
    sortedmax = np.sort(data[data > mu])
    gpdf = fittedPDF(data, mu, xi, sigma)
    pp_x = sm.ProbPlot(sortedmax)
    pp_x.ppplot(xlabel="Empirical", ylabel="Model", ax=axes[0], line='45')
    axes[0].set_title("Probability plot")

    prplot = sm.ProbPlot(sortedmax, genpareto, distargs=(xi,),
                         loc=mu, scale=sigma)
    prplot.qqplot(xlabel="Model", ylabel="Empirical", ax=axes[1], line='45')
    axes[1].set_title("Quantile plot")

    ax2 = axes[2]
    rp = np.array([1, 2, 5, 10, 20, 50, 100, 200,
                   500, 1000, 2000, 5000, 10000])
    rate = float(len(sortedmax)) / float(len(data))
    rval = returnLevels(rp, mu, xi, sigma, rate)

    emprp = empReturnPeriod(np.sort(data))
    ax2.semilogx(rp, rval, label="Fitted RP curve", color='r')
    ax2.scatter(emprp[emprp > 1], np.sort(data)[emprp > 1],
                color='b', label="Empirical RP", s=100)
    ax2.legend(loc=2)
    ax2.set_xlabel("Return period")
    ax2.set_ylabel("Return level")
    ax2.set_title("Return level plot")
    ax2.grid(True)
    maxbin = 4 * np.ceil(np.floor(data.max() / 4) + 1)
    sns.distplot(sortedmax, bins=np.arange(mu, maxbin, 2),
                 hist=True, axlabel='Wind speed (m/s)',
                 kde_kws={"label":"Empirical PDF"},
                 ax=axes[3])
    axes[3].plot(sortedmax, gpdf, color='r', label='Fitted PDF')
    axes[3].set_title("Density plot")
    axes[3].legend(loc=1)
    plt.tight_layout()
    plt.savefig(figfile)
    plt.close()
