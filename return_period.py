"""
:mod:`return_period` -- return levels for empirical/fitted distributions
========================================================================

.. module:: return_period
    :synopsis: Offer functions for calculating return periods or return
               levels using data and/or a fitted distribution.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""
from __future__ import division
import numpy as np
from distributions import empiricalPDF
from itertools import product
from scipy.stats import genpareto
from scipy.optimize import curve_fit

import logging
LOG = logging.getLogger(__name__)

OBS_PER_YEAR = 365.25


def returnLevels(intervals, mu, xi, sigma, rate, npyr=OBS_PER_YEAR):
    """
    Calculate return levels for specified intervals for a distribution with
    the given threshold, scale and shape parameters.

    :param intervals: :class:`numpy.ndarray` or float of return period intervals
              to evaluate return levels for.
    :param float mu: Threshold parameter (also called location).
    :param float chi: Shape parameter.
    :param float sigma: Scale parameter.
    :param float rate: Rate of exceedances (i.e. number of observations greater
                       than `mu`, divided by total number of observations).
    :param float npyr: Number of observations per year.

    :returns: return levels for the specified return intervals.

    """

    rp = mu + (sigma / xi) * (np.power(intervals * npyr * rate, xi) - 1.)
    return rp


def empiricalReturnPeriod(data, npyr=OBS_PER_YEAR):
    """
    Returns the empirically-based return interval (in years) for a set
    of observations.

    It is assumed the data are daily observations.

    The highest return period should be (approximately) len(``data``)/``npy``.

    :param data: :class:`numpy.ndarray` containing the observed values (with
                 missing values removed).
    :param float npy: Number of observations per year (default=365.25)

    :returns: Return periods for the observed data.
    :rtype: :class:`numpy.ndarray`
    """
    nobs = len(data)
    # Empirical return periods:
    emprp = 1. / (1. - np.arange(1, nobs + 1, 1) / (nobs + 1)) / npyr
    return emprp


def returnPeriodUncertainty(data, mu, xi, sigma, intervals):
    """
    Calculate uncertainty around a fit, holding threshold fixed

    :param data: :class:`numpy.ndarray` containing the observed values (with
                 missing values removed).
    :param float mu: Threshold parameter (also called location).
    :param float chi: Shape parameter.
    :param float sigma: Scale parameter.
    :param intervals: :class:`numpy.ndarray` or float of return period intervals
              to evaluate return level uncertainties for. 

    :returns: Array of standard deviation values for each return period, based
              on all permutations of shape and scale parameters with standard
              errors.
    :rtype: :class:`numpy.ndarray`

    """
    sortedmax = np.sort(data[data > mu])
    nobs = len(sortedmax)
    rate = float(nobs) / float(len(data))
    emppdf = empiricalPDF(data[data > mu])

    # Perform the curve fitting, holding ``mu`` fixed and allowing
    # ``xi`` and ``sigma`` to vary.
    try:
        popt, pcov = curve_fit(lambda x, xi, sigma:
                               genpareto.pdf(x, xi, loc=mu, scale=sigma),
                               sortedmax, emppdf, (xi, sigma))
    except RuntimeError as e:
        LOG.exception("Curve fitting failed: %s", e)
        return np.zeros(len(intervals))

    sd = np.sqrt(np.diag(pcov))

    svals = (sigma - sd[1], sigma, sigma + sd[1])
    mvals = (mu, mu, mu)
    xvals = (xi - sd[0], xi, xi + sd[0])

    rpvalues = np.array([returnLevels(intervals, m, xii, s, rate) for
                         (s, m, xii) in product(svals, mvals, xvals)])

    rpFitError = np.std(rpvalues, axis=0)

    return rpFitError
