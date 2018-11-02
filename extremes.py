"""
:mod:`extremes` -- routines for evaluating extreme value distributions
========================================================================

.. module:: extremes
    :synopsis: Offer functions for calculating return periods or return
               levels using data and/or a fitted distribution.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""
from __future__ import division, print_function
from itertools import product
import logging
from functools import wraps
import time

import numpy as np
from distributions import empiricalPDF
from scipy.stats import genpareto, scoreatpercentile
from scipy.optimize import curve_fit
import lmfit

LOG = logging.getLogger(__name__)

OBS_PER_YEAR = 365.25

__all__ = ['returnLevels', 'empReturnPeriod', 'nearestIndex',
           'gpdCalculateShape', 'gpdSelectThreshold', 'gpdAsymptote',
           'returnPeriodUncertainty']

def timer(f):
    """
    A simple timing decorator for the entire process.

    """
    @wraps(f)
    def wrap(*args, **kwargs):
        t1 = time.time()
        res = f(*args, **kwargs)

        tottime = time.time() - t1
        msg = "%02d:%02d:%02d " % \
          reduce(lambda ll, b : divmod(ll[0], b) + ll[1:],
                        [(tottime,), 60, 60])

        LOG.info("Time for {0}: {1}".format(f.func_name, msg) )
        return res

    return wrap

def returnLevels(intervals, mu, xi, sigma, rate, npyr=OBS_PER_YEAR):
    """
    Calculate return levels for specified intervals for a generalised pareto
    distribution with the given threshold, scale and shape parameters.

    :param intervals: :class:`numpy.ndarray` or float of recurrence intervals
              to evaluate return levels for.
    :param float mu: Threshold parameter (also called location).
    :param float xi: Shape parameter.
    :param float sigma: Scale parameter.
    :param float rate: Rate of exceedances (i.e. number of observations greater
                       than `mu`, divided by total number of observations).
    :param float npyr: Number of observations per year.

    :returns: return levels for the specified recurrence intervals.

    """

    rp = mu + (sigma / xi) * (np.power(intervals * npyr * rate, xi) - 1.)
    return rp


def empReturnPeriod(data, npyr=OBS_PER_YEAR):
    """
    Returns the empirically-based recurrence interval (in years) for a set
    of observations.

    It is assumed the data are daily observations. If the observations are not
    daily, there are two options: set the ``npyr`` variable, or backfill the
    ``data`` variable with zero values to match the assumed length of the
    record.

    The highest return period should be (approximately) len(``data``)/``npyr``.

    :param data: :class:`numpy.ndarray` containing the observed values (with
                 missing values removed).
    :param float npy: Number of observations per year (default=365.25)

    :returns: Recurrence intervals for the observed data.
    :rtype: :class:`numpy.ndarray`
    """
    nobs = len(data)
    # Empirical return periods:
    emprp = 1. / (1. - np.arange(1, nobs + 1, 1) / (nobs + 1)) / npyr
    return emprp


def returnPeriodUncertainty(data, mu, xi, sigma, intervals):
    """
    Calculate uncertainty around a fit, holding threshold fixed.

    :param data: :class:`numpy.ndarray` containing the observed values (with
                 missing values removed).
    :param float mu: Threshold parameter (also called location).
    :param float xi: Shape parameter.
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

    # We need an empirical PDF to serve as the dependant data in
    # the fitting routine.
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

def returnPeriodConfidence(data, mu, xi, sigma, intervals):
    """
    Calculate confidence intervals around a fit, using bootstrap resampling.

    :param data: :class:`numpy.ndarray` containing the observed values (with
                 missing values removed).
    :param float mu: Threshold parameter (also called location).
    :param float xi: Shape parameter.
    :param float sigma: Scale parameter.
    :param intervals: :class:`numpy.ndarray` or float of recurrence intervals
              to evaluate return level uncertainties for.

    :returns: 5th and 95th percentile values of return level for the given
              recurrence intervals.
    :rtype: :class:`numpy.ndarray`

    """

    sortedmax = np.sort(data[data > mu])
    nobs = len(sortedmax)
    rate = float(nobs) / float(len(data))
    rpvals = np.zeros((len(intervals), 1000))
    for i in range(1000):
        sample = np.random.choice(sortedmax, size=len(sortedmax)-1)
        emppdf = empiricalPDF(sample)
        try:
            popt, pcov = curve_fit(lambda x, mu, xi, sigma:
                                   genpareto.pdf(x, xi, loc=mu, scale=sigma),
                                   sample, emppdf, (mu, xi, sigma), maxfev=5000)
        except RuntimeError:
            pass
        else:
            rpvals[:, i] = returnLevels(intervals, popt[0], popt[1],
                                        popt[2], rate)

    rp95CI = np.percentile(rpvals, 95, axis=1)
    rp05CI = np.percentile(rpvals, 5, axis=1)

    return (rp05CI, rp95CI)

def nearestIndex(array, value):
    """
    Determine the index of the array element closest to the given value.

    :param array: input array containing values to check against.
    :type array: :class:`numpy.ndarray`
    :param value: value to find closest element of `array`.
    :type value: float or int

    :returns: index of `array` that has the closest value to `value`.
    :rtype: int
    """

    idx = (np.abs(array - value)).argmin()
    return idx

def gpdCalculateShape(mu, data):
    """
    Calculate the shape (and scale) parameter for a Generalised
    Pareto Distribution (GPD), given a threshold parameter.

    :param float mu: threshold parameter for the GPD distribution.
    :param data: :class:`numpy.ndarray` of data values to fit.

    :returns: tuple of GPD parameters
              (threshold, shape, scale | mu, xi, sigma)
    """
    nobs = len(data)
    nexc = len(data[data > mu])
    rate = float(nexc)/float(nobs)
    gpd = genpareto.fit(data[data > mu] - mu)

    return gpd

@timer
def gpdSelectThreshold(data, nexc=10):
    """
    Select an appropriate threshold for fitting a Generalised Pareto
    Distribution, using the approach described by Sanabria and Cechet (2007).

    The only constraint placed on the selection is that the shape
    parameter is negative (such that the distribution is bounded).

    :param data: :class:`numpy.ndarray` containing data values (with
                 missing values removed).
    :param int nexc: Minimum number of data points exceeding the threshold for
                     a fit to be performed.
    :returns: tuple of the shape, scale and threshold.

    References: Sanabria, L. A. and Cechet, R. P. (2007), *A Statistical
                Model of Severe Winds*. Geoscience Australia Record 2007/12.
    """

    sh = []
    sc = []
    t = []
    q1000list = []
    q10000list = []

    eps = -0.01
    datamax = data.max()
    nobs = len(data)
    if len(data.compress(data > 0)) == 0:
        LOG.warn("All data for threshold selection are zero - exiting")
        return 0, 0, 0
    startmu = data.compress(data > 0).max()/2
    for mu in np.arange(startmu, datamax, 0.05):
        numexceed = len(data[data > mu])
        rate = numexceed / nobs
        if numexceed < nexc:
            break

        pp = gpdCalculateShape(mu, data)
        q1000, q10000 = returnLevels(np.array([1000, 10000]),
                                     mu, pp[0], pp[2], rate)
        if np.isnan(q1000):
            continue

        if np.isnan(q10000):
            continue

        qdiff = np.abs(q10000 - q1000)

        # If the shape parameter is negative, the difference between
        # the 1000- and 10000-year return period values is less than
        # 12%, and the difference is positive, then we store
        # the threshold, shape and scale parameters, as well as
        # the 1000- and 10000-year return period values.
        # the 12% difference is to avoid numerical instability issues
        # in the fitted distribution for large return intervals.
        if pp[0] < eps and qdiff < 0.12*q10000 and qdiff > -eps:
            t.append(mu)
            sh.append(pp[0])
            sc.append(pp[2])
            q1000list.append(q1000)
            q10000list.append(q10000)

    # If there are no valid threshold values, then return zeros:
    if len(t) == 0:
        LOG.warn("No suitable shape parameters identified")
        return 0, 0, 0

    # Take the average value of the 1000- and 10000-year return period
    # values, then add 5% and find the integer value above that.
    # This gives the optimal threshold value, taking the higher threshold
    # when viewing the threshold for the different quantiles.
    Av1000 = np.mean(np.array(q1000list))
    Av10000 = np.mean(np.array(q10000list))
    Av1000 = np.ceil(Av1000 + 0.05*Av1000)
    Av10000 = np.ceil(Av10000 + 0.05*Av10000)

    idx1000 = nearestIndex(np.array(q1000list), Av1000)
    idx10000 = nearestIndex(np.array(q10000list), Av10000)

    u1000 = t[idx1000]
    u10000 = t[idx10000]

    if u1000 > u10000:
        shmax = sh[idx1000]
        scmax = sc[idx1000]
    else:
        shmax = sh[idx10000]
        scmax = sc[idx10000]

    return shmax, scmax, u1000

def gpdAsymptote(mu, xi, sigma):
    """
    Calculate the limiting value of a (bounded) GPD, based on the fitted
    values of mu, sigma and xi (as demonstrated in Coles,2001). Note xi is
    negative, so we take the absolute value of the relation.

    :param float mu: Location parameter of the fitted GPD.
    :param float xi: Shape parameter of the fitted GPD.
    :param float sigma: Scale parameter of the fitted GPD.

    :returns: asymptote value
    :rtype: float

    :raises ValueError: If sigma >= 0, which defines an unbounded GPD.
    """

    if sigma >= 0.0:
        raise ValueError("Shape parameter value defines unbounded distribution")
    limit = np.abs((mu - xi) / sigma)

    return limit

def residual(p, x, y):
    return genpareto.pdf(x, p['xi'], loc=p['mu'], scale=p['sig']) - y

def calcSD(pars, cov, rate, npyr, intervals):
    nsims = 1000
    rps = np.zeros((nsims, len(intervals)))
    for i in range(nsims):
        try:
            xi = pars[0] + np.random.normal(0, np.sqrt(cov[0,0]))
            mu = pars[1] + np.random.normal(0, np.sqrt(cov[1,1]))
            sig = pars[2] + np.random.normal(0, np.sqrt(cov[2,2]))
            rps[i, : ] = mu + (sig / xi) * (np.power(intervals * npyr * rate, xi) - 1.)
        except ValueError:
            rps[i, : ] = np.zeros(len(intervals))

    lowerrp = scoreatpercentile(rps, 5, axis=0)
    upperrp = scoreatpercentile(rps, 95, axis=0)
    return lowerrp, upperrp

def calculateUncertainty(wspd, intervals, xi, mu, sig):
    """
    :param wspd: :class:`numpy.array` of wind speeds
    :param float xi: initial guess for 
    """
    bins = np.arange(0.5, 100, 1)
    n, bins = np.histogram(wspd, bins, normed=True)
    centres = 0.5*(bins[1:]+bins[:-1])
    try:
        pars,cov = curve_fit(lambda x, xi, mu, sig: genpareto.pdf(x, xi, loc=mu, 
                                                                  scale=sig), 
                             centres, n, p0=[0, np.mean(wspd), np.mean(wspd)], 
                             maxfev=10000 )
    except RuntimeError as e:
        LOG.warn(e)
        return None, None, None

    gpd = genpareto.fit(wspd, floc=mu)
    npyr = 365.25

    p = lmfit.Parameters()
    p.add_many(('xi', gpd[0], True, -np.inf, 2.),
               ('mu', gpd[1]),
               ('sig', gpd[2]))

    mini = lmfit.Minimizer(residual, p, fcn_args=(centres, n),
                           nan_policy='omit')

    # first solve with Nelder-Mead
    out1 = mini.minimize(method='Nelder')
    out2 = mini.minimize(method='leastsq', params=out1.params)

    if hasattr(out2, 'covar'):

        cmu = out2.params['mu'].value
        cxi = out2.params['xi'].value
        csig = out2.params['sig'].value
        rate = len(wspd)/(npyr*10000)
        crp = cmu + (csig / cxi) * (np.power(intervals * npyr * rate, cxi) - 1.)
        lrp, urp = calcSD((cxi, cmu, csig), out2.covar, rate, npyr, intervals)
        return (cxi, cmu, csig), crp, lrp, urp

    else:
        LOG.warn("No covariance matrix from the minimizer")
        return None, None, None
