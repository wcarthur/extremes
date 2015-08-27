"""
:mod:`distributions` -- calculate empirical/fitted distributions for obs
========================================================================

.. module:: distributions
    :synopsis: Calculate empirical and/or fitted distributions for data

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""

import numpy as np
from scipy.stats import genpareto
from statsmodels.nonparametric.kde import KDEUnivariate

import logging
LOG = logging.getLogger(__name__)


def empiricalPDF(data):
    """
    Evaluate a probability density function using kernel density
    estimation for input data.

    :param data: :class:`numpy.ndarray` of data values.

    :returns: PDF values at the data points.
    """
    LOG.debug("Calculating empirical PDF")
    sortedmax = np.sort(data)
    kde = KDEUnivariate(sortedmax)
    kde.fit()
    try:
        res = kde.evaluate(sortedmax)
    except MemoryError:
        res = np.zeros(len(sortedmax))
    return res


def fittedPDF(data, mu, sigma, xi):
    """
    Calculate probability denisty function values given data and
    GPD fit parameters.

    :param data: :class:`numpy.ndarray` of data values.
    :param float mu: Location parameter of the fitted GPD.
    :param float sigma: Shape parameter of the fitted GPD.
    :param float xi: Scale parameter of the fitted GPD.

    :returns: :class:`numpy.ndarray` of PDF values at the data points.

    """

    LOG.debug("Calculating fitted GPD PDF")
    res = genpareto.pdf(np.sort(data[data > mu]),
                        sigma, loc=mu, scale=xi)

    return res


def generateDistributions(data, mu, sigma, xi):
    """
    Generate empirical and fitted PDF values for selected data, based on
    threshold, shape and scale parameters.

    :param data: :class:`numpy.ndarray` of data values.
    :param float mu: Location parameter of the fitted GPD.
    :param float sigma: Shape parameter of the fitted GPD.
    :param float xi: Scale parameter of the fitted GPD.

    """

    if mu > data.max():
        raise ValueError("Threshold greater than maximum data value")

    emppdf = empiricalPDF(data[data > mu], mu)
    gpdf = fittedPDF(data, mu, sigma, xi)

    return emppdf, gpdf
