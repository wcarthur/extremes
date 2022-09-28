import unittest
import numpy as np
from NumpyTestCase import NumpyTestCase

from extremes import gpdSelectThreshold, gpdAsymptote, ppoints

class TestPoints(NumpyTestCase):
    """
    Test ppoints function
    """
    def setUp(self):
        self.input = np.array(np.arange(1, 21))
        self.results = np.array([0.06097561, 0.15853659, 0.25609756,
                                 0.35365854, 0.45121951, 0.54878049,
                                 0.64634146, 0.74390244, 0.84146341,
                                 0.93902439])
        self.results1 = np.array([0.1, 0.26, 0.42, 0.58, 0.74, 0.9])
        self.results2 = np.array([0.025, 0.075, 0.125, 0.175, 0.225,
                                  0.275, 0.325, 0.375, 0.425,  0.475,
                                  0.525, 0.575, 0.625, 0.675, 0.725,
                                  0.775, 0.825, 0.875, 0.925, 0.975])

    def test_ppoints(self):
        self.numpyAssertAlmostEqual(ppoints(10), self.results)
        self.numpyAssertAlmostEqual(ppoints(6), self.results1)
        self.numpyAssertAlmostEqual(ppoints(self.input), self.results2)

class TestSelection(unittest.TestCase):
    """
    Test gpdSelectThreshold function
    """

    def setUp(self):
        self.emptyResult = (0, 0, 0)

    def test_allzeros(self):
        """Test all zero input gives zeros as result"""
        data = np.zeros(50)
        result = gpdSelectThreshold(data)

        self.assertEqual(type(result), type(self.emptyResult))
        self.assertEqual(len(result), len(self.emptyResult))
        self.assertTrue(np.alltrue(np.equal(result, self.emptyResult)))

    def test_nexceed(self):
        """Test that len(data) < nexc gives zeros as result"""
        data = np.ones(10)
        result = gpdSelectThreshold(data, nexc=20)
        self.assertEqual(type(result), type(self.emptyResult))
        self.assertEqual(len(result), len(self.emptyResult))
        self.assertTrue(np.alltrue(np.equal(result, self.emptyResult)))

class TestAsymptote(unittest.TestCase):

    def setUp(self):
        self.mus = np.array([])
        self.xis = np.array([])
        self.sigmas = np.array([])
        self.asyms = np.array([])

    def test_singledist(self):
        """Test function with single distribution parameters"""
        self.assertAlmostEqual(gpdAsymptote(27.43, 5.3089, -0.18813), 117.5841, 4)
        self.assertAlmostEqual(gpdAsymptote(30.92, 4.5458, -0.08184), 322.2654, 4)
        self.assertAlmostEqual(gpdAsymptote(18.74, 2.7587, -0.01000), 1598.130, 4)

    def test_unbounded(self):
        """Raise ValueError if the fitted distribution is unbounded"""
        self.assertRaises((ValueError,), gpdAsymptote, 20, 0.01, 6)
        self.assertRaises((ValueError,), gpdAsymptote, 20, 0.00, 6)

    def test_arrayinput(self):
        pass

if __name__ == '__main__':
    unittest.main()
