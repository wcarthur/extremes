import unittest
import numpy as np

from extremes import gpdSelectThreshold, gpdAsymptote

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
        """Test that len(data)<nexc gives zeros as result"""
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
