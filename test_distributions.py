import unittest
import numpy as np

from distributions import asymptote

class TestAsymptote(unittest.TestCase):

    def setUp(self):
        self.mus = np.array([])
        self.xis = np.array([])
        self.sigmas = np.array([])
        self.asyms = np.array([])

    def test_singledist(self):
        """Test function with single distribution parameters"""
        self.assertAlmostEqual(asymptote(27.43, 5.3089, -0.18813), 117.5841, 4)
        self.assertAlmostEqual(asymptote(30.92, 4.5458, -0.08184), 322.2654, 4)
        self.assertAlmostEqual(asymptote(18.74, 2.7587, -0.01000), 1598.130, 4)

    def test_unbounded(self):
        """Raise ValueError if the fitted distribution is unbounded"""
        self.assertRaises((ValueError,), asymptote, 20, 0.01, 6)
        self.assertRaises((ValueError,), asymptote, 20, 0.00, 6)

    def test_arrayinput(self):
        pass
    
if __name__ == '__main__':
    unittest.main()
