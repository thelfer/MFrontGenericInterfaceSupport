# -*- coding: utf-8 -*-

import os
try:
    import unittest2 as unittest
except ImportError:
    import unittest
import mgis.behaviour as mgis_bv


class ParameterTest(unittest.TestCase):
    def test_pass(self):
        yg = 150e9
        nu = 0.3
        eps = 1.e-14
        # path to the test library
        lib = os.environ['MGIS_TEST_BEHAVIOURS_LIBRARY']
        version = os.environ['MGIS_TEST_TFEL_VERSION']
        h = mgis_bv.Hypothesis.Tridimensional
        b = mgis_bv.load(lib, 'ParameterTest', h)
        self.assertTrue(b.behaviour == "ParameterTest",
                        "invalid behaviour name")
        self.assertTrue(b.hypothesis == h, "invalid hypothesis")
        self.assertTrue(b.source == "ParameterTest.mfront", "invalid source")
        self.assertTrue(b.tfel_version == version, "invalid TFEL version")
        self.assertTrue(len(b.params) == 6, "invalid number of parameters")
        self.assertTrue(b.params[0] == "YoungModulus",
                        "invalid first parameter")
        self.assertTrue(b.params[1] == "PoissonRatio",
                        "invalid second parameter")
        self.assertTrue(b.params[2] == "ParametersArray[0]",
                        "invalid third parameter")
        self.assertTrue(b.params[3] == "ParametersArray[1]",
                        "invalid fourth parameter")
        self.assertTrue(b.params[4] == "minimal_time_step_scaling_factor",
                        "invalid fifth parameter")
        self.assertTrue(b.params[5] == "maximal_time_step_scaling_factor",
                        "invalid sixth parameter")
        yg_v = mgis_bv.getParameterDefaultValue(b, "YoungModulus")
        nu_v = mgis_bv.getParameterDefaultValue(b, "PoissonRatio")
        self.assertTrue(
            abs(yg_v - yg) < eps * yg, "invalid 'YoungModulus' default value")
        self.assertTrue(
            abs(nu_v - nu) < eps * nu, "invalid 'PoissonRatio' default value")


if __name__ == '__main__':
    unittest.main()
