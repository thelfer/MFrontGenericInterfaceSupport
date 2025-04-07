# -*- coding: utf-8 -*-

import os
try:
    import unittest2 as unittest
except ImportError:
    import unittest
import mgis.behaviour as mgis_bv


class BoundsCheckTest(unittest.TestCase):
    def test_pass(self):

        eps = 1.e-14
        yg_min = 100.e9
        yg_max = 200.e9
        nu_min = -1
        nu_max = 0.5
        iv_pmin = 0.
        iv_pmax = 0.8
        iv_min = 0.2
        iv_max = 0.5
        ev_pmin = 0.
        ev_pmax = 500.
        ev_min = 200.
        ev_max = 400.
        # path to the test library
        lib = os.environ['MGIS_TEST_BEHAVIOURS_LIBRARY']
        version = os.environ['MGIS_TEST_TFEL_VERSION']
        h = mgis_bv.Hypothesis.Tridimensional
        b = mgis_bv.load(lib, 'BoundsCheckTest', h)
        self.assertTrue(b.behaviour == "BoundsCheckTest",
                        "invalid behaviour name")
        self.assertTrue(b.hypothesis == h, "invalid hypothesis")
        self.assertTrue(b.source == "BoundsCheckTest.mfront", "invalid source")
        self.assertTrue(b.tfel_version == version, "invalid TFEL version")
        # test on material properties
        self.assertTrue(mgis_bv.hasBounds(b, "YoungModulus"),
                        "'YoungModulus' shall have bounds")
        self.assertTrue(mgis_bv.hasLowerBound(b, "YoungModulus"),
                        "'YoungModulus' shall have a lower bound")
        bound = mgis_bv.getLowerBound(b, "YoungModulus")
        self.assertTrue(
            abs(bound - yg_min) < eps * yg_min,
            "invalid value for the 'YoungModulus' lower bound")
        bound = mgis_bv.getUpperBound(b, "YoungModulus")
        self.assertTrue(
            abs(bound - yg_max) < eps * yg_max,
            "invalid value for the 'YoungModulus' upper bound")
        self.assertTrue(mgis_bv.hasUpperBound(b, "YoungModulus"),
                        "'YoungModulus' shall have an upper bound")
        self.assertTrue(mgis_bv.hasPhysicalBounds(b, "YoungModulus"),
                        "'YoungModulus' shall have physical bounds")
        self.assertTrue(mgis_bv.hasLowerPhysicalBound(b, "YoungModulus"),
                        "'YoungModulus' shall have a lower physical bound")
        bound = mgis_bv.getLowerPhysicalBound(b, "YoungModulus")
        self.assertTrue(
            abs(bound) < eps * yg_max, "invalid value for the 'YoungModulus' "
            "physical lower bound")
        self.assertTrue(
            not mgis_bv.hasUpperPhysicalBound(b, "YoungModulus"),
            "'YoungModulus' shall not have an "
            "upper physical bound")
        self.assertTrue(not mgis_bv.hasBounds(b, "PoissonRatio"),
                        "'PoissonRatio' shall not have bounds")
        self.assertTrue(not mgis_bv.hasLowerBound(b, "PoissonRatio"),
                        "'PoissonRatio' shall not have a lower bound")
        self.assertTrue(not mgis_bv.hasUpperBound(b, "PoissonRatio"),
                        "'PoissonRatio' shall not have an upper bound")
        self.assertTrue(mgis_bv.hasPhysicalBounds(b, "PoissonRatio"),
                        "'PoissonRatio' shall have physical bounds")
        self.assertTrue(mgis_bv.hasLowerPhysicalBound(b, "PoissonRatio"),
                        "'PoissonRatio' shall have a lower physical bound")
        bound = mgis_bv.getLowerPhysicalBound(b, "PoissonRatio")
        self.assertTrue(
            abs(bound - nu_min) < eps,
            "invalid value for the 'PoissonRatio' physical "
            "lower bound")
        self.assertTrue(mgis_bv.hasUpperPhysicalBound(b, "PoissonRatio"),
                        "'PoissonRatio' shall have an upper physical bound")
        bound = mgis_bv.getUpperPhysicalBound(b, "PoissonRatio")
        self.assertTrue(
            abs(bound - nu_max) < eps * nu_max,
            "invalid value for the 'PoissonRatio' "
            "physical upper bound")
        # internal state variables
        self.assertTrue(mgis_bv.hasBounds(b, "StateVariable"),
                        "'StateVariable' shall have bounds")
        self.assertTrue(mgis_bv.hasLowerBound(b, "StateVariable"),
                        "'StateVariable' shall have a lower bound")
        bound = mgis_bv.getLowerBound(b, "StateVariable")
        self.assertTrue(
            abs(bound - iv_min) < eps * iv_min,
            "invalid value for the 'StateVariable' lower bound")
        bound = mgis_bv.getUpperBound(b, "StateVariable")
        self.assertTrue(
            abs(bound - iv_max) < eps * iv_max,
            "invalid value for the 'StateVariable' upper bound")
        self.assertTrue(mgis_bv.hasUpperBound(b, "StateVariable"),
                        "'StateVariable' shall have an upper bound")
        self.assertTrue(mgis_bv.hasPhysicalBounds(b, "StateVariable"),
                        "'StateVariable' shall have physical bounds")
        self.assertTrue(mgis_bv.hasLowerPhysicalBound(b, "StateVariable"),
                        "'StateVariable' shall have a lower physical bound")
        bound = mgis_bv.getLowerPhysicalBound(b, "StateVariable")
        self.assertTrue(
            abs(bound - iv_pmin) < eps,
            "invalid value for the 'StateVariable' "
            "physical lower bound")
        self.assertTrue(
            mgis_bv.hasUpperPhysicalBound(b, "StateVariable"),
            "'StateVariable' shall not have an "
            "upper physical bound")
        bound = mgis_bv.getUpperPhysicalBound(b, "StateVariable")
        self.assertTrue(
            abs(bound - iv_pmax) < eps * iv_pmax,
            "invalid value for the 'StateVariable' "
            "physical upper bound")
        # external state variable
        ev = "ExternalStateVariable"
        self.assertTrue(mgis_bv.hasBounds(b, ev),
                        "'ExternalStateVariable' shall have bounds")
        self.assertTrue(mgis_bv.hasLowerBound(b, ev),
                        "'ExternalStateVariable' shall have a lower bound")
        bound = mgis_bv.getLowerBound(b, ev)
        self.assertTrue(
            abs(bound - ev_min) < eps * ev_min,
            "invalid value for the 'ExternalStateVariable' "
            "lower bound")
        bound = mgis_bv.getUpperBound(b, ev)
        self.assertTrue(
            abs(bound - ev_max) < eps * ev_max,
            "invalid value for the 'ExternalStateVariable' "
            "upper bound")
        self.assertTrue(mgis_bv.hasUpperBound(b, ev),
                        "'ExternalStateVariable' shall have an upper bound")
        self.assertTrue(mgis_bv.hasPhysicalBounds(b, ev),
                        "'ExternalStateVariable' shall have physical bounds")
        self.assertTrue(
            mgis_bv.hasLowerPhysicalBound(b, ev),
            "'ExternalStateVariable' shall have a "
            "lower physical bound")
        bound = mgis_bv.getLowerPhysicalBound(b, ev)
        self.assertTrue(
            abs(bound - ev_pmin) < eps,
            "invalid value for the 'ExternalStateVariable' "
            "physical lower bound")
        self.assertTrue(
            mgis_bv.hasUpperPhysicalBound(b, ev),
            "'ExternalStateVariable' shall not have an "
            "upper physical bound")
        bound = mgis_bv.getUpperPhysicalBound(b, ev)
        self.assertTrue(
            abs(bound - ev_pmax) < eps * ev_pmax,
            "invalid value for the 'ExternalStateVariable' "
            "physical upper bound")


if __name__ == '__main__':
    unittest.main()
