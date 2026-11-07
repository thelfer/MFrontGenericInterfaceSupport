# -*- coding: utf-8 -*-

import os

try:
    import unittest2 as unittest
except ImportError:
    import unittest
import mgis.behaviour as mgis_bv


class ArgumentsCheckTest(unittest.TestCase):
    def __load(self, name, h=mgis_bv.Hypothesis.Tridimensional):
        lib = os.environ['MGIS_TEST_BEHAVIOURS_LIBRARY']
        return mgis_bv.load(lib, name, h)

    def test_behaviour_data_layout(self):
        d = mgis_bv.BehaviourData(self.__load('Elasticity'))
        gurson = self.__load('Gurson')
        self.assertRaises(RuntimeError, mgis_bv.integrate, d, gurson)
        self.assertRaises(RuntimeError, mgis_bv.integrate_debug, d, gurson)
        # same behaviour, another modelling hypothesis
        d = mgis_bv.BehaviourData(
            self.__load('Norton', mgis_bv.Hypothesis.PlaneStrain))
        self.assertRaises(RuntimeError, mgis_bv.integrate, d,
                          self.__load('Norton'))
        # the same behaviour loaded twice
        d = mgis_bv.BehaviourData(self.__load('Norton'))
        mgis_bv.setExternalStateVariable(d.s0, 'Temperature', 293.15)
        mgis_bv.setExternalStateVariable(d.s1, 'Temperature', 293.15)
        d.dt = 1
        self.assertEqual(mgis_bv.integrate(d, self.__load('Norton')), 1)


if __name__ == '__main__':
    unittest.main()
