# -*- coding: utf-8 -*-

import os
import array
import decimal
import fractions
import numpy

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

    def test_external_array_not_reallocated(self):
        # the object exporting the buffer can't be resized while it is used
        b = self.__load('Norton')
        n = 4
        extra = [0.] * 1000
        m = mgis_bv.MaterialDataManager(b, n)
        T = array.array('d', [293.15] * n)
        mgis_bv.setExternalStateVariable(
            m.s1, 'Temperature', T,
            mgis_bv.MaterialStateManagerStorageMode.ExternalStorage)
        self.assertRaises(BufferError, T.extend, extra)
        # the buffer is released once the manager no longer uses it
        mgis_bv.setExternalStateVariable(m.s1, 'Temperature', 293.15)
        T.extend(extra)
        vs = array.array('d', [0.] * n)
        m.useExternalArrayOfSpeedOfSounds(vs)
        self.assertRaises(BufferError, vs.extend, extra)
        i = mgis_bv.MaterialDataManagerInitializer()
        K = array.array('d', [0.] * (n * 36))
        i.bindTangentOperator(K)
        m = mgis_bv.MaterialDataManager(b, n, i)
        self.assertRaises(BufferError, K.extend, extra)

    def test_scalars_of_a_state(self):
        d = mgis_bv.BehaviourData(self.__load('Norton'))
        for v in (numpy.float32(293.), numpy.int32(293), numpy.array(293.),
                  decimal.Decimal(293), fractions.Fraction(293)):
            for k in ('Temperature', 0):
                d.s1.external_state_variables[0] = 0
                mgis_bv.setExternalStateVariable(d.s1, k, v)
                self.assertEqual(d.s1.external_state_variables[0], 293)


if __name__ == '__main__':
    unittest.main()
