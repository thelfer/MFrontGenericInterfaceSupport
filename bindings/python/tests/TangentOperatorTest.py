# -*- coding: utf-8 -*-

# tests of the tangent and prediction operators returned by the integrations
# of a material data manager and of a behaviour data, with the `TangentKind`
# behaviour whose operator identifies the type of integration

import os
import numpy
try:
    import unittest2 as unittest
except ImportError:
    import unittest
import mgis
import mgis.behaviour as mgis_bv

types = [
    mgis_bv.IntegrationType.PredictionWithTangentOperator,
    mgis_bv.IntegrationType.PredictionWithSecantOperator,
    mgis_bv.IntegrationType.PredictionWithElasticOperator,
    mgis_bv.IntegrationType.IntegrationWithElasticOperator,
    mgis_bv.IntegrationType.IntegrationWithSecantOperator,
    mgis_bv.IntegrationType.IntegrationWithTangentOperator,
    mgis_bv.IntegrationType.IntegrationWithConsistentTangentOperator
]


def expected_operator(t, k):
    """operator computed by the `TangentKind` behaviour for the type of
    integration t and the value k of the material property"""
    K = int(t) * k * numpy.eye(6)
    if int(t) > 0:
        K[0, 1] += k / 2
    return K


class TangentOperatorTest(unittest.TestCase):

    def test_material_data_manager(self):
        lib = os.environ['MGIS_TEST_BEHAVIOURS_LIBRARY']
        h = mgis_bv.Hypothesis.Tridimensional
        b = mgis_bv.load(lib, 'TangentKind', h)
        n = 4
        k = numpy.array([10., 20., 30., 40.])
        p = mgis.ThreadPool(2)
        for t in types:
            # mode 0: range [1, 3), mode 1: thread pool on all the points
            for mode in range(2):
                m = mgis_bv.MaterialDataManager(b, n)
                for s in [m.s0, m.s1]:
                    mgis_bv.setMaterialProperty(
                        s, 'k', k,
                        mgis_bv.MaterialStateManagerStorageMode.LocalStorage)
                    mgis_bv.setExternalStateVariable(s, 'Temperature', 293.15)
                for i in range(n):
                    m.s1.gradients[i, 0] = 1.e-3 * (i + 1)
                if mode == 0:
                    r = mgis_bv.integrate(m, t, 1., 1, 3)
                else:
                    r = mgis_bv.integrate(p, m, t, 1.)
                self.assertTrue(r == 1, 'integration failed')
                self.assertTrue(m.K.shape == (n, 6, 6),
                                'invalid shape of the tangent operators')
                for i in range(n):
                    # the points outside the range are left untouched
                    treated = (mode == 1) or (1 <= i < 3)
                    K = expected_operator(t, k[i]) if treated \
                        else numpy.zeros((6, 6))
                    self.assertTrue(
                        numpy.allclose(m.K[i], K, rtol=0, atol=1.e-12 * k[i]),
                        'invalid value of the operator')
                    # a prediction leaves the internal state variables
                    # untouched
                    a = 1 if (treated and (int(t) > 0)) else 0
                    self.assertTrue(
                        abs(m.s1.internal_state_variables[i, 0] - a) < 1.e-14,
                        'invalid value of the internal state variable')

    def test_behaviour_data(self):
        lib = os.environ['MGIS_TEST_BEHAVIOURS_LIBRARY']
        h = mgis_bv.Hypothesis.Tridimensional
        b = mgis_bv.load(lib, 'TangentKind', h)
        for t in types:
            d = mgis_bv.BehaviourData(b)
            for s in [d.s0, d.s1]:
                mgis_bv.setMaterialProperty(s, 'k', 7.)
                mgis_bv.setExternalStateVariable(s, 'Temperature', 293.15)
            d.s1.gradients[0] = 1.e-3
            d.dt = 1
            d.K[0, 0] = int(t)
            v = mgis_bv.make_view(d)
            self.assertTrue(mgis_bv.integrate(v, b) == 1, 'integration failed')
            self.assertTrue(
                numpy.allclose(d.K, expected_operator(t, 7.), rtol=0,
                               atol=1.e-12 * 7),
                'invalid value of the operator of the behaviour data')


if __name__ == '__main__':
    unittest.main()
