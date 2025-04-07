# -*- coding: utf-8 -*-

import os
import math
import numpy
try:
    import unittest2 as unittest
except ImportError:
    import unittest
import mgis.behaviour as mgis_bv
import mgis.model


class IntegrateTest5(unittest.TestCase):
    def test_pass(self):

        # path to the test library
        lib = os.environ['MGIS_TEST_MODELS_LIBRARY']
        # modelling hypothesis
        h = mgis_bv.Hypothesis.Tridimensional
        # loading the behaviour
        model = mgis.model.load(lib, 'ode_rk54', h)
        # default value of parameter A
        A = model.getParameterDefaultValue('A')
        # number of integration points
        nig = 100
        # material data manager
        m = mgis_bv.MaterialDataManager(model, nig)
        # index of x in the array of state variable
        o = mgis_bv.getVariableOffset(model.isvs, 'x', h)
        # time step increment
        dt = 0.1
        # setting the temperature
        T = 293.15 * numpy.ones(nig)
        # type of storage
        Ts = mgis_bv.MaterialStateManagerStorageMode.ExternalStorage
        mgis_bv.setExternalStateVariable(m.s1, 'Temperature', T, Ts)
        # Initial value of x
        for n in range(0, nig):
            m.s1.internal_state_variables[n][o] = 1
        # copy d.s1 in d.s0
        mgis_bv.update(m)
        # index of the first integration point
        ni = 0
        # index of the last integration point
        ne = nig - 1
        # values of x for the first integration point
        xi = [m.s0.internal_state_variables[ni][o]]
        # values of x for the last integration point
        xe = [m.s0.internal_state_variables[ne][o]]
        # integration
        for i in range(0, 10):
            it = mgis_bv.IntegrationType.IntegrationWithoutTangentOperator
            mgis_bv.integrate(m, it, dt, 0, m.n)
            mgis_bv.update(m)
            xi.append(m.s1.internal_state_variables[ni][o])
            xe.append(m.s1.internal_state_variables[ne][o])
        # checks
        # comparison criterion
        eps = 1.e-10
        t = 0
        for i in range(0, 11):
            x_ref = math.exp(-A * t)
            self.assertTrue(abs(xi[i] - x_ref) < eps)
            self.assertTrue(abs(xe[i] - x_ref) < eps)
            t = t + dt

        pass


if __name__ == '__main__':
    unittest.main()
