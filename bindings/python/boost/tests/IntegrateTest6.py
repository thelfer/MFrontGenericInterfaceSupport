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


class IntegrateTest6(unittest.TestCase):
    def test_pass(self):

        # path to the test library
        lib = os.environ['MGIS_TEST_MODELS_LIBRARY']
        # modelling hypothesis
        h = mgis_bv.Hypothesis.Tridimensional
        # loading the behaviour
        model = mgis.model.load(lib, 'ode_rk54', h)
        # default value of parameter A
        A = model.getParameterDefaultValue('A')
        # material data manager
        d = mgis_bv.BehaviourData(model)
        # index of x in the array of state variable
        o = mgis_bv.getVariableOffset(model.isvs, 'x', h)
        # time step increment
        d.dt = 0.1
        # type of storage
        mgis_bv.setExternalStateVariable(d.s1, 'Temperature', 293.15)
        # Initial value of x
        d.s1.internal_state_variables[o] = 1
        # copy d.s1 in d.s0
        mgis_bv.update(d)
        # values of  x
        xvalues = [d.s0.internal_state_variables[o]]
        # integration
        for i in range(0, 10):
            mgis_bv.integrate(d, model)
            mgis_bv.update(d)
            xvalues.append(d.s1.internal_state_variables[o])
        # checks
        # comparison criterion
        eps = 1.e-10
        t = 0
        for i in range(0, 11):
            x_ref = math.exp(-A * t)
            self.assertTrue(abs(xvalues[i] - x_ref) < eps)
            t = t + d.dt

        pass


if __name__ == '__main__':
    unittest.main()
