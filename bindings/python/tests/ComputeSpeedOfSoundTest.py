# -*- coding: utf-8 -*-

import os
import math
try:
    import unittest2 as unittest
except ImportError:
    import unittest
import mgis.behaviour as mgis_bv


class ComputeSpeedOfSoundTest(unittest.TestCase):
    def test_pass(self):

        # path to the test library
        lib = os.environ['MGIS_TEST_BEHAVIOURS_LIBRARY']
        # Young modulus and mass density
        yg = 150e9
        rho = 7850
        # modelling hypothesis
        h = mgis_bv.Hypothesis.Tridimensional
        # loading the behaviour
        b = mgis_bv.load(lib, 'ComputeSpeedOfSound', h)
        # material data manager
        m = mgis_bv.MaterialDataManager(b, 100)
        # initialize the external state variable and the mass density
        mgis_bv.setExternalStateVariable(m.s1, 'Temperature', 293.15)
        mgis_bv.setMassDensity(m.s1, rho)
        # copy d.s1 in d.s0
        mgis_bv.update(m)
        m.s1.gradients[:, 0] = 5.e-5
        # integration with the computation of the speed of sound
        opts = mgis_bv.BehaviourIntegrationOptions()
        opts.integration_type = \
            mgis_bv.IntegrationType.IntegrationWithoutTangentOperator
        opts.compute_speed_of_sound = True
        r = mgis_bv.integrate(m, opts, 180, 0, m.n)
        self.assertTrue(r.exit_status == 1, 'integration failed')
        v_ref = math.sqrt(yg / rho)
        self.assertTrue(
            abs(m.speed_of_sound[0] - v_ref) < 1.e-12 * v_ref,
            'invalid value for the speed of sound')

        pass

    def test_behaviour_data(self):
        lib = os.environ['MGIS_TEST_BEHAVIOURS_LIBRARY']
        yg = 150e9
        rho = 7850
        h = mgis_bv.Hypothesis.Tridimensional
        b = mgis_bv.load(lib, 'ComputeSpeedOfSound', h)
        d = mgis_bv.BehaviourData(b)
        d.s1.mass_density = rho
        d.s0.external_state_variables[:] = 293.15
        d.s1.external_state_variables[:] = 293.15
        d.s1.gradients[0] = 5.e-5
        d.dt = 180
        # integration without tangent operator, with the speed of sound
        d.K[0] = 100
        v = mgis_bv.make_view(d)
        r = mgis_bv.integrate(v, b)
        self.assertTrue(r == 1, 'integration failed')
        v_ref = math.sqrt(yg / rho)
        self.assertTrue(
            abs(d.speed_of_sound - v_ref) < 1.e-12 * v_ref,
            'invalid value for the speed of sound')


if __name__ == '__main__':
    unittest.main()
