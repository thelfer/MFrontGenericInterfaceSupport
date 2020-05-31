# -*- coding: utf-8 -*-

import os
import math
try:
    import unittest2 as unittest
except ImportError:
    import unittest
import mgis.behaviour as mgis_bv

class IntegrateTest2(unittest.TestCase):

    def test_pass(self):

        # path to the test library
        lib = os.environ['MGIS_TEST_BEHAVIOURS_LIBRARY']
        # reference values
        pref = [0,
                1.3523277308229e-11,
                1.0955374667213e-07,
                5.5890770166084e-06,
                3.2392193670428e-05,
                6.645865307584e-05,
                9.9676622883138e-05,
                0.00013302758358953,
                0.00016635821069889,
                0.00019969195920296,
                0.00023302522883648,
                0.00026635857194317,
                0.000299691903777,
                0.0003330252373404,
                0.00036635857063843,
                0.00039969190397718,
                0.00043302523730968,
                0.00046635857064314,
                0.00049969190397646,
                0.00053302523730979,
                0.00056635857064313]
        # modelling hypothesis
        h = mgis_bv.Hypothesis.Tridimensional
        # loading the behaviour        
        b = mgis_bv.load(lib,'Norton',h)
        # number of integration points
        nig = 100
        # material data manager
        m = mgis_bv.MaterialDataManager(b,nig)
        # index of the equivalent viscplastic strain in the array of
        # state variable
        o = mgis_bv.getVariableOffset(b.isvs,
                                      'EquivalentViscoplasticStrain',h)
        # strain increment per time step
        de = 5.e-5
        # time step increment
        dt = 180  
        # setting the temperature
        mgis_bv.setExternalStateVariable(m.s1,'Temperature', 293.15)
        # copy d.s1 in d.s0
        mgis_bv.update(m)
        # index of the first integration point
        ni = 0
        # index of the last integration point
        ne = nig-1
        # values of the equivalent plastic strain
        # for the first integration point
        pi = [m.s0.internal_state_variables[ni][o]]
        # values of the equivalent plastic strain
        # for the last integration point
        pe = [m.s0.internal_state_variables[ne][o]]
        # setting the gradient at the end of the first time step
        for i in range(0,nig):
            m.s1.gradients[i][0] = de
        ## integration
        for i in range(0,20):
            it = mgis_bv.IntegrationType.IntegrationWithoutTangentOperator
            mgis_bv.integrate(m, it, dt, 0, m.n)
            mgis_bv.update(m)
            for p in range(0,nig):
                m.s1.gradients[p][0] += de
            pi.append(m.s0.internal_state_variables[ni][o])
            pe.append(m.s0.internal_state_variables[ne][o])
        ## checks
        # comparison criterion
        eps = 1.e-12
        for i in range(0,21):
            self.assertTrue(abs(pi[i]-pref[i])<eps)
            self.assertTrue(abs(pe[i]-pref[i])<eps)
            
        pass
    
if __name__ == '__main__':
    unittest.main()
