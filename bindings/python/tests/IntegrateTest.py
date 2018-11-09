# -*- coding: utf-8 -*-

import os
import math
try:
    import unittest2 as unittest
except ImportError:
    import unittest
import mgis.behaviour as mgis_bv

class IntegrateTest(unittest.TestCase):

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
        # comparison criterion
        eps = 1.e-12
        
        b = mgis_bv.load(lib,'Norton',
                         mgis_bv.Hypothesis.Tridimensional)
        d = mgis_bv.BehaviourData(b)
        o = mgis_bv.getVariableOffset(b.isvs, 'EquivalentViscoplasticStrain', b.hypothesis)
        
        # strain increment per time step
        de = 5.e-5
        # time step
        d.dt = 180  
	
        # setting the temperature
        mgis_bv.setExternalStateVariable(d.s1, 'Temperature', 293.15)
        
        # copy d.s1 in d.s0
        mgis_bv.update(d)
        d.s1.gradients[0] = de
        
        # equivalent plastic strain
        p = [d.s0.internal_state_variables[o]]
        
        # integrate the behaviour
        for i in range(0,20):
            mgis_bv.integrate(d, b)
            mgis_bv.update(d)
            d.s1.gradients[0] += de
            p.append(d.s1.internal_state_variables[o])
        
        # check-in results
        for i in range(0,20):
            self.assertTrue(abs(p[i]-pref[i])<eps)
        
        pass

def main(args):
    pass
    
if __name__ == '__main__':
    unittest.main()
