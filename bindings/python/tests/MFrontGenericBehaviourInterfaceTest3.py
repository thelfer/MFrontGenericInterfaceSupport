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
        version = os.environ['MGIS_TEST_TFEL_VERSION']
        h = mgis_bv.Hypothesis.Tridimensional
        o = mgis_bv.FiniteStrainBehaviourOptions()
        o.stress_measure = mgis_bv.FiniteStrainBehaviourOptionsStressMeasure.PK1
        o.tangent_operator = mgis_bv.FiniteStrainBehaviourOptionsTangentOperator.DPK1_DF
        b = mgis_bv.load(o, lib,'FiniteStrainSingleCrystal',h)
        self.assertTrue(len(b.gradients) == 1, 'invalid number of gradients')
        self.assertTrue(b.gradients[0].name == 'DeformationGradient',
                        'invalid gradient name')
        self.assertTrue(b.gradients[0].getType() == 'Tensor',
                        'invalid gradient type')
        self.assertTrue(b.gradients[0].type == mgis_bv.VariableType.TENSOR,
                        'invalid gradient type')
        self.assertTrue(b.gradients[0].type == mgis_bv.VariableType.Tensor,
                        'invalid gradient type')
        self.assertTrue(len(b.thermodynamic_forces) == 1,
                        'invalid number of thermodynamic_forces')
        self.assertTrue(b.thermodynamic_forces[0].name == 'FirstPiolaKirchhoffStress',
                        'invalid thermodynamic force  name')
        self.assertTrue(b.thermodynamic_forces[0].getType() == 'Tensor',
                        'invalid thermodynamic force  type')
        self.assertTrue(b.thermodynamic_forces[0].type == mgis_bv.VariableType.TENSOR,
                        'invalid thermodynamic force  type')
        
if __name__ == '__main__':
    unittest.main()
