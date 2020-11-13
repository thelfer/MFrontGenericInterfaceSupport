# -*- coding: utf-8 -*-

import os
try:
    import unittest2 as unittest
except ImportError:
    import unittest
import mgis.behaviour as mgis_bv


class IntegrateTest(unittest.TestCase):

    def test_pass(self):

        stress_options = mgis_bv.FiniteStrainBehaviourOptionsStressMeasure
        to_options = mgis_bv.FiniteStrainBehaviourOptionsTangentOperator
        # path to the test library
        lib = os.environ['MGIS_TEST_BEHAVIOURS_LIBRARY']
        version = os.environ['MGIS_TEST_TFEL_VERSION']
        h = mgis_bv.Hypothesis.Tridimensional
        o = mgis_bv.FiniteStrainBehaviourOptions()
        o.stress_measure = stress_options.PK1
        o.tangent_operator = to_options.DPK1_DF
        b = mgis_bv.load(o, lib, 'FiniteStrainSingleCrystal', h)
        self.assertTrue(b.tfel_version == version, "invalid TFEL version")
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
        pk1 = b.thermodynamic_forces[0]
        self.assertTrue(pk1.name == 'FirstPiolaKirchhoffStress',
                        'invalid thermodynamic force  name')
        self.assertTrue(pk1.getType() == 'Tensor',
                        'invalid thermodynamic force  type')
        self.assertTrue(pk1.type == mgis_bv.VariableType.TENSOR,
                        'invalid thermodynamic force  type')


if __name__ == '__main__':
    unittest.main()
