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
        self.assertTrue(b.behaviour == 'FiniteStrainSingleCrystal',
                        'invalid behaviour name')
        self.assertTrue(b.hypothesis == h, 'invalid hypothesis')
        self.assertTrue(b.source == 'FiniteStrainSingleCrystal.mfront',
                        'invalid source')
        self.assertTrue(b.tfel_version == version, 'invalid TFEL version')
        self.assertTrue(b.getBehaviourType() == 'StandardFiniteStrainBehaviour',
                        'invalid behaviour type')
        self.assertTrue(b.btype == mgis_bv.BehaviourType.STANDARDFINITESTRAINBEHAVIOUR,
                        'invalid behaviour type')
        self.assertTrue(b.btype == mgis_bv.BehaviourType.StandardFiniteStrainBehaviour,
                        'invalid behaviour type')
        self.assertTrue(b.getKinematic() == 'F_CAUCHY',
                        'invalid kinematic value')
        self.assertTrue(b.kinematic == mgis_bv.BehaviourKinematic.FINITESTRAINKINEMATIC_F_CAUCHY,
                        'invalid kinematic value')
        self.assertTrue(b.kinematic == mgis_bv.BehaviourKinematic.FiniteStrainKinematic_F_Cauchy,
                        'invalid kinematic value')
        self.assertTrue(b.getSymmetry() == 'Orthotropic',
                        'invalid behaviour symmetry')
        self.assertTrue(b.symmetry == mgis_bv.BehaviourSymmetry.ORTHOTROPIC,
                        'invalid behaviour symmetry')
        self.assertTrue(b.symmetry == mgis_bv.BehaviourSymmetry.Orthotropic,
                        'invalid behaviour symmetry')
        self.assertTrue(len(b.gradients) == 1, 'invalid number of gradients')
        F = b.gradients[0]
        self.assertTrue(F.name == 'DeformationGradient',
                        'invalid gradient name')
        self.assertTrue(F.getType() == 'Tensor',
                        'invalid gradient type')
        self.assertTrue(F.type == mgis_bv.VariableType.TENSOR,
                        'invalid gradient type')
        self.assertTrue(F.type == mgis_bv.VariableType.Tensor,
                        'invalid gradient type')
        self.assertTrue(len(b.thermodynamic_forces) == 1,
                        'invalid number of thermodynamic_forces')
        pk1 = b.thermodynamic_forces[0]
        self.assertTrue(pk1.name == 'FirstPiolaKirchhoffStress',
                        'invalid flux name')
        self.assertTrue(pk1.getType() == 'Tensor',
                        'invalid flux type')
        self.assertTrue(pk1.type == mgis_bv.VariableType.TENSOR,
                        'invalid flux type')
        self.assertTrue(pk1.type == mgis_bv.VariableType.Tensor,
                        'invalid flux type')
        self.assertTrue(len(b.mps) == 16,
                        'invalid number of material properties')

        def check_mp_name(mp, n):
            self.assertTrue(mp.name ==n,
                            "invalid material property name, " +
                            "expected '" + n + "'")

        for mp in b.mps:
            self.assertTrue(mp.getType() == 'Scalar',
                            "invalid material property type '" + mp.name + "'")
            self.assertTrue(mp.type == mgis_bv.VariableType.SCALAR,
                            "invalid material property type '" + mp.name + "'")
            self.assertTrue(mp.type == mgis_bv.VariableType.Scalar,
                            "invalid material property type '" + mp.name + "'")
        check_mp_name(b.mps[0], 'YoungModulus1')
        check_mp_name(b.mps[1], 'YoungModulus2')
        check_mp_name(b.mps[2], 'YoungModulus3')
        check_mp_name(b.mps[3], 'PoissonRatio12')
        check_mp_name(b.mps[4], 'PoissonRatio23')
        check_mp_name(b.mps[5], 'PoissonRatio13')
        check_mp_name(b.mps[6], 'ShearModulus12')
        check_mp_name(b.mps[7], 'ShearModulus23')
        check_mp_name(b.mps[8], 'ShearModulus13')
        check_mp_name(b.mps[9], 'm')
        check_mp_name(b.mps[10], 'K')
        check_mp_name(b.mps[11], 'C')
        check_mp_name(b.mps[12], 'R0')
        check_mp_name(b.mps[13], 'Q')
        check_mp_name(b.mps[14], 'b')
        check_mp_name(b.mps[15], 'd1')
        self.assertTrue(len(b.isvs) == 37,
                        'invalid number of internal state variables')
        
        def check_iv_name(iv, n, i):
            vn = n + '[' + str(i) + ']'
            self.assertTrue(iv.name == vn,
                            "invalid internal state variable name, " +
                            "expected '" + vn + "'")
            self.assertTrue(iv.getType() == 'Scalar',
                            'invalid type for internal ' +
                            'state variable \'' + vn + '\'')
            self.assertTrue(iv.type == mgis_bv.VariableType.SCALAR,
                            'invalid type for internal ' +
                            'state variable \'' + vn + '\'')
            self.assertTrue(iv.type == mgis_bv.VariableType.Scalar,
                            'invalid type for internal ' +
                            'state variable \'' + vn + '\'')
        for i in range(0, 12):
            check_iv_name(b.isvs[i], 'g', i)
            check_iv_name(b.isvs[13 + i], 'p', i)
            check_iv_name(b.isvs[25 + i], 'a', i)
        self.assertTrue(b.isvs[12].name == 'Fe',
                        "invalid name for internal state variable 'Fe'")
        self.assertTrue(b.isvs[12].getType() == "Tensor",
                        "invalid type for internal state variable 'Fe'")
        self.assertTrue(b.isvs[12].type == mgis_bv.VariableType.TENSOR,
                        "invalid type for internal state variable 'Fe'")
        self.assertTrue(b.isvs[12].type == mgis_bv.VariableType.Tensor,
                        "invalid type for internal state variable 'Fe'")
        self.assertTrue(len(b.esvs) == 1,
                        'invalid number of external state variables')
        self.assertTrue(b.esvs[0].name == 'Temperature',
                        'invalid name for the first external state variable')
        self.assertTrue(b.esvs[0].getType() == 'Scalar',
                        'invalid type for the first external state variable')
        self.assertTrue(b.esvs[0].type == mgis_bv.VariableType.SCALAR,
                        'invalid type for the first external state variable')
        self.assertTrue(b.esvs[0].type == mgis_bv.VariableType.Scalar,
                        'invalid type for the first external state variable')


if __name__ == '__main__':
    unittest.main()
