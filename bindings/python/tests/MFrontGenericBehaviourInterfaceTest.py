# -*- coding: utf-8 -*-

import os
try:
    import unittest2 as unittest
except ImportError:
    import unittest
import mgis.behaviour as mgis_bv


class IntegrateTest(unittest.TestCase):

    def test_pass(self):
        btype = mgis_bv.BehaviourType
        bkinematic = mgis_bv.BehaviourKinematic
        # path to the test library
        lib = os.environ['MGIS_TEST_BEHAVIOURS_LIBRARY']
        version = os.environ['MGIS_TEST_TFEL_VERSION']
        h = mgis_bv.Hypothesis.Tridimensional
        b = mgis_bv.load(lib, 'Gurson', h)
        self.assertTrue(b.behaviour == "Gurson", "invalid behaviour name")
        self.assertTrue(b.hypothesis == h, "invalid hypothesis")
        self.assertTrue(b.source == "Gurson.mfront", "invalid source")
        self.assertTrue(b.tfel_version == version, "invalid TFEL version")
        self.assertTrue(len(b.mps) == 0,
                        "invalid number of material properties")
        self.assertTrue(b.getBehaviourType() == "StandardStrainBasedBehaviour",
                        "invalid behaviour type")
        self.assertTrue(b.btype == btype.StandardStrainBasedBehaviour,
                        "invalid behaviour type")
        self.assertTrue(b.btype == btype.STANDARDSTRAINBASEDBEHAVIOUR,
                        "invalid behaviour type")
        self.assertTrue(b.getKinematic() == "SmallStrainKinematic",
                        "invalid kinematic value")
        self.assertTrue(b.kinematic == bkinematic.SmallStrainKinematic,
                        "invalid kinematic value")
        self.assertTrue(b.kinematic == bkinematic.SMALLSTRAINKINEMATIC,
                        "invalid kinematic value")
        self.assertTrue(b.symmetry == mgis_bv.BehaviourSymmetry.Isotropic,
                        "invalid behaviour symmetry")
        self.assertTrue(len(b.gradients) == 1,
                        "invalid number of gradients")
        F = b.gradients[0]
        self.assertTrue(F.name == "Strain",
                        "invalid gradient name")
        self.assertTrue(F.type == mgis_bv.VariableType.Stensor,
                        "invalid gradient type")
        self.assertTrue(F.type == mgis_bv.VariableType.STENSOR,
                        "invalid gradient type")
        self.assertTrue(len(b.thermodynamic_forces) == 1,
                        "invalid number of thermodynamic_forces")
        pk1 = b.thermodynamic_forces[0]
        self.assertTrue(pk1.name == "Stress",
                        "invalid flux name")
        self.assertTrue(pk1.type == mgis_bv.VariableType.Stensor,
                        "invalid flux type")
        self.assertTrue(pk1.type == mgis_bv.VariableType.STENSOR,
                        "invalid flux type")
        self.assertTrue(len(b.isvs) == 4,
                        "invalid number of internal state variables")
        self.assertTrue(len(b.internal_state_variables) == 4,
                        "invalid number of internal state variables")
        self.assertTrue(b.isvs[0].name == "ElasticStrain",
                        "invalid name for the first internal state variable")
        self.assertTrue(b.isvs[0].type == mgis_bv.VariableType.STENSOR,
                        "invalid type for the first internal state variable")
        self.assertTrue(b.isvs[1].name == "EquivalentPlasticStrain",
                        "invalid name for the second internal state variable")
        self.assertTrue(b.isvs[1].type == mgis_bv.VariableType.SCALAR,
                        "invalid type for the second internal state variable")
        self.assertTrue(b.isvs[1].type == mgis_bv.VariableType.Scalar,
                        "invalid type for the second internal state variable")
        self.assertTrue(b.isvs[1].getType() == 'Scalar',
                        "invalid type for the second internal state variable")
        self.assertTrue(b.isvs[2].name == "MatrixEquivalentPlasticStrain",
                        "invalid name for the third internal state variable")
        self.assertTrue(b.isvs[2].type == mgis_bv.VariableType.SCALAR,
                        "invalid type for the third internal state variable")
        self.assertTrue(b.isvs[3].name == "Porosity",
                        "invalid name for the third internal state variable")
        self.assertTrue(b.isvs[3].type == mgis_bv.VariableType.SCALAR,
                        "invalid type for the fourth internal state variable")
        # external state variables
        self.assertTrue(len(b.esvs) == 1,
                        "invalid number of external state variables")
        self.assertTrue(len(b.external_state_variables) == 1,
                        "invalid number of external state variables")
        self.assertTrue(b.esvs[0].name == "Temperature",
                        "invalid name for the first external state variable")
        self.assertTrue(b.esvs[0].type == mgis_bv.VariableType.SCALAR,
                        "invalid type for the first external state variable")
        # parameters
        params = b.parameters
        self.assertTrue(len(params) == 15,
                        "invalid number of parameters")
        self.assertTrue('epsilon' in params,
                        "'epsilon is not a parameter")
        self.assertTrue(len(b.iparams) == 0,
                        "invalid number of integer parameters")
        self.assertTrue(len(b.integer_parameters) == 0,
                        "invalid number of integer parameters")
        self.assertTrue(len(b.unsigned_short_parameters) == 1,
                        "invalid number of unsigned short parameters")
        self.assertTrue(b.usparams[0] == 'iterMax',
                        "invalid name for the first unsigned short parameter")


if __name__ == '__main__':
    unittest.main()
