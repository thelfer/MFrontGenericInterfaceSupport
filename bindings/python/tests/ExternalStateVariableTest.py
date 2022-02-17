# -*- coding: utf-8 -*-

import os
import numpy

try:
    import unittest2 as unittest
except ImportError:
    import unittest
import mgis.behaviour as mgis_bv

class ExternalStateVariableTest(unittest.TestCase):

    def __get_behaviour(self, h):
        lib = os.environ['MGIS_TEST_BEHAVIOURS_LIBRARY']
        return mgis_bv.load(lib, 'TensorialExternalStateVariableTest', h)
    
    def test_behaviour(self):
        h = mgis_bv.Hypothesis.Tridimensional
        b = self.__get_behaviour(h)
        # path to the test library
        version = os.environ['MGIS_TEST_TFEL_VERSION']
        self.assertTrue(b.behaviour == "TensorialExternalStateVariableTest",
                        "invalid behaviour name")
        self.assertTrue(b.hypothesis == h, "invalid hypothesis")
        self.assertTrue(b.source == "TensorialExternalStateVariableTest.mfront",
                        "invalid source")
        self.assertTrue(b.tfel_version == version, "invalid TFEL version")
        self.assertTrue(mgis_bv.getArraySize(b.esvs, b.hypothesis) == 55,
                        "invalid array size for the external state variables")
        esvnames = ["Temperature", "v_esv",     "v2_esv[0]", "v2_esv[1]", "s_esv",
                    "s2_esv[0]",   "s2_esv[1]", "t_esv",     "t2_esv[0]", "t2_esv[1]"]
        esvtypes = [mgis_bv.VariableType.SCALAR, mgis_bv.VariableType.VECTOR,
                    mgis_bv.VariableType.VECTOR, mgis_bv.VariableType.VECTOR,
                    mgis_bv.VariableType.STENSOR, mgis_bv.VariableType.STENSOR,
                    mgis_bv.VariableType.STENSOR, mgis_bv.VariableType.TENSOR,
                    mgis_bv.VariableType.TENSOR,  mgis_bv.VariableType.TENSOR]
        esvsizes = [1, 3, 3, 3, 6, 6, 6, 9, 9, 9]
        for i in range(0, 10):
            self.assertTrue(b.esvs[i].name == esvnames[i],
                            "invalid external state variable name")
            self.assertTrue(b.esvs[i].type == esvtypes[i],
                            "invalid external variable type")
            self.assertTrue(mgis_bv.getVariableSize(b.esvs[i], b.hypothesis) == esvsizes[i],
                            "invalid external variable type")

        pass

    def test_behaviour_data(self):
        eps = 1.e-14
        h = mgis_bv.Hypothesis.Tridimensional
        b = self.__get_behaviour(h)
        d = mgis_bv.BehaviourData(b)
        self.assertTrue(len(d.s0.external_state_variables) == 55,
                        "invalid array size for the external state variables")
        self.assertTrue(len(d.s1.external_state_variables) == 55,
                        "invalid array size for the external state variables")
        v_esv_values = numpy.asarray([1, 2, 3], dtype=numpy.float64)
        print(v_esv_values)
        mgis_bv.setExternalStateVariable(d.s1, "v_esv", v_esv_values)
        mgis_bv.integrate(d, b)
        for i in range(0, 3):
            self.assertTrue(abs(d.s1.external_state_variables[i + 1] -
                                v_esv_values[i]) < eps,
                            "invalid external state variable value")
            self.assertTrue(abs(d.s1.internal_state_variables[i] -
                                v_esv_values[i]) <  eps,
                            "invalid internal state variable value")
        pass

    def test_material_data_manager(self):
        eps = 1.e-14
        dt = 0.1
        h = mgis_bv.Hypothesis.Tridimensional
        b = self.__get_behaviour(h)
        d = mgis_bv.MaterialDataManager(b, 2)
        mgis_bv.setMaterialProperty(d.s0, "YoungModulus", 150e9)
        mgis_bv.setMaterialProperty(d.s0, "PoissonRatio", 0.3)
        mgis_bv.setMaterialProperty(d.s1, "YoungModulus", 150e9)
        mgis_bv.setMaterialProperty(d.s1, "PoissonRatio", 0.3)
        # 
        zeros = numpy.zeros(9)
        # v_esv is uniform
        v_esv_values = numpy.asarray([1, 2, 3], dtype=numpy.float64)
        mgis_bv.setExternalStateVariable(d.s0, "v_esv", zeros[0:3],
                                         mgis_bv.MaterialStateManagerStorageMode.EXTERNAL_STORAGE)
        mgis_bv.setExternalStateVariable(d.s1, "v_esv", v_esv_values,
                                         mgis_bv.MaterialStateManagerStorageMode.LOCAL_STORAGE)
        # v2_esv[0] is not uniform
        v2_esv0_values = numpy.asarray([1, 2, 3, 7, 6, 5], dtype=numpy.float64)
        mgis_bv.setExternalStateVariable(d.s0, "v2_esv[0]",zeros[0:3],
                                         mgis_bv.MaterialStateManagerStorageMode.EXTERNAL_STORAGE)
        mgis_bv.setExternalStateVariable(d.s1, "v2_esv[0]", v2_esv0_values,
                                         mgis_bv.MaterialStateManagerStorageMode.EXTERNAL_STORAGE)
        v2_esv1_values =  numpy.asarray([9, 8, 2, 3, 6, 4], dtype=numpy.float64)
        mgis_bv.setExternalStateVariable(d.s0, "v2_esv[1]", zeros[0:3],
                                         mgis_bv.MaterialStateManagerStorageMode.EXTERNAL_STORAGE)
        mgis_bv.setExternalStateVariable(d.s1, "v2_esv[1]", v2_esv1_values,
                                         mgis_bv.MaterialStateManagerStorageMode.EXTERNAL_STORAGE)
        for s in [d.s0, d.s1]:
            mgis_bv.setExternalStateVariable(s, "Temperature", zeros[0:1],
                                             mgis_bv.MaterialStateManagerStorageMode.EXTERNAL_STORAGE)
            mgis_bv.setExternalStateVariable(s, "s_esv", zeros[0:6],
                                             mgis_bv.MaterialStateManagerStorageMode.EXTERNAL_STORAGE)
            mgis_bv.setExternalStateVariable(s, "s2_esv[0]", zeros[0:6],
                                             mgis_bv.MaterialStateManagerStorageMode.EXTERNAL_STORAGE)
            mgis_bv.setExternalStateVariable(s, "s2_esv[1]", zeros[0:6],
                                             mgis_bv.MaterialStateManagerStorageMode.EXTERNAL_STORAGE)
            mgis_bv.setExternalStateVariable(s, "t_esv", zeros[0:9],
                                             mgis_bv.MaterialStateManagerStorageMode.EXTERNAL_STORAGE)
            mgis_bv.setExternalStateVariable(s, "t2_esv[0]", zeros[0:9],
                                             mgis_bv.MaterialStateManagerStorageMode.EXTERNAL_STORAGE)
            mgis_bv.setExternalStateVariable(s, "t2_esv[1]", zeros[0:9],
                                             mgis_bv.MaterialStateManagerStorageMode.EXTERNAL_STORAGE)
        # checking if the external storage do work as expected
        v2_esv1_values[3] = -1
        mgis_bv.integrate(d, mgis_bv.IntegrationType.INTEGRATION_NO_TANGENT_OPERATOR, dt, 0, 2)
        for i in range (0, 3):
            self.assertTrue(abs(d.s1.internal_state_variables[0, i] - v_esv_values[i]) < eps,
                            "invalid internal state variable value")
            self.assertTrue(abs(d.s1.internal_state_variables[1, i] - v_esv_values[i]) < eps,
                            "invalid internal state variable value")
            self.assertTrue(abs(d.s1.internal_state_variables[0, i + 3] -
                                v2_esv0_values[i]) < eps,
                            "invalid internal state variable value")
            self.assertTrue(abs(d.s1.internal_state_variables[1, i + 3] -
                                v2_esv0_values[i + 3]) < eps,
                            "invalid internal state variable value")
            self.assertTrue(abs(d.s1.internal_state_variables[0, i + 6] -
                                v2_esv1_values[i]) < eps,
                             "invalid internal state variable value")
            self.assertTrue(abs(d.s1.internal_state_variables[1, i + 6] -
                                v2_esv1_values[i + 3]) < eps,
                            "invalid internal state variable value")

if __name__ == '__main__':
    unittest.main()
