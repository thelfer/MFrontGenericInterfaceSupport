# -*- coding: utf-8 -*-

import os
import numpy

try:
    import unittest2 as unittest
except ImportError:
    import unittest
import mgis.behaviour as mgis_bv


class PostProcessingTest(unittest.TestCase):
    def __get_behaviour(self, h):
        lib = os.environ['MGIS_TEST_BEHAVIOURS_LIBRARY']
        return mgis_bv.load(lib, 'PostProcessingTest', h)

    def test_behaviour(self):
        h = mgis_bv.Hypothesis.Tridimensional
        b = self.__get_behaviour(h)
        # path to the test library
        version = os.environ['MGIS_TEST_TFEL_VERSION']
        self.assertTrue(b.behaviour == "PostProcessingTest",
                        "invalid behaviour name")
        self.assertTrue(b.hypothesis == h, "invalid hypothesis")
        self.assertTrue(b.source == "PostProcessingTest.mfront",
                        "invalid source")
        self.assertTrue(b.tfel_version == version, "invalid TFEL version")
        postprocessings = b.getPostProcessingsNames()
        self.assertTrue(
            len(postprocessings) == 1, "invalid number of post-processings")
        self.assertTrue("PrincipalStrain" in postprocessings,
                        "invalid post-processing names")
        o = b.getPostProcessingOutputs("PrincipalStrain")[0]
        self.assertTrue(o.getType() == "Vector3D", "invalid output type")
        self.assertTrue(
            mgis_bv.getVariableTypeSymbolicRepresentation(o.type_identifier)
            == "tvector<3, real>", "invalid output type identifier")
        pass

    def test_behaviour_data(self):
        eps = 1.e-14
        e = numpy.asarray([1.3e-2, 1.2e-2, 1.4e-2, 0., 0., 0.],
                          dtype=numpy.float64)
        e2 = numpy.asarray([1.2e-2, 1.3e-2, 1.4e-2, 0., 0., 0.],
                           dtype=numpy.float64)
        h = mgis_bv.Hypothesis.Tridimensional
        b = self.__get_behaviour(h)
        m = mgis_bv.MaterialDataManager(b, 3)
        mgis_bv.setMaterialProperty(m.s1, "YoungModulus", 150e9)
        mgis_bv.setMaterialProperty(m.s1, "PoissonRatio", 0.3)
        mgis_bv.setExternalStateVariable(m.s1, "Temperature", 293.15)
        mgis_bv.update(m)
        # a different strain at each integration point
        for j in range(0, 3):
            m.s1.gradients[j, :] = (j + 1) * e
        outputs = numpy.empty(shape=(3, 3), dtype=numpy.float64)
        mgis_bv.executePostProcessing(outputs.reshape(9), m, "PrincipalStrain")
        for j in range(0, 3):
            for i in range(0, 3):
                self.assertTrue(
                    abs(outputs[j, i] - (j + 1) * e2[i]) < (j + 1) * eps,
                    f"invalid output value ({outputs[j, i]})")
        # range excluding the first integration point
        outputs[:] = -1
        mgis_bv.executePostProcessing(outputs.reshape(9), m, "PrincipalStrain",
                                      1, 3)
        for j in range(0, 3):
            for i in range(0, 3):
                o = -1 if j == 0 else (j + 1) * e2[i]
                self.assertTrue(
                    abs(outputs[j, i] - o) < (j + 1) * eps,
                    f"invalid output value ({outputs[j, i]})")
        # the results can not be written in a converted copy
        self.assertRaises(RuntimeError, mgis_bv.executePostProcessing,
                          numpy.empty(6, dtype=numpy.float32), m,
                          "PrincipalStrain")


if __name__ == '__main__':
    unittest.main()
