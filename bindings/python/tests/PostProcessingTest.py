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
        pass

    def test_behaviour_data(self):
        eps = 1.e-14
        e = numpy.asarray([1.3e-2, 1.2e-2, 1.4e-2, 0., 0., 0.],
                          dtype=numpy.float64)
        e2 = numpy.asarray([1.2e-2, 1.3e-2, 1.4e-2, 0., 0., 0.],
                           dtype=numpy.float64)
        h = mgis_bv.Hypothesis.Tridimensional
        b = self.__get_behaviour(h)
        m = mgis_bv.MaterialDataManager(b, 2)
        mgis_bv.setMaterialProperty(m.s1, "YoungModulus", 150e9)
        mgis_bv.setMaterialProperty(m.s1, "PoissonRatio", 0.3)
        mgis_bv.setExternalStateVariable(m.s1, "Temperature", 293.15)
        mgis_bv.update(m)
        m.s1.gradients[0:] = e
        m.s1.gradients[1:] = e
        print(m.s1.gradients)
        outputs = numpy.empty(shape=(2, 3), dtype=numpy.float64)
        mgis_bv.executePostProcessing(outputs.reshape(6), m, "PrincipalStrain")
        for i in range(0, 3):
            self.assertTrue(
                abs(outputs[0, i] - e2[i]) < eps,
                f"invalid output value ({outputs[0, i]} vs {e2[i]})")
            self.assertTrue(
                abs(outputs[1, i] - e2[i]) < eps,
                f"invalid output value ({outputs[1, i]} vs {e2[i]})")


if __name__ == '__main__':
    unittest.main()
