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
        self.assertTrue(len(postprocessings) == 1, "invalid number of post-processings")
        self.assertTrue("PrincipalStrain" in postprocessings, "invalid post-processing names")
        pass

if __name__ == '__main__':
    unittest.main()
