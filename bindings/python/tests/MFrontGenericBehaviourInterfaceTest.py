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
        b = mgis_bv.load(lib,'Gurson',h)
        self.assertTrue(b.behaviour == "Gurson", "invalid behaviour name")
        self.assertTrue(b.hypothesis == h, "invalid hypothesis")
        self.assertTrue(b.source == "Gurson.mfront", "invalid source")
        self.assertTrue(b.tfel_version == version, "invalid TFEL version")
    
if __name__ == '__main__':
    unittest.main()
