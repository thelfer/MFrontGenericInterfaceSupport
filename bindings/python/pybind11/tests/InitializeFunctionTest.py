# -*- coding: utf-8 -*-

import os
import numpy

try:
    import unittest2 as unittest
except ImportError:
    import unittest
import mgis.behaviour as mgis_bv


class InitalizeFunctionTest(unittest.TestCase):
    def __get_behaviour(self, h):
        lib = os.environ['MGIS_TEST_BEHAVIOURS_LIBRARY']
        return mgis_bv.load(lib, 'InitializeFunctionTest', h)

    def test_behaviour(self):
        h = mgis_bv.Hypothesis.Tridimensional
        b = self.__get_behaviour(h)
        # path to the test library
        version = os.environ['MGIS_TEST_TFEL_VERSION']
        self.assertTrue(b.behaviour == "InitializeFunctionTest",
                        "invalid behaviour name")
        self.assertTrue(b.hypothesis == h, "invalid hypothesis")
        self.assertTrue(b.source == "InitializeFunctionTest.mfront",
                        "invalid source")
        self.assertTrue(b.tfel_version == version, "invalid TFEL version")
        ifcts = b.getInitializeFunctionsNames()
        self.assertTrue(
            len(ifcts) == 2, "invalid number of initialization functions")
        self.assertTrue("StressFromInitialPressure" in ifcts,
                        "invalid initialization functions names")
        self.assertTrue("ElasticStrainFromInitialStress" in ifcts,
                        "invalid initialization functions names")
        pass

    def test_behaviour_data(self):
        pr = -1.2e5
        eps = -pr * 1.e-14
        h = mgis_bv.Hypothesis.Tridimensional
        b = self.__get_behaviour(h)
        d = mgis_bv.BehaviourData(b)
        mgis_bv.setExternalStateVariable(d.s0, "Temperature", 293.15)
        mgis_bv.setExternalStateVariable(d.s1, "Temperature", 293.15)
        inputs = numpy.asarray([pr], dtype=numpy.float64)
        v = mgis_bv.make_view(d)
        mgis_bv.executeInitializeFunction(v, b, "StressFromInitialPressure",
                                          inputs)
        # Due to restrictions of the BehaviourDataView class (the
        # initial state is assumed immutable), only the state at the
        # end of the time step is initialized. Calling update is thus
        # required in most cases
        mgis_bv.update(d)
        s0 = d.s0.thermodynamic_forces
        s1 = d.s1.thermodynamic_forces
        for i in range(0, 3):
            self.assertTrue(
                abs(s0[i] - pr) < eps,
                "invalid stress value at the beginning of the time step")
            self.assertTrue(
                abs(s1[i] - pr) < eps,
                "invalid stress value at the end of the time step")
        for i in range(3, 6):
            self.assertTrue(
                abs(s0[i]) < eps,
                "invalid stress value at the beginning of the time step")
            self.assertTrue(
                abs(s1[i]) < eps,
                "invalid stress value at the end of the time step")

    def test_behaviour_data2(self):
        E = 200e9
        nu = 0.3
        sxx = 150e6
        eps = 10 * sxx * 1e-14
        h = mgis_bv.Hypothesis.Tridimensional
        b = self.__get_behaviour(h)
        d = mgis_bv.BehaviourData(b)
        d.s0.thermodynamic_forces[:] = [sxx, 0, 0, 0, 0, 0]
        mgis_bv.setExternalStateVariable(d.s0, "Temperature", 293.15)
        mgis_bv.setExternalStateVariable(d.s1, "Temperature", 293.15)
        v = mgis_bv.make_view(d)
        mgis_bv.executeInitializeFunction(v, b,
                                          "ElasticStrainFromInitialStress")
        # Due to restrictions of the BehaviourDataView class (the
        # initial state is assumed immutable), only the state at the
        # end of the time step is initialized. Calling update is thus
        # required in most cases
        mgis_bv.update(d)
        eel_values = [sxx / E, -nu * sxx / E, -nu * sxx / E, 0, 0, 0]
        eel0 = d.s0.internal_state_variables
        eel1 = d.s1.internal_state_variables
        for i in range(0, 6):
            self.assertTrue(
                abs(eel0[i] - eel_values[i]) < eps,
                "invalid elastic strain value at the beginning of the time step"
            )
            self.assertTrue(
                abs(eel1[i] - eel_values[i]) < eps,
                "invalid elastic strain value at the end of the time step")

    def test_data_manager(self):
        """
        Call an initialize function with uniform inputs
        """
        pr = -1.2e5
        eps = -pr * 1.e-14
        h = mgis_bv.Hypothesis.Tridimensional
        b = self.__get_behaviour(h)
        d = mgis_bv.MaterialDataManager(b, 2)
        mgis_bv.setExternalStateVariable(d.s0, "Temperature", 293.15)
        mgis_bv.setExternalStateVariable(d.s1, "Temperature", 293.15)
        inputs = numpy.asarray([pr], dtype=numpy.float64)
        mgis_bv.executeInitializeFunction(d, "StressFromInitialPressure",
                                          inputs)
        # Due to restrictions of the BehaviourDataView class (the
        # initial state is assumed immutable), only the state at the
        # end of the time step is initialized. Calling update is thus
        # required in most cases
        mgis_bv.update(d)
        s0 = d.s0.thermodynamic_forces
        s1 = d.s1.thermodynamic_forces
        for n in range(0, 2):
            for i in range(0, 3):
                self.assertTrue(
                    abs(s0[n][i] - pr) < eps,
                    "invalid stress value at the beginning of the time step")
                self.assertTrue(
                    abs(s1[n][i] - pr) < eps,
                    "invalid stress value at the end of the time step")
            for i in range(3, 6):
                self.assertTrue(
                    abs(s0[n][i]) < eps,
                    "invalid stress value at the beginning of the time step")
                self.assertTrue(
                    abs(s1[n][i]) < eps,
                    "invalid stress value at the end of the time step")

    def test_data_manager2(self):
        """
        Call an initialize function with non uniform inputs
        """
        pr = [-1.2e5, 0.7e5]
        eps = -pr[0] * 1.e-14
        h = mgis_bv.Hypothesis.Tridimensional
        b = self.__get_behaviour(h)
        d = mgis_bv.MaterialDataManager(b, 2)
        mgis_bv.setExternalStateVariable(d.s0, "Temperature", 293.15)
        mgis_bv.setExternalStateVariable(d.s1, "Temperature", 293.15)
        inputs = numpy.asarray(pr, dtype=numpy.float64)
        mgis_bv.executeInitializeFunction(d, "StressFromInitialPressure",
                                          inputs)
        # Due to restrictions of the BehaviourDataView class (the
        # initial state is assumed immutable), only the state at the
        # end of the time step is initialized. Calling update is thus
        # required in most cases
        mgis_bv.update(d)
        s0 = d.s0.thermodynamic_forces
        s1 = d.s1.thermodynamic_forces
        for n in range(0, 2):
            for i in range(0, 3):
                self.assertTrue(
                    abs(s0[n][i] - pr[n]) < eps,
                    "invalid stress value at the beginning of the time step")
                self.assertTrue(
                    abs(s1[n][i] - pr[n]) < eps,
                    "invalid stress value at the end of the time step")
            for i in range(3, 6):
                self.assertTrue(
                    abs(s0[n][i]) < eps,
                    "invalid stress value at the beginning of the time step")
                self.assertTrue(
                    abs(s1[n][i]) < eps,
                    "invalid stress value at the end of the time step")

    def test_material_data_manager3(self):
        E = 200e9
        nu = 0.3
        sxx = 150e6
        eps = 10 * sxx * 1e-14
        h = mgis_bv.Hypothesis.Tridimensional
        b = self.__get_behaviour(h)
        d = mgis_bv.MaterialDataManager(b, 2)
        d.s0.thermodynamic_forces[:] = [[sxx, 0, 0, 0, 0, 0],
                                        [sxx, 0, 0, 0, 0, 0]]
        mgis_bv.setExternalStateVariable(d.s0, "Temperature", 293.15)
        mgis_bv.setExternalStateVariable(d.s1, "Temperature", 293.15)
        mgis_bv.executeInitializeFunction(d, "ElasticStrainFromInitialStress")
        # Due to restrictions of the BehaviourDataView class (the
        # initial state is assumed immutable), only the state at the
        # end of the time step is initialized. Calling update is thus
        # required in most cases
        mgis_bv.update(d)
        eel_values = [sxx / E, -nu * sxx / E, -nu * sxx / E, 0, 0, 0]
        eel0 = d.s0.internal_state_variables
        eel1 = d.s1.internal_state_variables
        for n in range(0, 2):
            for i in range(0, 6):
                self.assertTrue(
                    abs(eel0[n][i] - eel_values[i]) < eps,
                    "invalid elastic strain value at the beginning of the time step"
                )
                self.assertTrue(
                    abs(eel1[n][i] - eel_values[i]) < eps,
                    "invalid elastic strain value at the end of the time step")


if __name__ == '__main__':
    unittest.main()
