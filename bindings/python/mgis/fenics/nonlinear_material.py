#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MFrontNonlinearMaterial class

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.f
"""
import mgis.behaviour as mgis_bv
from .gradient_flux import Var
import dolfin
import subprocess
import os

mgis_hypothesis = {"plane_strain": mgis_bv.Hypothesis.PlaneStrain,
                   "plane_stress": mgis_bv.Hypothesis.PlaneStress,
                   "3d": mgis_bv.Hypothesis.Tridimensional,
                   "axisymmetric": mgis_bv.Hypothesis.Axisymmetrical}


class MFrontNonlinearMaterial:
    """
    This class handles the loading of a MFront behaviour through MFrontGenericInterfaceSupport.
    """
    def __init__(self, path, name, hypothesis="3d",
                 material_properties={}, parameters={}):
        """
        Parameters
        -----------

        path : str
            path to the 'libMaterial.so' library containing MFront material laws
        name : str
            name of the MFront behaviour
        hypothesis : {"plane_strain", "3d", "axisymmetric"}
            modelling hypothesis
        material_properties : dict
            a dictionary of material properties. The dictionary keys must match
            the material property names declared in the MFront behaviour. Values
            can be constants or functions.
        parameters : dict
            a dictionary of parameters. The dictionary keys must match the parameter
            names declared in the MFront behaviour. Values must be constants.
        """
        self.path = path
        self.name = name
        # Defining the modelling hypothesis
        self.hypothesis = mgis_hypothesis[hypothesis]
        self.material_properties = material_properties
        # Loading the behaviour
        try:
            self.load_behaviour()
        except:
            cwd = os.getcwd()
            install_path = "/".join(path.split("/")[:-2])+"/"
            os.chdir(install_path)
            print("Behaviour '{}' has not been found in '{}'.".format(self.name, self.path))
            print("Attempting to compile '{}.mfront' in '{}'...".format(self.name, install_path))
            subprocess.run(["mfront", "--obuild", "--interface=generic", self.name+".mfront"])
            os.chdir(cwd)
            self.load_behaviour()
        self.update_parameters(parameters)

    def load_behaviour(self):
        self.is_finite_strain = mgis_bv.isStandardFiniteStrainBehaviour(self.path, self.name)
        if self.is_finite_strain:
            # finite strain options
            bopts = mgis_bv.FiniteStrainBehaviourOptions()
            bopts.stress_measure = mgis_bv.FiniteStrainBehaviourOptionsStressMeasure.PK1
            bopts.tangent_operator = mgis_bv.FiniteStrainBehaviourOptionsTangentOperator.DPK1_DF
            self.behaviour = mgis_bv.load(bopts, self.path, self.name, self.hypothesis)
        else:
            self.behaviour = mgis_bv.load(self.path, self.name, self.hypothesis)

    def set_data_manager(self, ngauss):
        # Setting the material data manager
        self.data_manager = mgis_bv.MaterialDataManager(self.behaviour, ngauss)
        self.update_material_properties()

    def update_parameters(self, parameters):
        for (key, value) in parameters.items():
            self.behaviour.setParameter(key, value)

    def update_material_properties(self, material_properties=None):
        if material_properties is not None:
            self.material_properties = material_properties
        for s in [self.data_manager.s0, self.data_manager.s1]:
            for (key, value) in self.material_properties.items():
                if type(value) in [int, float]:
                    mgis_bv.setMaterialProperty(s, key, value)
                else:
                    if isinstance(value, dolfin.Function):
                        value = value.vector().get_local()
                    mgis_bv.setMaterialProperty(s, key, value, mgis_bv.MaterialStateManagerStorageMode.LocalStorage)

    def update_external_state_variables(self, external_state_variables):
        for s in [self.data_manager.s0, self.data_manager.s1]:
            for (key, value) in external_state_variables.items():
                if type(value) in [int, float]:
                    mgis_bv.setExternalStateVariable(s, key, value)
                elif isinstance(value, dolfin.Constant):
                    mgis_bv.setExternalStateVariable(s, key, float(value))
                else:
                    if isinstance(value, dolfin.Function):
                        values = value.vector().get_local()
                    elif isinstance(value, Var):
                        value.update()
                        values = value.function.vector().get_local()
                    else:
                        print(isinstance(value, Var))
                        print(value)
                        raise NotImplementedError("{} type is not supported for external state variables".format(type(value)))
                    mgis_bv.setExternalStateVariable(s, key, values, mgis_bv.MaterialStateManagerStorageMode.LocalStorage)

    def get_parameter(self, name):
        return self.behaviour.getParameterDefaultValue(name)

    def get_parameter_names(self):
        return self.behaviour.params

    def get_material_property_names(self):
        return [svar.name for svar in self.behaviour.mps]

    def get_external_state_variable_names(self):
        return [svar.name for svar in self.behaviour.external_state_variables]

    def get_internal_state_variable_names(self):
        return [svar.name for svar in self.behaviour.internal_state_variables]

    def get_gradient_names(self):
        return [svar.name for svar in self.behaviour.gradients]

    def get_flux_names(self):
        return [svar.name for svar in self.behaviour.thermodynamic_forces]

    def get_material_property_sizes(self):
        return [mgis_bv.getVariableSize(svar, self.hypothesis) for svar in self.behaviour.mps]

    def get_external_state_variable_sizes(self):
        return [mgis_bv.getVariableSize(svar, self.hypothesis) for svar in self.behaviour.external_state_variables]

    def get_internal_state_variable_sizes(self):
        return [mgis_bv.getVariableSize(svar, self.hypothesis) for svar in self.behaviour.internal_state_variables]

    def get_gradient_sizes(self):
        return [mgis_bv.getVariableSize(svar, self.hypothesis) for svar in self.behaviour.gradients]

    def get_flux_sizes(self):
        return [mgis_bv.getVariableSize(svar, self.hypothesis) for svar in self.behaviour.thermodynamic_forces]

    def get_tangent_block_names(self):
        return [(t[0].name, t[1].name) for t in self.behaviour.tangent_operator_blocks]

    def get_tangent_block_sizes(self):
        return [tuple([mgis_bv.getVariableSize(tt, self.hypothesis) for tt in t]) \
                for t in self.behaviour.tangent_operator_blocks]