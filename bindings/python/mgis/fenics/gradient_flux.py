#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Gradient and Flux classes

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.f
"""
from dolfin import *
import numpy as np
import ufl
from .utils import local_project, symmetric_tensor_to_vector, \
                nonsymmetric_tensor_to_vector, get_quadrature_element, \
                vector_to_tensor

class QuadratureFunction:
    """An abstract class for functions defined at quadrature points"""
    def __init__(self, name, shape):
        self.shape = shape
        self.name = name

    def initialize_function(self, mesh, quadrature_degree):
        self.quadrature_degree = quadrature_degree
        self.mesh = mesh
        self.dx = Measure("dx", metadata={"quadrature_degree": quadrature_degree})
        We = get_quadrature_element(mesh.ufl_cell(), quadrature_degree, self.shape)
        self.function_space = FunctionSpace(mesh, We)
        self.function = Function(self.function_space, name=self.name)

    def update(self, x):
        self.function.vector().set_local(x)

    def project_on(self, space, degree, as_tensor=False):
        """
        Returns the projection on a standard CG/DG space

        Parameters
        ----------

        space: str
            FunctionSpace type ("CG", "DG",...)
        degree: int
            FunctionSpace degree
        as_tensor: bool
            Returned as a tensor if True, in vector notation otherwise
        """
        fun = self.function
        if as_tensor:
            fun = vector_to_tensor(self.function)
            shape = ufl.shape(fun)
            V = TensorFunctionSpace(self.mesh, space, degree, shape=shape)
        elif self.shape == 1:
            V = FunctionSpace(self.mesh, space, degree)
        else:
            V = VectorFunctionSpace(self.mesh, space, degree, dim=self.shape)
        v = Function(V, name=self.name)
        v.assign(project(fun, V,
                         form_compiler_parameters={"quadrature_degree": self.quadrature_degree}))
        return v

class Gradient(QuadratureFunction):
    """
    This class provides a representation of MFront gradient objects. Its main
    purpose is to provide the corresponding UFL expression, linking MFront and
    FEniCS concepts. It also handles:
        - the reshaping from UFL tensorial representation to MFront
          vectorial conventions
        - the symbolic expression of the gradient variation
          (directional derivative)
        - the representation as a Quadrature function

    This class is intended for internal use only. Gradient objects must be
    declared by the user using the registration concept.
    """
    def __init__(self, variable, expression, name, symmetric=None):
        self.variable = variable
        if symmetric is None:
            self.expression = expression
        elif symmetric:
            if self.variable.geometric_dimension() == 2:
                self.expression = as_vector([symmetric_tensor_to_vector(expression)[i] for i in range(4)])
            else:
                self.expression = symmetric_tensor_to_vector(expression)
        else:
            if self.variable.geometric_dimension() == 2:
                self.expression = as_vector([nonsymmetric_tensor_to_vector(expression)[i] for i in range(5)])
            else:
                self.expression = nonsymmetric_tensor_to_vector(expression)
        shape = ufl.shape(self.expression)
        if len(shape) == 1:
            self.shape = shape[0]
        elif shape == ():
            self.shape = 0
        else:
            self.shape = shape
        self.name = name

    def __call__(self, v):
        return ufl.replace(self.expression, {self.variable: v})

    def variation(self, v):
        """ Directional derivative in direction v """
        # return ufl.algorithms.expand_derivatives(ufl.derivative(self.expression, self.variable, v))
        deriv = sum([ufl.derivative(self.expression, var, v_) for (var, v_) in zip(split(self.variable), split(v))])
        return ufl.algorithms.expand_derivatives(deriv)

    def initialize_function(self, mesh, quadrature_degree):
        super().initialize_function(mesh, quadrature_degree)
        self.update()

    def update(self, x=None):
        if x is None:
            self._evaluate_at_quadrature_points(self.expression)
        else:
            if isinstance(x, np.ndarray):
                self.function.vector().set_local(x)
            else:
                self._evaluate_at_quadrature_points(x)

    def _evaluate_at_quadrature_points(self, x):
        local_project(x, self.function_space, self.dx, self.function)


class Var(Gradient):
    """ A simple variable """
    # def __init__(self, variable, name):
    #     return Gradient.__init__(self, variable, variable, name)

    def _evaluate_at_quadrature_points(self, x):
        local_project(x, self.function_space, self.dx, self.function)

class QuadratureFunctionTangentBlocks(QuadratureFunction):
    """An abstract class for Flux and InternalStateVariables"""
    def initialize_tangent_blocks(self, variables):
        self.variables = variables
        values = [Function(FunctionSpace(self.mesh,
                          get_quadrature_element(self.mesh.ufl_cell(),
                          self.quadrature_degree, dim=(self.shape, v.shape))),
                           name="d{}_d{}".format(self.name, v.name))
                          for v in self.variables]
        keys = [v.name for v in self.variables]
        self.tangent_blocks = dict(zip(keys, values))

class Flux(QuadratureFunctionTangentBlocks):
    pass

class InternalStateVariable(QuadratureFunctionTangentBlocks):
    pass
