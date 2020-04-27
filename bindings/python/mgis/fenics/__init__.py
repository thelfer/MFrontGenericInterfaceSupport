#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mgis.fenics module

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.f
"""
from dolfin import *
parameters["form_compiler"]["representation"] = 'quadrature'
import warnings
from ffc.quadrature.deprecation import QuadratureRepresentationDeprecationWarning
warnings.simplefilter("ignore", QuadratureRepresentationDeprecationWarning)

from .nonlinear_material import MFrontNonlinearMaterial
from .nonlinear_problem import MFrontNonlinearProblem, MFrontOptimisationProblem, \
                                list_predefined_gradients
from .gradient_flux import *