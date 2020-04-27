#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 16:35:18 2020

@author: bleyerj
"""
from dolfin import *
import ufl

def local_project(v, V, dx, u=None):
    """
    projects v on V with custom quadrature scheme dedicated to
    FunctionSpaces V of `Quadrature` type

    if u is provided, result is appended to u
    """
    v_ = TestFunction(V)
    dv = TrialFunction(V)
    a_proj = inner(dv, v_)*dx
    b_proj = inner(v, v_)*dx
    solver = LocalSolver(a_proj, b_proj, LocalSolver.SolverType.Cholesky)
    solver.factorize()
    if u is None:
        u = Function(V)
        solver.solve_local_rhs(u)
        return u
    else:
        solver.solve_local_rhs(u)
        return

def symmetric_tensor_to_vector(T, T22=0):
    """ Return symmetric tensor components in vector form notation following MFront conventions
        T22 can be specified when T is only (2,2)"""
    if ufl.shape(T)==(2, 2):
        return as_vector([T[0, 0], T[1, 1], T22, sqrt(2)*T[0, 1]])
    elif ufl.shape(T)==(3, 3):
        return as_vector([T[0, 0], T[1, 1], T[2, 2], sqrt(2)*T[0, 1], sqrt(2)*T[0, 2], sqrt(2)*T[1, 2]])
    elif len(ufl.shape(T)) == 1:
        return T
    else:
        raise NotImplementedError

def nonsymmetric_tensor_to_vector(T, T22=0):
    """ Return nonsymmetric tensor components in vector form notation following MFront conventions
        T22 can be specified when T is only (2,2) """
    if ufl.shape(T)==(2, 2):
         return as_vector([T[0, 0], T[1, 1], T22, T[0, 1], T[1, 0]])
    elif ufl.shape(T)==(3, 3):
        return as_vector([T[0, 0], T[1, 1], T[2, 2], T[0, 1], T[1, 0], T[0, 2], T[2, 0], T[1, 2], T[2, 1]])
    elif len(ufl.shape(T)) == 1:
        return T
    else:
        raise NotImplementedError

def axi_grad(r, v):
    """
        Axisymmetric gradient in cylindrical coordinate (er, etheta, ez) for:
        * a scalar v(r, z)
        * a 2d-vectorial (vr(r,z), vz(r, z))
        * a 3d-vectorial (vr(r,z), 0, vz(r, z))
    """
    if ufl.shape(v)==(3,):
        return as_matrix([[v[0].dx(0), -v[1]/r, v[0].dx(1)],
                           [v[1].dx(0), v[0]/r, v[1].dx(1)],
                           [v[2].dx(0), 0, v[2].dx(1)]])
    elif ufl.shape(v)==(2,):
        return as_matrix([[v[0].dx(0), v[0].dx(1), 0],
                           [v[1].dx(0), v[1].dx(1), 0],
                           [0, 0, v[0]/r]])
    elif ufl.shape(v)==():
        return as_vector([v.dx(0), 0, v.dx(1)])
    else:
        raise NotImplementedError

def symmetric_gradient(g):
    """ Return symmetric gradient components in vector form"""
    return symmetric_tensor_to_vector(sym(g))

def transformation_gradient(g, dim=3):
    """ Return transformation gradient components in vector form"""
    return nonsymmetric_tensor_to_vector(Identity(dim)+g, T22=1)

def gradient(g):
    """ Return displacement gradient components in vector form"""
    return nonsymmetric_tensor_to_vector(g)

def get_quadrature_element(cell, degree, dim=0):
    if dim in [0, 1, (), (0,), (1,), (0, 0), (0, 1), (1, 0), (1, 1)]:
        return FiniteElement("Quadrature", cell, degree=degree,
                            quad_scheme='default')
    elif type(dim) == int:
        return VectorElement("Quadrature", cell, degree=degree,
                              dim=dim, quad_scheme='default')
    elif len(dim)==1 or (len(dim)==2 and 0 in dim):
        d = [dd for dd in list(dim) if dd != 0][0]
        return VectorElement("Quadrature", cell, degree=degree,
                              dim=d, quad_scheme='default')
    elif type(dim) == tuple:
        return TensorElement("Quadrature", cell, degree=degree, shape=dim, quad_scheme='default')
    else:
        raise ValueError("Wrong shape for dim=", dim)