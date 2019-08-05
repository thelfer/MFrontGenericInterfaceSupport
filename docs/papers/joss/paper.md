---
title: The `MFrontGenericInterfaceSupport` project
tags:
  - MFront
  - Mechanical solvers
  - Mechanical behaviours
authors:
 - name: Thomas Helfer
   orcid: 0000-0003-2460-5816
   affiliation: 1
 - name: Tero Frondelius
   orcid: 0000-0003-0872-7098
   affiliation: 2
affiliations:
 - name: French Alternative Energies and Atomic Energy Commission
   index: 1
 - name: Wärtsilä
   index: 2
date: 2 August 2019
bibliography: bibliography.bib
---

# Introduction

The behaviour of solid materials is modelled using so-called
constitutive equations which describes how the internal state of the
material evolves with time. The knowledge of those internal state
variables allows the computation of local thermodynamic forces which
affects the material equilibrium at the structural scale. Du to the
large number of phenomena that needs to be described (plasticity,
viscoplasticy, damage, etc...), computational mechanics is one of most
demanding domain for advanced constitutive equations.

The ability to easily integrate user defined constitutive equations
plays a major role in the versatility of solvers.

The `MFront` code generator has been designed to simplify the process of
implementing constitutive equations
[@helfer_introducing_2015;@cea_mfront_2019]. From a source file,
`MFront` generates `C++` code specific to many well-established (mostly
thermomechanical) solvers through dedicated interfaces: `Cast3M`,
`Code_Aster`, etc...

The aim of this paper is to describe a side project of `MFront` called
`MFrontGenericInterfaceSupport` which aims at proving tools (functions,
classes, bindings, etc…) to handle behaviours written using `MFront`'
`generic` interface [@helfer_mgis_2019]. Those tools are meant to be
used by (FEM, FFT, etc.) solver developers. Permissive licences have
been chosen to allow integration in open-source and proprietary codes.

# Overview

The expected workflow is meant to be:

1. The user implements its constitutive equations using one of the
  domain specific languages provided by `MFront` and generates a set of
  `C++` files using the `generic` interface and compiles those sources
  into a shared library.
2. In the solver input files, the user indicates that he/she wants to
  use this behaviour.

At this stage, the solver must load the shared library and the function
implementing the behaviour. During the pre-processing stage, i.e. before
actually running the simulation, the solver must:

- Check if the behaviour is consistent with the computation to be
  performed. For example, one shall not use a finite strain behaviours
  in a small strain analysis.
- Check if all the information required by the behaviour are available.
  In practice, a behaviour can require that the solver provides a set of
  values describing the material to be treated: those value are called
  material properties. In general, those material properties must be
  given by the user in the input file. For multi-physics modelling, a
  behaviour may also need to know how

computation the memory holding the internal state variables at each
  integration points.
- allocate the memory holding the internal state variables at each
  integration points.

# Available bindings

# Usage

## `FEniCS`

## `OpenGeoSys`

## `XPer`

## `JuliaFEM`

# Conclusions
