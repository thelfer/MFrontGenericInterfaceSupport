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

The behaviour of solid materials is modelled using so-called
constitutive equations which describes how the internal state of the
material evolves with time. The knowledge of those internal state
variables allows the computation of local thermodynamic forces which
affects the material equilibrium at the structural scale. Du to the
large number of phenomena that needs to be described (plasticity,
viscoplasticy, damage, etc...), computational mechanics is one of most
demanding domain for advanced constitutive equations. The ability to
easily integrate user defined constitutive equations plays a major role
in the versatility of solvers.

The `MFront` code generator has been designed to simplify the process of
implementing constitutive equations. From a source file, `MFront`
generates `C++` code specific to many well-established (mostly
thermomechanical) solvers through dedicated interfaces: `Cast3M`,
`Code_Aster`, etc... 

The aim of this paper is to describe a side project of `MFront` called
`MFrontGenericInterfaceSupport` which aims at proving tools (functions,
classes, bindings, etc…) to handle behaviours written using `MFront`'
`generic` interface. Those tools are meant to be used by (FEM, FFT,
etc.) solver developers. Permissive licences have been chosen to allow
integration in open-source and proprietary codes. 

# Overview

# Available bindings

# Usage

## `FEniCS`

## `OpenGeoSys`

## `XPer`

## `JuliaFEM`

# Conclusions
