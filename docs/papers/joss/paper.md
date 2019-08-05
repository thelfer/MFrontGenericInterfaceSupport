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

<!--
pandoc -f markdown_strict --bibliography=bibliography.bib --filter pandoc-citeproc paper.md -o paper.pdf
-->

# Introduction

The behaviour of solid materials is modelled using so-called
constitutive equations which describes how the internal state of the
material evolves with time. The knowledge of those internal state
variables allows the computation of local thermodynamic forces which
affects the material equilibrium at the structural scale.

Du to the large number of phenomena that needs to be described
(plasticity, viscoplasticy, damage, etc...), computational mechanics is
one of most demanding domain for advanced constitutive equations.

The ability to easily integrate user defined constitutive equations
plays a major role in the versatility of (mechanical) solvers^[Here we
use the term solver to emphasize that the numerical method used to
discretize the equilibrium equations is not significant here (Finite
Element Method (FEM), Fast Fourier Transform (FFT), etc.).].

The `MFront` open-source code generator has been designed to simplify
the process of implementing constitutive equations
[@helfer_introducing_2015;@cea_mfront_2019]. From a source file,
`MFront` generates `C++` code specific to many well-established (mostly
thermomechanical) solvers through dedicated interfaces: `Cast3M`,
`Code_Aster`, etc... Those sources are then compiled into shared
libraries.

`MFront` recently introduced a so-called `generic` interface. The aim of
this paper is to describe the `MFrontGenericInterfaceSupport` project
(denoted `MGIS` in the following), which aims at proving tools (functions,
classes, bindings, etc…) to handle behaviours generated using `MFront`'
`generic` interface [@helfer_mgis_2019]. Those tools are meant to
alleviate the work required by solvers' developers. Permissive licences
have been chosen to allow integration in open-source and proprietary
codes.

This paper is divided in three parts:

1. Section 1 gives a brief overiew of `MGIS`.
2. Section 2 describes the various bindings available.
3. Section 3 describes some example of usage in various open-source
  solvers: `FEniCS`, `OpenGeoSys` and `JuliaFEM`.

# Overview

The aims of the `MFrontGenericInterfaceSupport` project is twofold:

1. At the pre-processing state, allow retrieving meta data about a
  particular behaviours and proper memory allocation. At the
  post-processing stage, ease access to internal state variables. This
  is described in Section @sec:prepost.
2. During computations, simplify the integration the behaviour at
  integration points^[Integration points is used here as a generic
  placeholder. When using FFT for solving the equilibrium equations, the
  integration points are voxels. When using FEM, the integrations points
  are the usual Gauss points of the elements. etc.] and the update of
  the internal state variables.

## Preprocessing and post-processing stages{#sec:prepost}

### Retrieving meta-data

# Available bindings

# Example of usage

## `FEniCS`

## `OpenGeoSys`

## `JuliaFEM`

# Conclusions


# Acknowledgements

This research was conducted in the framework of the `PLEIADES`
project,which is supported financially by the CEA (Commissariat à
l’Energie Atomique et aux Energies Alternatives), EDF (Electricité de
France) and Framatome.Acknowledgements

# References