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
   orcid: 0000-0003-2288-0902
   affiliation: "2, 3"
 - name: Thomas Nagel
   orcid: 0000-0001-8459-4616
   affiliation: 4
 - name: Dmitri Naumov
   affiliation: 5
affiliations:
 - name: French Alternative Energies and Atomic Energy Commission
   index: 1
 - name: Wärtsilä
   index: 2
 - name: University of Oulu
   index: 3
 - name: Technische Universität Bergakademie Freiberg
   index: 4
 - name: Helmholtz Centre for Environmental Research -- UFZ
   index: 5
date: 2 August 2019
bibliography: bibliography.bib
---

<!--
pandoc -f markdown_strict --bibliography=bibliography.bib --filter pandoc-citeproc paper.md -o paper.pdf
-->

# Introduction

The behaviour of solid materials is modelled using so-called
constitutive equations which describe how the internal state of the
material evolves with changing external conditions. Those state variables can describe many
microstructural aspects of the material: grain size, dislocation
density, hardening state, etc. or be considered purely phenomenologically.

The knowledge of those internal state variables allows the computation
of local thermodynamic forces which affects the material equilibrium at
the structural scale.
<!-- refine later -->

Due to the large number of phenomena that needs to be described
(plasticity, viscoplasticy, damage, etc...) for a tremendous range of natural
and man-made materials, computational mechanics is
one of most demanding domains for advanced constitutive equations.

The ability to easily integrate user-defined constitutive equations
plays a major role in the versatility of (mechanical) solvers^[Here we
use the term solver to emphasize that the numerical method used to
discretize the equilibrium equations is not significant here (Finite
Element Method (FEM), Fast Fourier Transform (FFT), etc.).].

The `MFront` open-source code generator has been designed to simplify
the process of implementing constitutive equations
[@helfer_introducing_2015;@cea_mfront_2019]. From a source file,
`MFront` generates `C++` code specific to many well-established (mostly
thermo-mechanical) solvers through dedicated interfaces: `Cast3M`,
`Code_Aster`, etc... Those sources are then compiled into shared
libraries.

`MFront` recently introduced a so-called `generic` interface. The aim of
this paper is to describe the `MFrontGenericInterfaceSupport` project
(denoted `MGIS` in the following), which aims at proving tools
(functions, classes, bindings, etc.) to handle constitutive equations
generated using `MFront`' `generic` interface [@helfer_mgis_2019]. Those
tools are meant to alleviate the work required by solvers' developers.
Permissive licences have been chosen to allow integration in open-source
and proprietary codes.

This paper is divided in three parts:

1. Section 1 gives a brief overiew of `MGIS`.
2. Section 2 describes the various bindings available.
3. Section 3 describes some example of usage in various open-source
  solvers: `FEniCS`, `OpenGeoSys` and `JuliaFEM`.

# Overview

The aims of the `MFrontGenericInterfaceSupport` project is twofold:

1. At the pre-processing state, allow retrieving meta data about a
  particular behaviours and proper memory allocation. At the
  post-processing stage, ease access to internal state variables.
2. During computations, simplify the integration the behaviour at
  integration points^[Integration points is used here as a generic
  placeholder. When using FFT for solving the equilibrium equations, the
  integration points are voxels. When using FEM, the integrations points
  are the usual Gauss points of the elements. etc.] and the update of
  the internal state variables.

## Preprocessing and post-processing stages{#sec:prepost}

When dealing with user defined behaviours, most solvers, including
`Abaqus/Standard` for example, require that the user describe the
behaviour in the input file (e.g. number of internal state variables) and
take care of the consistency of the behaviour with the hypothesis made
during the computation (e.g. a finite strain behaviour must be used in a
finite strain analysis based on the appropriate deformation and stress measures
as well as reference configurations). This is an error-prone task and may lead to spurious
and/or undefined behaviour in case of conceptual or coding errors. Most of the time, the 

In `MGIS`, the logic is quite different: the user must only declare the
shared library, the behaviour and the modelling hypothesis
(tridimensional, plane strain, etc.). With this information, the
library can retrieve various metadata from the shared library which
fully describe how to interact with the behaviour. The solver using
`MGIS` can then check if the behaviour is consistent with the
computations to be performed.

The metadata can also be used to allocate the memory required to store
the state of the material at each integration point. Note that `MGIS`
has been designed to allow the following storage types of the state:

- an `MGIS` data structure per integration point. While this causes
  memory fragmentation, this is the most frequent choice. The memory is
  automatically allocated by `MGIS`.
- an `MGIS` data structure that stores the states of an arbitrary number
  of integration points. `MGIS` can allocate the memory associated with
  the state of all specified integrations points or borrow memory allocated
  by the solver.

For post-processing a set of functions are provided to retrieve
information about the state of the material. For example, one can
retrieve the value of a state variable from the previous data
structures.

## Computation stage

At each time step, the constitutive equations must be integrated to get
the state of the material at the end of the time step based on alterations of 
corresponding driving forces. As most phenomena
are nonlinear, an iterative scheme is required at the equilibrium scale
to find the local loading of the material: the integration of the
constitutive equations is thus called several times with different
estimates of the loading of the material.

`MGIS` provides a function to integrate the constitutive equations at one
integration point or on a set of integration points^[This strongly
depends on the data structure chosen to store the internal state
variables.].

The integration of the constitutive equations at different integration
points are usually independent: thus, when handling a set of integration points,
`MGIS` can parallelize the integrations using a granularity chosen by
the solver.

# Main language and available bindings

`MGIS` is written in `C++-11`. The `C++` API is described in another
report, see [@helfer_brief_2019].

The following bindings are available:

- `python`.
- `Julia`.
- `Fortran 2003`.
- `C`.

# Example of usage

## `FEniCS`

!["Large strain elasto-plastic modelling of a notched
bar"](img/FEniCS.png "Large strain elasto-plastic modelling of a notched
bar")

`FEniCS` is a popular open-source computing platform for solving partial
differential equations [@logg_automated_2012;@alnaes_fenics_2015]. Non
linear mechanics computations combining `FEniCS` at the equilibrium
scale and `MFront` to describe the constitutive equations can be
performed through the `python` bindings of `MGIS` as demonstrated by
Bleyer et al (see [@bleyer_elasto-plastic_2019;@bleyer_fenics_2019]).
Extensions to finite strain elasto-plasticity as been recently added as
shown on Figure 1 which models a tensile test on a notched bar^[This
case is adapted from a non-regression test of `Code_Aster` finite
element solver, see @edf_ssna303_2011 for details].

## `OpenGeoSys`

[@kolditz_opengeosys:_2012]
[@Bilke2019]

OpenGeoSys (OGS) is a scientific open-source initiative for the numerical simulation of thermo-hydro-mechanical/
chemical (THMC) processes in porous and fractured media, inspired by FEFLOW and ROCKFLOW concepts
and continuously developed since the mid-eighties, see Kolditz (1990), Wollrath (1990), Kroehn (1991) and Helmig
(1993).

The OGS framework is targeting applications in environmental geoscience, e.g., in the fields of contaminant
hydrology (Walther et al, 2014), water resources and waste management (Kalbacher et al, 2012), geotechnical
applications (Xu et al, 2013; Böttcher et al, 2017; Nagel et al, 2017), geothermal energy systems (Hein et al,
2016) and energy storage (Nagel et al, 2016; Wang et al, 2017; Lehmann et al, 2017).

The most recent version, OpenGeoSys-6 (OGS-6)
(Naumov et al, 2018, Bilke et al 2019), is a fundamental re-implementation of the multi-physics code OpenGeoSys-4/5 (Kolditz
and Bauer, 2004; Wang and Kolditz, 2006) using advanced methods in software engineering and architecture
with a focus on code quality, modularity, performance and comprehensive documentation.

<!--
Particular emphasis is put on 
the implementation of advanced numerical methods for the propagation of discontinuities, such as enriched
finite element function spaces (Watanabe et al, 2012), non-local formulations (Parisio et al, 2018) and phase-field
models for fracture (Yoshioka et al, 2018) with the ability to utilize HPC platforms (Wang et al, 2015, 2017).
-->

To simplify the implementation of new constitutive models for solid
phases, `MFront` has been integrated through the `C` bindings of `MGIS`.

!["Elasto-plastic modelling of a cyclically loaded cavity in a cohesive-frictional material."](img/MCAS_disc_hole_cyclic_show_axes.png "Elasto-plastic modelling of a cyclically loaded cavity in a cohesive-frictional material."){#fig:OpenGeoSys}

Figure @fig:OpenGeoSys shows the results of a a cyclically loaded cavity
in a cohesive-frictional material modelled by a non-associated plastic
behaviour based on the Mohr-Coulomb criterion.


<!--
Current tests include elastic (isotropic and anisotropic),
elasto-plastic (see Fig. 2) and visco-plastic materials.
-->



## `JuliaFEM`

# Conclusions


# Acknowledgements

This research was conducted in the framework of the `PLEIADES`
project,which is supported financially by the CEA (Commissariat à
l’Energie Atomique et aux Energies Alternatives), EDF (Electricité de
France) and Framatome.Acknowledgements


This project uses code extracted from the following projects:

- https://github.com/bitwizeshift/string_view-standalone by Matthew
  Rodusek
- https://github.com/mpark/variant: by Michael Park
- https://github.com/progschj/ThreadPool by Jakob Progsch and Václav
  Zeman
- https://github.com/martinmoene/span-lite by Martin Moene
- https://bitbucket.org/fenics-apps/fenics-solid-mechanics/ by
  Kristian B. Ølgaard and Garth N. Wells.

# References
