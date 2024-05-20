---
title: The MFrontGenericInterfaceSupport project
author: Thomas Helfer
date: 2020
lang: en-EN
numbersections: true
link-citations: true
colorlinks: true
figPrefixTemplate: "$$i$$"
tabPrefixTemplate: "$$i$$"
secPrefixTemplate: "$$i$$"
eqnPrefixTemplate: "($$i$$)"
bibliography: bibliography.bib
---

<!-- <div id="slideshow"> -->
<!--   <ul class="slides"> -->
<!--     <li> -->
<!--       <iframe width="560" height="315" src="https://www.youtube.com/embed/juWMIkJ64iE" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"> -->
<!--       </iframe> -->
<!--     </li> -->
<!--   </ul> -->
<!--   <span class="arrow previous"></span> -->
<!--   <span class="arrow next"></span> -->
<!-- </div> -->
<!-- <script src="http://ajax.googleapis.com/ajax/libs/jquery/1.4.2/jquery.min.js"></script> -->
<!-- <script src="js/slideshow.js"></script> -->

![Principle of MGIS](img/mgis.svg){width=100%}

This project aims at proving tools (functions, classes, bindings,
etc...) to handle behaviours written using `MFront` generic interface.
For information about `MFront`, see <http://tfel.sourceforge.net>.

Those tools are meant to be used by (`FEM`, `FFT`, etc.) solver
developers. This tools are *not* linked to the `TFEL` libraries.
Permissive licences have been chosen to allow integration in open-source
and proprietary codes.

# The `MFrontGenericInterface` `C++` library

The project is build around the `MFrontGenericInterface` library. This
library provides two main functions:

- the `mgis::behaviour::load` functions loads `MFront` behaviours from
  external shared libraries and retrieve all relevant meta data
  function. Those relevant information are stored in the
  `mgis::behaviour::Behaviour` class.
- the `mgis::behaviour::integrate` integrates the behaviour over one
  time step. The data associated with an integration point are handled
  by the `mgis::behaviour::BehaviourData` class which contains the state
  of the integration point at the beginning and at the end of the time
  step.

The library also supports handling a group of integration points though
the `mgis::behaviour::MaterialStateManager` class.

An introduction to the `C++` library may be found [here](bindings-cxx.html)

# Projects using `MGIS` or in which `MGIS` can be coupled

- [`code-aster`](https://code-aster.org)
- [`OpenGeoSys`](https://www.opengeosys.org/). See this
  [talk](https://github.com/thelfer/tfel-doc/blob/master/MFrontUserDays/FifthUserDays/OpenGeoSys-Nagel-MFrontUserDays-2019.pdf)
  for details.
- [`FEniCS`](https://fenicsproject.org/). Coupling with `FEniCS` can be
  made using the `python` bindings (See this
  [talk](https://github.com/thelfer/tfel-doc/blob/master/MFrontUserDays/FifthUserDays/OpenGeoSys-Nagel-MFrontUserDays-2019.pdf)
  and this
  [tutorial](https://github.com/thelfer/tfel-doc/blob/master/MFrontUserDays/FifthUserDays/OpenGeoSys-Nagel-MFrontUserDays-2019.pdf)
  for details). However, we recommend using the higher level
  [`mgis.fenics`](mgis_fenics.html) module.
- [`MFEM`](https://mfem.org/) through the
  [`MFEM-MGIS`](https://github.com/thelfer/mfem-mgis/). See [this
  talk](https://github.com/thelfer/tfel-doc/blob/master/MFrontUserDays/SeventhUserDays/mfem-mgis.pdf)
  for details.
- [`MoFEM`](http://mofem.eng.gla.ac.uk/mofem/html/). See [this
  talk](https://github.com/thelfer/tfel-doc/blob/master/MFrontUserDays/SixthUserDays/Talk6-MoFEM_MFront.pdf)
  for details.
- [`Europlexus`](https://europlexus.jrc.ec.europa.eu/)
- `XPer`. See this
  [talk](https://github.com/thelfer/tfel-doc/blob/master/MFrontUserDays/FifthUserDays/Xper-Perales-MFrontUserDays-2019.pdf)
  for details.
- [`MBDyn`](https://www.mbdyn.org/). See [this demo](https://www.youtube.com/watch?v=I8HENx5mszA).
- The (yet to be released) solver `MANTA` developed at CEA [@jamond_manta_2022].
- [`JuliaFEM`](http://www.juliafem.org/) through the
  [`MFrontInterface.jl`](https://github.com/JuliaFEM/MFrontInterface.jl)
  package.
- [`Kratos Multiphysics`](https://github.com/KratosMultiphysics/Kratos)
  through the experimental `MGISApplication`
- An example using `MGIS` and the [`DUNE`
  platform](https://www.dune-project.org/) is available
  [here](https://github.com/thelfer/dune-mgis).

# Bindings

## Existing bindings

The following bindings are available:

- `c` binding
- `python` binding: this binding also provide an interface to `FEniCS` which is 
  discussed [here along with a collection of commented demos](mgis_fenics.html)
- `fortran` binding
- `julia` binding

## Future bindings (contributors are welcomed)

The following bindings are under consideration:

- `octave` binding

# Versions, branches

## Branches

- the `master` branch follows the evolution of the `master` branch of
  the `TFEL` project
- the `rliv-1.1` follows the evolution of the 3.3.x series of the `TFEL`
  project.
- the `rliv-1.0` follows the evolution of the 3.2.x series of the `TFEL`
  project. Note that this branch is **not** compatible with
  `TFEL-3.2.0`.

# Acknowledgement

This project uses code extracted from the following projects:

- <https://github.com/bitwizeshift/string_view-standalone> by Matthew
  Rodusek
- <https://github.com/mpark/variant>: by Michael Park
- <https://github.com/progschj/ThreadPool> by Jakob Progsch and Václav
  Zeman
- <https://github.com/martinmoene/span-lite> by Martin Moene


