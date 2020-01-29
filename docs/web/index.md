% The MFrontGenericInterfaceSupport project
% Thomas Helfer
% 20/11/2018

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

# Bindings

## Existing bindings

The following bindings are available:

- `c` binding
- `python` binding 
- `fortran` binding
- `julia` binding
- `fenics` bindings (experimental). Those bindings are strongly inspired
  by the `fenics-solid-mechanics` project. Those bindings are currently
  quite limited as mostly serve as a proof of concept. Note that `MGIS`
  can also be used in `FEniCS` through the `python` interface. This is
  discussed [here](FEniCSBindings.html).

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

## Release notes

- [Version 1.1](release-notes-1.1.html)

# Installation

The installation process is discussed [here](install.html).

# Acknowledgement

This project uses code extracted from the following projects:

- <https://github.com/bitwizeshift/string_view-standalone> by Matthew
  Rodusek
- <https://github.com/mpark/variant>: by Michael Park
- <https://github.com/progschj/ThreadPool> by Jakob Progsch and Václav
  Zeman
- <https://github.com/martinmoene/span-lite> by Martin Moene
