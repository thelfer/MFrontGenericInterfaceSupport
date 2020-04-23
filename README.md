# MFrontGenericInterfaceSupport

This project aims at providing tools (functions, classes, bindings,
etc...) to handle behaviours written using `MFront` generic interface.
For information about `MFront`, see <http://tfel.sourceforge.net>.

Those tools are meant to be used by (`FEM`, `FFT`, etc.) solver
developers. This tools are *not* linked to the `TFEL` libraries.
Permissive licences have been chosen to allow integration in open-source
and proprietary codes.

This project is described in this paper:
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02003/status.svg)](https://doi.org/10.21105/joss.02003)

The official website can be found here:
<https://thelfer.github.io/mgis/web/index.html>.

## The `MFrontGenericInterface` `C++` library

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

## Bindings

### Existing 

The following bindings are available:

- `c` bindings.
- `python` bindings .
- `fortran03` bindings.
- `fenics` bindings (under current work in the `master` branch). Those
  bindings are strongly inspired by the `fenics-solid-mechanics`
  project. Those bindings are currently quite limited as mostly serve
  as a proof of concept. Note that `MGIS` can also be used in `FEniCS`
  through the `python` interface. This is discussed here:
  <https://thelfer.github.io/mgis/web/FEniCSBindings.html>.
- `julia` bindings (experimental)

### Future bindings (contributors are welcomed)

The following bindings are under consideration:

- `octave` binding

# Versions, branches, tags

- the `master` branch follows the evolution of the `master` branch of
  the `TFEL` project
- the `rliv-1.1` follows the evolution of the 3.3.x series of the `TFEL`
  project.
- the `rliv-1.0` follows the evolution of the 3.2.x series of the `TFEL`
  project. Note that this branch is **not** compatible with
  `TFEL-3.2.0`.

- tag `MFrontGenericInterfaceSupport-1.1` is meant to be build against
  `TFEL` 3.3.0
- tag `MFrontGenericInterfaceSupport-1.0` is meant to be build against
  `TFEL` 3.2.0
- tag `MFrontGenericInterfaceSupport-1.0.1` is meant to be build
  against `TFEL` 3.2.1

# Acknowledgement

This project uses code extracted from the following projects:

- https://github.com/bitwizeshift/string_view-standalone by Matthew
  Rodusek
- https://github.com/mpark/variant: by Michael Park
- https://github.com/progschj/ThreadPool by Jakob Progsch and Václav
  Zeman
- https://github.com/martinmoene/span-lite by Martin Moene
- https://bitbucket.org/fenics-apps/fenics-solid-mechanics/ by
  Kristian B. Ølgaard and Garth N. Wells.
