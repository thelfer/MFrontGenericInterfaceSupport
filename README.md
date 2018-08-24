# MFrontGenericInterfaceSupport

This project aims at proving tools (functions, classes, bindings,
etc...) to handle behaviours written using `MFront` generic interface.

Those tools are meant to be used by (`FEM`, `FFT`, etc.) solver
developers. This tools are *not* linked to the `TFEL` libraries.
Permissive licences have been chosen to allow integration in open-source
and proprietary codes.

## The `MFrontGenericInterface` `C++` library

The project is build around the `MFrontGenericInterface` library. This
library provides two main functions:

- load `MFront` behaviours from external shared libraries and retrieve
  meta data about the behaviour through the `mgis::behaviour::load`
  function. All the relevant information about the behaviour are stored
  in the `mgis::behaviour::Description` class.
- provide generic function to handle material data (such as internal
  state variables) at integration points and call `MFront` behaviours.

## Bindings

The following bindings are under consideration:

- `python` binding
- `c` binding
- `fortran90` binding
- `octave` binding

# Acknowledgement

This project uses code extracted from the following projects:

- https://github.com/bitwizeshift/string_view-standalone by Matthew
  Rodusek
- https://github.com/mpark/variant: by Michael Park
- https://github.com/progschj/ThreadPool by Jakob Progsch and Václav
  Zeman
- https://github.com/martinmoene/span-lite by Martin Moene