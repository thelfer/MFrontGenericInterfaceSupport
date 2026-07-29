---
title: MFrontGenericInterfaceSupport Version 3.1
author: Thomas Helfer
date: 2025
lang: en-EN
numbersections: true
documentclass: article
from: markdown+tex_math_single_backslash
geometry:
  - margin=2cm
papersize: a4
link-citations: true
colorlinks: true
figPrefixTemplate: "$$i$$"
tabPrefixTemplate: "$$i$$"
secPrefixTemplate: "$$i$$"
eqnPrefixTemplate: "($$i$$)"
bibliography: bibliography.bib
---

This version is meant to be used with `TFEL` Version 5.1. This version
is the first that can optionnaly use directly the `TFEL` libraries.

This version was released on December 15, 2025. This version was
released on December, 15 2025 along with:

- [Version 1.1.2](release-notes-1.1.2.html)
- [Version 1.2.3](release-notes-1.2.3.html)
- [Version 2.0.1](release-notes-2.0.1.html)
- [Version 2.1.1](release-notes-2.1.1.html)
- [Version 2.2.2](release-notes-2.2.2.html)
- [Version 3.0.2](release-notes-3.0.2.html)

and inherits from all fixes from those releases.

# The `MGIS/Function` library

The main new feature of `MGIS` is the `MGIS/Function`[^functions]
library, which has been introduced to provide standard mechanical
post-processings based on the `TFEL/Math` and `TFEL/Material` libraries.
`MGIS/Function` allows reusing the large amount of documented and
verified functionalities provided by those libraries at the structural
scale, including: rotation from/to material frame, conversion between
stress measures, Hencky and Green-Lagrange strains, Tresca, von Mises,
Hosford, Hill, Barlat equivalent stresses, principal values, Lode's
angle, etc.

[^functions] A function denote a set of values associated with the
elements (nodes, quadrature points, cells) of a discretized space.
Functions are often called `fields` in many solvers.

The following snippet computes the Cauchy stress in the global frame
from the first Piola-Kirchhoff stress known in the material frame.

~~~~{.cxx}
const auto evaluator = pk1_function | as_tensor<3>                    |
                       from_pk1_to_cauchy (F_function | as_tensor<3>) |
                       rotate_backwards (R_function | as_matrix<3,3>);
const auto ok = assign(ctx, sig | as_stensor<3>, evaluator);
~~~~

The library is developed in `C++-20` with high-quality requirements and
strives to be `constexpr` friendly, memory-safe, thread-safe,
exception-safe, and `GPU-friendly`.

The data structures of the functions depend on the targeted solver, as
well as the implementation of the `assign`, which depends on the
programming model chosen by the targeted solver (`OpenMP`, `CUDA`,
`SYCL`, `Kokkos`, etc.).

For standard data structures, such as the one used by `Numpy`'s array,
so-called `views` are provided allowing direct usage of the library. In
this case, an implementation of the `assign` algorithm is provided by
the parallel version of the Standard library. Depending on the compiler
and implementation of the standard library, various parallel programming
models are used (multithreading and GPU-offloading).

The `MGIS/Function` library is described on [this page](functions.html).

# Improvements

## The `integrate_debug` functions

The `integrate_debug` function can be used to generate debug files to
analyse integration failures. This function is very customizable, see
[this page](behaviour-integration-failure-analysis.html) for details.

### Example of usage

This snippet shows that the `integrate_debug` functions can be used as
drop-in replacement for the `integrate` function.

~~~~{.cxx}
const auto r = integrate_debug(v, b);
~~~~

## `getDatabase` and `loadFromDatabase`

> **Note**
>
> This feature requires `MGIS` to be compiled againt the `TFEL` libraries.

The `mgis::getDatabase` function returns a global instance of an
`MFrontDatabase` object. See [this
page](https://thelfer.github.iotfel/web/tfel-mfront-database.html) for
more details. This instance can be populated using a library, a
directory or a list of directories specified in an environment variable.

The function `mgis::behaviour::loadFromDatabase` relies on the previous
database to load a behaviour matching the given criteria:

- the name (required) of the behaviour of the loaded,
- the material (optional)

### Example of usage

The following example registers all the material knownledge available in
the shared libraries found in directories listed in an environment
variable and lists all the behaviours compiled with the `generic`
interface:

~~~~{.python}
import mgis
db = mgis.getDatabase()
db.analyseDirectoriesListedInEnvironmentVariable('MFRONT_GALLERY_LIBRARIES',
                                                 ignore_errors = True)
for e in db.getEntryPoints(type = 'behaviour',
                           interface = 'generic',
                           material = 'concrete'):
    print(f"- {e.library}: {e.name}")
~~~~

The following snippet populates the database using all the libraries of
a given directory and loads a behaviour whose name of which contains
`Elasticity` and which is associated with the material wood.

~~~~{.python}
import mgis
import mgis.behaviour as mgis_bv
db = mgis.getDatabase()
db.analyseDirectory('/home/th202608/codes/MFrontGallery/master/install/lib/',
                    ignore_errors = True)
b = mgis_bv.loadFromDatabase(name='\\w+Elasticity\\w+', material='Wood',
                             hypothesis = mgis.behaviour.Hypothesis.Tridimensional)
~~~~

# Python bindings

Python bindings are now generated using the
[`pybind11`](https://github.com/pybind/pybind11) library.

## New tutorials

Two tutorials show how to use the python bindings to interact with a
finite-strain mechanical behaviour compiled with `MFront` and develop an
alternative to `MTest`:

- [case of a finite strain behaviour](bindings-python-finite_strain.html)
- [case of a small strain behaviour](bindings-python-small_strain.html)

# Tested configurations

+---------------------------------------+----------------------------+------------+----------------+----------------------+-----------------------+
|     System                            | Compiler                   | Build type | Options        | Generator            | Percentage of success |
+:=====================================:+:==========================:+:==========:+:==============:+:====================:+:=====================:+
| Ubuntu 22.04.5 LTS  (Jammy Jellyfish) | gcc 11.4.0                 | Release    |                | Unix Makefiles       | 100% (76/76)          |
+---------------------------------------+----------------------------+------------+----------------+----------------------+-----------------------+
| Ubuntu 22.04.5 LTS  (Jammy Jellyfish) | gcc 11.4.0                 | Debug      |                | Unix Makefiles       | 100% (76/76)          |
+---------------------------------------+----------------------------+------------+----------------+----------------------+-----------------------+
| Ubuntu 22.04.5 LTS  (Jammy Jellyfish) | gcc 11.4.0                 | Release    | fast-math      | Unix Makefiles       | 100% (76/76)          |
+---------------------------------------+----------------------------+------------+----------------+----------------------+-----------------------+
| Ubuntu 22.04.5 LTS  (Jammy Jellyfish) | gcc 11.4.0                 | Release    |                | Ninja                | 100% (76/76)          |
+---------------------------------------+----------------------------+------------+----------------+----------------------+-----------------------+
| Ubuntu 22.04.5 LTS  (Jammy Jellyfish) | gcc 14.2.0                 | Release    |                | Unix Makefiles       | 100% (76/76)          |
+---------------------------------------+----------------------------+------------+----------------+----------------------+-----------------------+
| Ubuntu 22.04.5 LTS  (Jammy Jellyfish) | gcc 15.1.0                 | Release    |                | Unix Makefiles       | 100% (76/76)          |
+---------------------------------------+----------------------------+------------+----------------+----------------------+-----------------------+
| Ubuntu 22.04.5 LTS  (Jammy Jellyfish) | llvm 14.0.0                | Release    |                | Unix Makefiles       | 100% (76/76)          |
+---------------------------------------+----------------------------+------------+----------------+----------------------+-----------------------+
| Ubuntu 22.04.5 LTS  (Jammy Jellyfish) | llvm 19.1.0                | Release    |                | Unix Makefiles       | 100% (76/76)          |
+---------------------------------------+----------------------------+------------+----------------+----------------------+-----------------------+
| Ubuntu 22.04.5 LTS  (Jammy Jellyfish) | llvm 21.0.0                | Release    |                | Unix Makefiles       | 100% (32/32)          |
+---------------------------------------+----------------------------+------------+----------------+----------------------+-----------------------+
| Ubuntu 22.04.5 LTS  (Jammy Jellyfish) | llvm 21.0.0                | Release    | libc++         | Unix Makefiles       | 100% (32/32)          |
+---------------------------------------+----------------------------+------------+----------------+----------------------+-----------------------+
| Ubuntu 22.04.5 LTS  (Jammy Jellyfish) | Intel 2024.2               | Release    |                | Unix Makefiles       |  98% (61/62)          |
+---------------------------------------+----------------------------+------------+----------------+----------------------+-----------------------+
| Ubuntu 22.04.5 LTS  (Jammy Jellyfish) | NVHPC 25.3                 | Release    |                | Unix Makefiles       |  98% (61/62)          |
+---------------------------------------+----------------------------+------------+----------------+----------------------+-----------------------+
| Windows 11                            | Visual Studio 17.4         | Release    |                |                      | 100% (46/46)          |
+---------------------------------------+----------------------------+------------+----------------+----------------------+-----------------------+

## Github actions

Github actions allows to check compilation on Ubuntu (latest) `gcc` 13.3.0.

## Continuous integration provide by the `PLEIADES` project

+---------------------------------------+---------------+------------+-----------+----------------------+-----------------------+
|     System                            | Compiler      | Build type | Options   | Generator            | Percentage of success |
+:=====================================:+:=============:+:==========:+:=========:+:====================:+:=====================:+
| Ubuntu 24.04   LTS  (Noble Numbat)    | gcc  13.3.0   | Release    |           | Unix Makefiles       | 100% (62/62)          |
+---------------------------------------+---------------+------------+-----------+----------------------+-----------------------+
| Ubuntu 22.04.5 LTS  (Jammy Jellyfish) | gcc  11.4.0   | Release    |           | Unix Makefiles       | 100% (62/62)          |
+---------------------------------------+---------------+------------+-----------+----------------------+-----------------------+
| Debian 12.12 (bookworm)               | gcc  12.2.0   | Release    |           | Unix Makefiles       | 100% (62/62)          |
+---------------------------------------+---------------+------------+-----------+----------------------+-----------------------+
| Debian 11.11 (bullseye)               | gcc  11.3.0   | Release    |           | Unix Makefiles       | 100% (62/62)          |
+---------------------------------------+---------------+------------+-----------+----------------------+-----------------------+
| Debian 10.13 (buster)                 | gcc  11.3.0   | Release    |           | Unix Makefiles       | 100% (62/62)          |
+---------------------------------------+---------------+------------+-----------+----------------------+-----------------------+
| Debian  9.13 (stretch)                | gcc  11.3.0   | Release    |           | Unix Makefiles       | 100% (62/62)          |
+---------------------------------------+---------------+------------+-----------+----------------------+-----------------------+

# Acknowledgements

The authors are grateful to the many contributors to the `TFEL/MFront`
project. This research was conducted in the framework of the PLEIADES
project, which was supported financially by the CEA (Commissariat à
l’Énergie Atomique et aux Énergies Alternatives), EDF (Électricité de
France) and Framatome. Work on `MGIS/Function` was performed as part of
the EURATOM OperaHPC Project co-funded by the European Union.

# Issues fixed

## Issue 148: [doc] better document `rdt` usage

For more details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/148>.