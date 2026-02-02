---
title: MFrontGenericInterfaceSupport Version 3.2
author: Thomas Helfer
date: 2026
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

This version is meant to be used with `TFEL` Version 5.2.

# Documentation

## New `cmake` options

### HDF5 support

The following options control the support of the `HDF5` library:

- `enable-hdf5-support`: enable HDF5 support for save/restore
  operations. Note that if this support is not explicitely requested,
  `MGIS` still tries by default to support `HDF5` but configuration will
  not fail if the `HDF5` library is not found. See
  `enable-hdf5-automatic-support` for details.
- `enable-hdf5-automatic-support`: if set, `MGIS` tries by default to
  support `HDF5` even if enable-hdf5-support` is not set.

## `Doxygen` documentation

The doxygen documentation is now online:
<https://thelfer.github.io/mgis/doxygen/index.html>.

# New features

## Scripts to define environment variables for `MGIS` to work properly

Depending on the system and compilation options, some of following
variables shall be set for `MGIS` to work properly: `MGISHOME`, `PATH`,
`LD_LIBRARY_PATH` and `PYTHONPATH`.

`MGIS` now installs automatically the following files in the installation
directory (refered to `<install_prefix>` in the following):

- `<install_prefix>/share/mgis/env/env.sh` for `UNIX` systems and the
  `bash` shell. This file shall be used as follows:

  ~~~~{.sh}
  $ source <install_prefix>/share/mgis/env/env.sh
  ~~~~
- `<install_prefix>\share\mgis\env\env.ps1` for `PowerShell`
  shell under `Windows`. This file shall be used as follows:

  ~~~~{.sh}
  $ .\<install_prefix>\share\mgis\env\env.ps1
  ~~~~
- `<install_prefix>\share\mgis\env\env.bat` for the historical `cmd`
  shell under `Windows`. This file shall be used as follows:

  ~~~~{.sh}
  $ call <install_prefix>\share\mgis\env\env.bat
  ~~~~

> **Note**
>
> Those variables are not required if `MGIS` is installed
> system-wide (for instance in `/usr/local`) and that the `MGIS`'s
> binaries are not relocated (i.e. moved to a different directory than
> the one specified during the compilation process as the installation
> directory).

> **Note**
>
> If `MGIS` is built with `TFEL` support, the `TFEL` environment
> shall be properly set.

# New features of the `MGIS/Behaviour` library

## Save/restore operations in `MaterialStateManager` and `MaterialDataManager`

If `HDF5` support is enabled, the functions `save` and `restore` allows
saving and restoring the values stored in a `MaterialStateManager` or in
a `MaterialDataManager` to an `HDF5` file.

## Control on variables updated/reverted in `MaterialStateManager`

By default, material properties, mass density and external state
variables are updated by the `update` and `revert` functions.

This can now be controlled by passing a value of the
`MaterialStateManager::UpdatePolicy` type to the functions
`setMaterialProperty`, `setMassDensity` or `setExternalStateVariable`.

### Example of usage

~~~~{.cxx}
setExternalStateVariable(m.s1, "Temperature", T1,
                         MaterialStateManager::UPDATE);
~~~~

# New features of the `MGIS/Function` library

## Functions using a strided memory access

The following classes have been introduced:

- `StridedCoalescedMemoryAccessTensorView`
- `StridedCoalescedMemoryAccessCompositeTensorsView`

### `StridedCoalescedMemoryAccessTensorView`

`StridedCoalescedMemoryAccessTensorView` is a tensorial function view
which stores its components in non interleaved manner using the
following scheme:

~~~~
| <------- Component 1 ---------> |....| <----- Component Nc ---------> |
+-------++-------++------++-------+----+-------++------++------++-------+
| Elt 1 || Elt 2 || .... || Elt N |....| Elt 1 ||Elt 2 || .... || Elt N |
+-------++-------++------++-------+----+-------++------++------++-------+
~~~~

#### Example of usage

~~~~{.cxx}
constexpr auto ne = size_type{2};
auto space = BasicLinearSpace{ne};
std::array<const real, 4 * ne> values = {1, 10, 2, 20, 3, 30, 4, 40};
const auto f = StridedCoalescedMemoryAccessTensorView<
    BasicLinearSpace, tfel::math::stensor<2, real>, false>{space, values};
const auto e1 = f(0);
// e1 = {1, 2, 3, 4}
const auto e2 = f(1);
// e2 = {10, 20, 30, 40}
~~~~
  
## `StridedCoalescedMemoryAccessCompositeTensorsView`

`StridedCoalescedMemoryAccessCompositeTensorsView` allows retrieving
scalar or tensorial objects which are stored in a non interleaved
manner.

#### Example of usage

~~~~{.cxx}
using CompositeFunctionView = 
    StridedCoalescedMemoryAccessCompositeTensorsView<BasicLinearSpace, 4,
                                                     false>;
constexpr auto ne = size_type{2};
auto space = BasicLinearSpace{ne};
std::array<const real, 4 * ne> values = {1, 10, 2, 20, 3, 30, 4, 40};
const auto f = CompositeFunctionView{space, values};
const auto e1 = f.get<0, tfel::math::stensor<2, real>>(0);
// e1 = {1, 2, 3, 4}
const auto e2 = f.get<0, tfel::math::stensor<2, real>>(1);
// e2 = {10, 20, 30, 40}
~~~~

# Issues fixed

## Issue 210: [documentation] Deployment of the doxygen documentation

For more details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/210>

## Issue 209: Add basic support for `save`/`restore` operations in  `MaterialDataManager`

For more details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/209>

## Issue 201: Add the ability to update only the state variables of a `MaterialStateManager`

For more details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/201>

## Issue 200: Create environment file in the installation directory

For more details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/200>

## Issue 196: [MGIS/Function] Add function view with strided coalesent memory access

For more details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/196>￼

# Acknowledgements

The authors are grateful to the many contributors to the `TFEL/MFront`
project. This research was conducted in the framework of the PLEIADES
project, which was supported financially by the CEA (Commissariat à
l’Énergie Atomique et aux Énergies Alternatives), EDF (Électricité de
France) and Framatome. Work on `MGIS/Function` was performed as part of
the EURATOM OperaHPC Project co-funded by the European Union.
