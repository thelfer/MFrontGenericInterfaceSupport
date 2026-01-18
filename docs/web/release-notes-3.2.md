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

# Acknowledgements

The authors are grateful to the many contributors to the `TFEL/MFront`
project. This research was conducted in the framework of the PLEIADES
project, which was supported financially by the CEA (Commissariat à
l’Énergie Atomique et aux Énergies Alternatives), EDF (Électricité de
France) and Framatome. Work on `MGIS/Function` was performed as part of
the EURATOM OperaHPC Project co-funded by the European Union.

# Issues fixed
