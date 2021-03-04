---
title: MFrontGenericInterfaceSupport Version 1.2 
author: Thomas Helfer
date: 2020
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

Version 1.2 of `MFrontGenericInterfaceSupport` is compatible with the
Version 3.4 of `TFEL/MFront`.

# Known incompatibilities

The `stored_energies` and `dissipated_energies` are now automatically
allocated only of the behaviours by the `MaterialStateManager` class
only if the behaviours is able to compute them.

# New functionalities

## Orthotropic behaviours

For orthotropic behaviours, the `Behaviour` structure exposes \(6\)
function pointers:

- `rotate_gradients_ptr`: pointer to a function implementing the
  rotation of the gradients from the global frame to the material frame.
- `rotate_array_of_gradients_ptr`: pointer to a function implementing
  the rotation of an array of gradients from the global frame to the
  material frame.
- `rotate_thermodynamic_forces_ptr`: pointer to a function implementing
  the rotation of the thermodynamic forces from the material frame to
  the global frame.
- `rotate_array_of_thermodynamic_forces_ptr`: pointer to a function
  implementing the rotation of an array of thermodynamic forces from the
  material frame to the global frame
- `rotate_tangent_operator_blocks_ptr`: pointer to a function
  implementing the rotation of the tangent operator blocks from the
  material frame to the global frame.
- `rotate_array_of_tangent_operator_blocks_ptr`: pointer to a function
  implementing the rotation of an array of tangent operator blocks from
  the material frame to the global frame.

Those functions takes pointer to the raw memory. The callee is
responsible of the consistency of the data.

> **In place transformations**
>
> All those functions take two parameters: the pointer to the rotated
>
> data on output and the pointer to the original data on input. In place
> transformations is allowed, i.e. those pointers can be equal.

> **The rotation matrix argument**
>
> All those functions takes the rotation matrix from the global frame to
> the material frame as last argument. If required, i.e. for
> thermodynamic forces and tangent operator blocks, this matrix is
> transposed internally to have the inverse transformation.
>
> The rotation matrix is given as a \(3\times3\) matrix, packed in
> an \(9\) continuous array in `C`-like column-major storage.
>
> No checks are made to ensure that the columns of the matrix makes
> and orthonormal basis of \(\mathcal{R}^{3}\). In \(1D\), this matrix
> is discarded an no operation is performed. In \(2D\), only the
> upper-left part of the matrix is used.

For convenience (and debugging), the call to those functions pointers
are mapped into the following free functions: `rotateGradients`,
`rotateArrayOfGradients`, `rotateThermodynamicForces` and
`rotateArrayOfThermodynamicForces`. Those functions perform additional
consistency checks (compared to the functions exposed by the `Behaviour`
class) which might hurt performances, especially when dealing with one
integration point only. Each of these functions is overloaded twice for
in-place operations and out-of-place operations.

### Example

The following example shows how to rotate the gradients of a small
strain strain behaviour in generalised plane strain:

~~~~{.cxx}
const std::array<real, 9> r = {0, 1, 0, 1, 0, 0, 0, 0, 1};
const std::array<real, 4> ge = {1, 0, 0, 0};
std::array<real, 4> me;
rotateGradients(me, b, ge, r);
~~~~

## Update to the `C` bindings

The following functions are now available:

- `mgis_bv_get_space_dimension`: this functions returns the space
  dimension associated with an hypothesis.
- `mgis_bv_get_stensor_size`: this functions returns the number of
  components of a symmetric tensor for the given hypothesis.
- `mgis_bv_get_tensor_size`: this function returns the number of
  components of a tensor for the given hypothesis.
- `mgis_bv_get_variable_size`: this function returns the size of a
  variable (i.e. the number of components) for the given hypothesis.

## Fortran bindings

The following functions are now available in the `mgis_behaviour`
module:

- `get_space_dimension`: this functions returns the space
  dimension associated with an hypothesis.
- `get_stensor_size`: this functions returns the number of
  components of a symmetric tensor for the given hypothesis.
- `get_tensor_size`: this function returns the number of
  components of a tensor for the given hypothesis.
- `get_variable_size`: this function returns the size of a
  variable (i.e. the number of components) for the given hypothesis.

# Issues solved

## Issues #54: Inform the calling code about `@DissipatedEnergy` and/or `@InternalEnergy`

The `Behaviour` class now exposes two new boolean data members:

- `computesStoredEnergy`: if true, the behaviour computes the stored energy
- `computesDissipatedEnergy`: if false, the behaviour computes the dissipated energy

In the `C` bindings, the `mgis_bv_behaviour_computes_stored_energy`  and `mgis_bv_behaviour_computes_dissipated_energy` functions are now available.

In the `fortran` bindings, the functions `behaviour_computes_stored_energy` and `behaviour_computes_dissipated_energy` are now available in the `mgis_behaviour` module.

In the `python` bindings, the `Behaviour` class now exposes two read only properties: `computesStoredEnergy` and `computesDissipatedEnergy`.

The `MaterialDataManager` constructor now only allocates the memory associated with the stored and disspated energies only if the behaviour computes those energies.

For details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/54>


## Issue #33: Function for checking if the behaviour is a Finite Strain one

The `mgis::behaviour::isStandardFiniteStrainBehaviour` has been added to
check if a behaviour is a finite strain behaviour and if its kinematic
is also standard (i.e. is of the `F-Cauchy` kind although the stress
measure can be chosen when loading the behaviour).

This function is exposed as:

- `mgis_bv_is_standard_finite_strain_behaviour` in the `C`' bindings.
- `is_standard_finite_strain_behaviour` in the `mgis_behaviour` module
  in the `Fortran`' bindings.
- `isStandardFiniteStrainBehaviour` in in the `mgis.behaviour` module in
  the `Python`' bindings.

For details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/33>

## Issue #32: Better `const` correctness

This issue follows this evolution in the generic interface of MFront:
<https://sourceforge.net/p/tfel/tickets/212/>

The state at the beginning of the time step is now described in a
structure called `mgis_bv_InitialStateView`, the fields of which are all
const.

The following fields of the `mgis_bv_StateView` are now `const`:

- `gradients`
- `material_properties`
- `external_state_variables`

For details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/32>

## Issue #30: Variables array size depending on type and modelling hypothesis

The `get_variable_size` function is now available in the
`mgis_behaviour` module. This function returns the size of a variable
(i.e. the number of components) for the given hypothesis.

For details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/30>

## Issue #28: Missing block jacobian info in `python` bindings

The tangent operator blocks are now available in the `python` bindings
under the `tangent_operator_blocks` property of the `Behaviour` class.

Those blocks are accessible as an array of tuples of instances of the
`Variable` class.

This feature can be used as follows:

~~~~{.python}
import mgis.behaviour as mgis_bv

h = mgis_bv.Hypothesis.Tridimensional
b = mgis_bv.load('src/libBehaviour.so','StationaryHeatTransfer2',h)
for t in b.tangent_operator_blocks:
    print('d{}_d{}'.format(t[0].name,t[1].name))
~~~~

For details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/28>

## Issue #27: Get internal state variables array size via fortran bindings

For details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/27>

## Issue #25: Exporting MGIS build tree for use by external projects

The following files are now generated and exported:

- `MFrontGenericInterfaceConfig.cmake`: configuration file for the
  `MFrontGenericInterface` library which contains the `C++` core library
  of `MGIS`.
- `MFrontGenericInterface-cConfig.cmake`: configuration file for the
  `MFrontGenericInterface` library which contains the `c` binding of
  `MGIS`.
- `MFrontGenericInterface-fortranConfig.cmake`: configuration file for
  the `MFrontGenericInterface-fortran` library which contains the
  `fortran` bindings of `MGIS`.

Those files can be used as follows:

~~~~{.cmake}
find_package (MFrontGenericInterface REQUIRED)
~~~~

The previous instruction imports the `mgis::MFrontGenericInterface` target, which can be used as follows:

~~~~{.cmake}
  target_link_libraries(HybridHighOrder
	PRIVATE mgis::MFrontGenericInterface)
~~~~

For details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/25>

# References