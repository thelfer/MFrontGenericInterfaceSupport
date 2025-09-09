---
title: MFrontGenericInterfaceSupport Version 2.1
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

# Known incompatibilities

## Integral values of `Variable::VECTOR` and `Variable::STENSOR`

The integral values of `Variable::VECTOR` and `Variable::STENSOR` have
changed to match `TFEL` conventions.

## Removed function in `C++`

The following functions were removed:

- `getUniformMaterialProperty`
- `getNonUniformMaterialProperty`
- `getUniformExternalStateVariable`
- `getNonUniformExternalStateVariable`

## `C` bindings

### Integral values of `MGIS_BV_VECTOR` and `MGIS_BV_STENSOR`

The integral values of `MGIS_BV_VECTOR` and `MGIS_BV_STENSOR` have
changed to match `TFEL` conventions.

### Removed functions

The following functions were removed:

- `mgis_bv_material_state_manager_get_uniform_material_property`
- `mgis_bv_material_state_manager_get_non_uniform_material_property`
- `mgis_bv_material_state_manager_get_uniform_external_state_variable`
- `mgis_bv_material_state_manager_get_non_uniform_external_state_variable`

## `fortran` bindings

### `mgis_behaviour` module

#### Integral values of `VECTOR` and `STENSOR` 

The integral values of `VECTOR` and `STENSOR` have changed to match
`TFEL` conventions.

#### Removed functions

The following functions were removed from the `mgis_behaviour` module:

- `material_state_manager_get_uniform_material_property`
- `material_state_manager_get_non_uniform_material_property`
- `material_state_manager_get_uniform_external_state_variable`
- `material_state_manager_get_non_uniform_external_state_variable`

# New features

## Support for extended types {#sec:mgis:2.1:extended_types}

In previous versions, only four types of variables were supported
(scalars, vectors, symmetric tensors and unsymmetric tensors). The size
of those variables is a function on the modelling hypothesis. The type
of the variables was given by the `Variable::Type` enumeration which
could hold the values `SCALAR`, `VECTOR`, `STENSOR`, and `TENSOR`.

The `Variable::Type` enumeration may now hold the following values:
`SCALAR`, `VECTOR`, `VECTOR_1D`, `VECTOR_2D`, `VECTOR_3D`, `STENSOR`,
`STENSOR_1D`, `STENSOR_2D`, `STENSOR_3D`, `TENSOR`, `TENSOR_1D`,
`TENSOR_2D`, `TENSOR_3D`, `HIGHER_ORDER_TENSOR` and `ARRAY`.

The `Variable` class exposes an integer named `type_identifier`
describing more precisely the type of the variable. The meaning of this
identifier is fully described in [this
page](https://thelfer.github.io/tfel/web/mfront-types.html). 

The `getVariableTypeSymbolicRepresentation` returns a symbolic
representation of a object using a `C++`-like representation from a type
identifier, as follows:

~~~~{.bash}
$ python3
>>> import mgis.behaviour
>>> print(mgis.behaviour.getVariableTypeSymbolicRepresentation(780))
derivative_type<stensor<N, real>, tensor<N, real>>
~~~~

## Support for behaviours' initialize functions {#sec:mgis:2.1:initialize_functions}

Since version 4.1, `MFront` behaviours can declare initialize functions
which are meant to initialize the state of the material.

The available initialize functions are described by the
`initialize_functions` member of the `Behaviour` class which associates
the name of the initialize function and a small data structure
containing:

- the pointer to the initialize function,
- the list of inputs of the initialize function, if any.

The `getInitializeFunctionVariablesArraySize` function returns the size
of an array able to contain the inputs an initialize function for an
integration point.

The `allocateInitializeFunctionVariables` functions return an array able
to store the inputs of an initialize function for a given integration
point or a set of integrations points.

The `executeInitializeFunction` functions execute an initialization
function on a unique integration point or a set of integration points.

> **Note about the result of the initialize function**
>
> The `BehaviourDataView` class imposes that the initial state is immutable.
> Thus, initialize functions can thus only initialize the state of the material
> at the end of the time step. In most case, the call to the selected
> initialize functions shall be follow to a call to the `update` function.

> **Note about the material properties and the external state variables**
>
> The material properties and the external state variables must be set
> before calling the `executeInitializeFunction` functions. Only the
> values at the beginning of the time step are passed to the behaviour.

### Example of usage

~~~~{.cxx}
auto d = BehaviourData{b};
// initialize the material properties and the external state variables
...
// calling an initialize function which requires an input
auto inputs = allocateInitializeFunctionVariables(b, "StressFromInitialPressure");
inputs[0] = pr;
auto v = make_view(d);
executeInitializeFunction(v, b, "StressFromInitialPressure", inputs);
~~~~

## Support for behaviours' post-processings {#sec:mgis:2.1:postprocessings}

Since version 4.1, `MFront` behaviours can declare postprocessings which
are meant to process the state of the material after the behaviour
integration.

The available postprocessings are described by the `postprocessings`
member of the `Behaviour` class which associates the name of the
postprocessing and a small data structure containing:

- the pointer to the postprocessing,
- the list of outputs of the postprocessing.

The `getPostProcessingVariablesArraySize` function returns the size
of an array able to contain the outputs an postprocessing for an
integration point.

The `allocatePostProcessingVariables` functions return an array able to
store the outputs of a postprocessing for a given integration point or a
set of integrations points.

The `executePostProcessing` functions execute a postprocessing on a
unique integration point or a set of integration points.

### Example of usage

~~~~{.cxx}
auto m = MaterialDataManager{b, 2u};
// initialize the state and perform the behaviour integration
...
// execute the post-processing
auto outputs = allocatePostProcessingVariables(m, "PrincipalStrain");
executePostProcessing(outputs, m, "PrincipalStrain");
~~~~

## Utility function to extract the value of an internal state variable {#sec:mgis:2.1:extractInternalStateVariable}

The `extractInternalStateVariable` function can now be used to extract
the value of an internal state variable in a pre-allocated buffer.

### Example of usage

~~~~{.cxx}
auto m = MaterialDataManager{b, 2u};
...
auto pr = std::vector<real>(2u);
extractInternalStateVariable(pr, m.s1, "HydrostaticPressure");
~~~~

# Issues solved

## Issue #95: Add an utility function to extract the value of an internal state variable

This feature is described in Section @sec:mgis:2.1:extractInternalStateVariable.

For more details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/95>.

## Issue #83: Add support for extended variable types 

This feature is described in depth in Section @sec:mgis:2.1:extended_types.

For more details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/83>.

## Issue #82: Add support for behaviour post-processings

This feature is described in depth in Section @sec:mgis:2.1:postprocessings.

For more details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/82>.

## Issue #78: Add support for behaviours which does not declare the temperature as the first external state variable

For more details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/78>.

## Issue #74: Support for tensorial external state variables

For more details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/74>.

## Issue #50: new feature request: initialize initial stress in `mgis.fenics` 

This feature is described in depth in Section @sec:mgis:2.1:initialize_functions.

For more details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/50>.