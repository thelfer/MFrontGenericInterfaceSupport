% MFrontGenericInterfaceSupport Version 1.2 
% Thomas Helfer
% 17/12/2019

Version 1.2 of `MFrontGenericInterfaceSupport` is compatible with the
Version 3.4 of `TFEL/MFront`.

# New functionalities

## C bindings

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