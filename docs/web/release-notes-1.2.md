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

## Issue #30: Variables array size depending on type and modelling hypothesis

The `get_variable_size` function is now available in the
`mgis_behaviour` module. This function returns the size of a variable
(i.e. the number of components) for the given hypothesis.

For details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/30>

# References