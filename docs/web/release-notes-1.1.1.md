% MFrontGenericInterfaceSupport Version 1.1.1
% Thomas Helfer

Version 1.1.1 is mostly a bug fix release of the 1.1
release.

# Issues solved

## Issue #46: Memory leakage in python bindings

The trouble is in the `mgis_convert_to_span` function which relies on
the `PyArray_FROM_OTF` function to retrieve the data associated with an
`numpy` array. I misread the documentation of this function and though
it only retrieves the pointer of the `numpy` object, but it does more:

- if its argument is a `numpy` array it increases the internal counter.
- if it is not, it tries to convert the object in a new `numpy` array.

In both cases, the counter of the returned object must be decreased at
the end of the function, which was not done. Hence the memory leak.

For details, see
<https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/46>.

## Issue #45: Stored and dissipated energies are not updated properly when updating the `MaterialDataManager`

The stored and dissipated energies were not updated when updating the
`MaterialDataManager`.

For details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/45>.
