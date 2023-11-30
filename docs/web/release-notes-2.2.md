---
title: MFrontGenericInterfaceSupport Version 2.2
author: Thomas Helfer
date: 2023
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

# New features

## Declaration of the mass density in the `MaterialStateManager` class and associated free functions {#sec:mgis:2.2:mass_density}

The `MaterialStateManager` class now have a data member named `mass_density`.

The following free functions are available in `C++`:

- `setMassDensity` which sets the mass density using either a scalar
  value (uniform case) or a set of values (non uniform case).
- `isMassDensityDefined` which returns if the mass density has been
  defined.
- `isMassDensityUniform` which returns if the mass density is uniform.
  This function throws if the mass density has not been defined.

### `python` bindings

The `MaterialStateManager` class now have a `setMassDensity` member. A
corresponding free function is also available.

### `C` bindings

The following functions are now available:
`mgis_bv_material_state_manager_set_uniform_scalar_mass_density`,
`mgis_bv_material_state_manager_set_non_uniform_mass_density`,
`mgis_bv_material_state_manager_is_mass_density_defined`,
`mgis_bv_material_state_manager_is_mass_density_uniform`,
`mgis_bv_material_state_manager_get_uniform_mass_density` and
`mgis_bv_material_state_manager_get_non_uniform_mass_density`.

### `fortran` bindings

The following functions are now available:
`material_state_manager_set_uniform_scalar_mass_density`,
`material_state_manager_is_mass_density_defined` and
`material_state_manager_is_mass_density_uniform`.

## Support for retrieving the author, the date, the validator and the build identifier

The `BehaviourDescription` class now exposes the following data members:

- `author`: author of the `MFront` file,
- `date`: date associated with the `MFront` file,
- `validator`: validator of the `MFront` file,
- `build_id`: build identifier defined when compiling the `MFront` file.

## `C` bindings

The following functions are now available:

- `mgis_bv_behaviour_get_author`, which returns the author of the `MFront` file.
- `mgis_bv_behaviour_get_date`, which returns the date associated with the `MFront` file.
- `mgis_bv_behaviour_get_validator`, which returns the validator of the
  `MFront` file.
- `mgis_bv_behaviour_get_build_id`, which returns the build identifier
  defined when compiling the `MFront` file.

## `fortran` bindings

The `mgis_behaviour` module now exposes the following functions:

- `behaviour_get_author`, which returns the author of the `MFront` file.
- `behaviour_get_date`, which returns the date associated with the
  `MFront` file.
- `behaviour_get_validator`, which returns the validator of the `MFront`
  file.
- `behaviour_get_build_id`, which returns the build identifier defined
  when compiling the `MFront` file.

## `python` bindings

The `BehaviourDescription` class now exposes the following properties:

- `author`: author of the `MFront` file,
- `date`: date associated with the `MFront` file,
- `validator`: validator of the `MFront` file,
- `build_identifier` (or `build_id`): build identifier defined when compiling the `MFront` file.

# Issues solved

## Issue #125: Add support for retrieving the build indentifier and the validator

For more details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/125>.

## Issue #121: Pass the mass density from the `MaterialDataManager`

This feature is described in Section @sec:mgis:2.2:mass_density.

For more details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/121>.

