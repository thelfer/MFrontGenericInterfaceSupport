---
title: MFrontGenericInterfaceSupport Version 2.2
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

# Issues solved

## Issue #121: Pass the mass density from the `MaterialDataManager`

This feature is described in Section @sec:mgis:2.2:mass_density.

For more details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/121>.
