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

## Removed function in `C++`

The following functions were removed:

- `getUniformMaterialProperty`
- `getNonUniformMaterialProperty`
- `getUniformExternalStateVariable`
- `getNonUniformExternalStateVariable`

## Removed function in `C`

The following functions were removed:

- `mgis_bv_material_state_manager_get_uniform_material_property`
- `mgis_bv_material_state_manager_get_non_uniform_material_property`
- `mgis_bv_material_state_manager_get_uniform_external_state_variable`
- `mgis_bv_material_state_manager_get_non_uniform_external_state_variable`

## Removed function in `fortran`

The following functions were removed from the `mgis_behaviour` module:

- `material_state_manager_get_uniform_material_property`
- `material_state_manager_get_non_uniform_material_property`
- `material_state_manager_get_uniform_external_state_variable`
- `material_state_manager_get_non_uniform_external_state_variable`

# Issues solved

## Issue #78: Add support for behaviours which does not declare the temperature as the first external state variable

For more details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/78>.

## Issue #74: Support for tensorial external state variables

For more details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/74>.