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
