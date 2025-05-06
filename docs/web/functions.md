---
title: Functions
author: Thomas Helfer
date: 2025
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

Functions, called fields in other contexts, describe values on attached
to the elements of a discretization space.

# Discretization spaces

In order to adapt to the structures used by as many solvers as possible,
spaces in `mgis` are defined using a set of `C++` concepts. Let
`SpaceType` be a type.

- `SpaceType` satisfies the `SpaceConcept` concept if:
  - it exposes a `size` method. This method shall return the total
    number of elements in the discretization space.
  - the `SpaceTraits<SpaceType>` shall expose an integral type named
    `size_type`.
- `SpaceType` satisfies the `LinearSpaceConcept` concept if:
  - it satisfies the `SpaceConcept`
- the `QuadratureSpace` concept

# Point-major or field major storage

see [petscsection_2025] for a discussion.

# The `Function` class

# The `FunctionView` class
