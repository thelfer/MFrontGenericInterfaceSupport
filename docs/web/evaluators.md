---
title: Evaluators and modifiers
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

# Evaluators

## Unary modifiers

- `view`

### Tensors

#### Views

- `as_stensor`
- `as_tensor`
- `as_fsarray`
- `as_tvector`
- `as_matrix`

#### Modifiers on tensor-valuated evaluators

- `eigen_values`
- `trace`
- `det`

### Mechanical evaluators

- `hydrostatic_stress`
- `vmis`
- `principal_stress`

## Binary modifiers

#### Modifiers on tensor-valuated evaluators

- `rotate`
- `rotate_backwards`
