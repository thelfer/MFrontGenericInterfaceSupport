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
- `negate`

#### Modifiers's generators

- `multiply_by_scalar`
- `divide_by_scalar`

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

### Modifiers dedicated to mechanics

- `hydrostatic_stress`
- `vmis`
- `principal_stress`

## Binary modifiers

-  `add`
-  `substract`
-  `multiply`
-  `divide`
-  `mean_value`
-  `inner_product`

### Binary modifiers on tensor-valuated evaluators

- `rotate`
- `rotate_backwards`

### Binary modifiers dedicated to mechanics

- `from_pk1_to_cauchy`
