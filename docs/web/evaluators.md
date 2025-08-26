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

Evaluators and modifiers are two of the main concepts of `MGIS/Function`
(see [this page](functions.html) for an introduction).

# Unary modifiers

- `view`: the `view` modifier takes a dynamic function or evaluator and
  return an evaluator with a fixed number of components. The `view`
  modifier takes the expected number of components as argument.

  ~~~~{.cxx}
  // take a view on some values with a runtime size
  const auto f = FunctionEvaluator<BasicLinearSpace>(space, values, 1);
  // make a view with a fixed number of components and apply the `transform` modifier
  const auto e = view<1>(f) | transform([](const real x) { return x * x; });
  ~~~~
- `negate`: the `view` modifier takes an evaluator which returns
  opposite values of the initial evaluator.

  ~~~~{.cxx}
  const auto f =
      FunctionEvaluator<BasicLinearSpace>(space, values, 2) | as_tvector<2>;
  const auto n = f | negate;
  ~~~~

### Modifiers's generators

- `multiply_by_scalar`: the `multiply_by_scalar` takes an real number
  \(alpha\) and returns a modifier which takes an evaluator and returns
  an evaluator which multiplies the values of the original evaluator bu
  \(\alpha\):

  ~~~~{.cxx}
  const auto e = f | multiply_by_scalar(4);
  ~~~~
- `divide_by_scalar`: the `divide_by_scalar` takes an real number
  \(alpha\) and returns a modifier which takes an evaluator and returns
  an evaluator which divides the values of the original evaluator bu
  \(\alpha\):

  ~~~~{.cxx}
  const auto e = f | divide_by_scalar(4);
  ~~~~

## Tensors

This section is devoted to modifiers related to tensors.

### Views

- `as_scalar`:  create an evaluator which returns a scalar.


  ~~~~{.cxx}
  const auto s = f | as_scalar;
  ~~~~
- `as_stensor`: create an evaluator which returns symmetric tensorial
  values. This
  modifier takes an integer as template parameter, the space dimension.

  ~~~~{.cxx}
  const auto trace_values = f2 | as_stensor<2> | trace;
  ~~~~
- `as_tensor`: create an evaluator which returns unsymmetric tensorial
  values.  This
  modifier takes an integer as template parameter, the space dimension.
- `as_fsarray`: create an evaluator which returns an (finite size)
  array. This modifier takes an integer as template parameter, the size
  of the array.
- `as_tvector`: create an evaluator which returns vectorial values. This
  modifier takes an integer as template parameter, the size of the
  vector.
- `as_matrix`: create an evaluator which returns matrix values. This
  modifier takes two integers as template parameters: the number of rows
  and the number of columns.

### Modifiers on tensor-valuated evaluators

- `eigen_values`: modifier creating an evaluator of the eigen values.
  The evaluator on input must return symmetric tensor values.
- `trace`:
- `det`:

## Modifiers dedicated to mechanics

- `hydrostatic_stress`: modifier creating an evaluator of the
  hydrostatic pressure \(p\) defined as:

  \[
  p = \mathrm{trace}\left(\underline{\sigma}\right/3
  \]
- `vmis`: modifier creating an evaluator of the von Mises stress.
- `principal_stress`: modifier creating an evaluator of the principal
  stress. This is an alias for the `eigen_values` modifier

# Binary modifiers

-  `add`:
-  `substract`:
-  `multiply`:
-  `divide`:
-  `mean_value`:
-  `inner_product`:

## Binary modifiers on tensor-valuated evaluators

- `rotate`:
- `rotate_backwards`:

## Binary modifiers dedicated to mechanics

- `from_pk1_to_cauchy`:
