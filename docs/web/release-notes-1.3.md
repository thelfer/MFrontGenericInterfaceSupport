---
title: MFrontGenericInterfaceSupport Version 1.3 
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

Version 1.3 of `MFrontGenericInterfaceSupport` is compatible with the
Version 4.0 of `TFEL/MFront`.

# Known incompatibilities

The deformation gradients are automatically initialized to identity by
the `MaterialStateManager` classes unless the gradients are initialized
using an externally allocated memory.
