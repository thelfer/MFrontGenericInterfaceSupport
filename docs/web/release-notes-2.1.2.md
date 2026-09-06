---
title: MFrontGenericInterfaceSupport Version 2.1.2
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

## Issue 252: Fix conditional tests  using `math_errhandling`: this may is defined on MacOS by a macro calling an intrinsic which is not `constexpr`

For more details, see <https://github.com/thelfer/MFrontGenericInterfaceSupport/issues/252>
