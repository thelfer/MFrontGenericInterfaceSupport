---
title: MFrontGenericInterfaceSupport Version 1.2.1 
author: Thomas Helfer, Jérémy Bleyer
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

# Issues solved

## Issue #66: `mgis.fenics`: Missing update of external state variables when seen as UFL expression

Declaring external state variables as `UFL` expressions were not properly updated between time steps.
