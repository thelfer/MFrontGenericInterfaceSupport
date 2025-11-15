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

Functions, called fields in other contexts, describe a mapping between
elements of a (discretization) space to values. Those elements can be
associated to quadrature points (also called integration points),
vertices (nodes), facets, cells, voxels, depending on the discretization
method considered.

# An overview of `MGIS/Function` {#sec:mgis:function:overview}

Let us consider a function `pk1` returning the first Piola-Kirchhoff
stress in the material frame, `MGIS/Function` allows to create an
**evaluator** of the Cauchy stress in the global frame as follows:

~~~~{.cxx}
const auto Ft  = F | as_tensor<3>;
const auto sig_ev = pk1 | as_tensor<3> | from_pk1_to_cauchy(Ft) | rotate_backwards(R);
~~~~

where `F` is an evaluator the deformation gradient in the material
frame, `R` is the rotation matrix from the global frame to the material
frame (a function returning a different rotation matrix for each element
can also be used).

This evaluator, named `sig_ev`, only describes the operations to be
performed to calculate the Cauchy stress and does not make any
calculation at this stage. `as_tensor<3>`, `from_pk1_to_cauchy` and the
object returned by `rotate_backwards(R)` are called *modifiers* in
`MGIS/Function`. The available operators are described in [this
page](evaluators.html). The pipe operator `|` is used to chain the
operations.

Hence, the following operations are chained:

- interpret the values returned by `pk1` as tridimensional tensors
  provided by the `TFEL/Math` library,
- compute the Cauchy stress in the material frame,
- rotate the Cauchy stress in the global frame by using the transpose
  of the rotation matrix `R`.

`sig_ev `exposes the same space as the `pk1` function (accessible from the
`getSpace` method) and the same call operators (`4` signatures are
possible, as described below).

For instance, if the values of `pk1` can be retrieved by using indexes
(that are not necessarily integers) associated with the elements of the
space, then the values of the Cauchy stress can be retrieved by the same
indexes. Hence, if i is such an index, `pk1(i)` returns the first
Piola-Kirchhoff stress in the material frame, and `sig_ev(i)` returns the
Cauchy stress in the global frame.

For spaces with very simple structure, `MGIS/Function` provides
algorithms, such as the `assign` algorithm which allows to assign the
values of an evaluator to a function. Let `sig` be a function, the
following code assigns the values calculated by `sig_ev` 

Those concepts now being defined, `MGIS/Function` can now be defined as
set of algorithms, evaluators and modifiers used to operate and
calculate values on functions. `MGIS/Function` also provides basic
implementations of spaces and functions.

# Discretization spaces

The description of a discretization space is dependent of the solver
considered.

A space must provide a way to iterate over its elements:

- in many cases, iterating over all the elements require a simple loop.
  Such a case is called an *element space* in `MGIS/Function` and each
  element is associated with a single index. The index type is often an
  integer, but this is sometimes not the case. `MGIS/Function` assumes
  that a function defined on an element space can return the values
  associated with an element using its index.
- quadrature points deserve a special treatment as they are generally
  considered as part of a cell (or finite element). In such a case,
  iterating over all integrations points require a loop over the finite
  elements and an inner loop over the integration points of this
  element. A quadrature point is then associated with two indexes, one
  for the cell, the other being an index inside the element. Such a
  space is called *quadrature space* in `MGIS/Function`. `MGIS/Function`
  assumes that a function defined on a quadrature space can return the
  values associated with an quadrature point using its indexes.

> **Notes**
> 
> A space describing quadrature points can be either an *element space*
> and/or a quadrature space.

## `C++` Concepts defined `MGIS/Function`

In order to adapt to the structures used by as many solvers as possible,
spaces in `mgis` are defined using a set of `C++` concepts:
`SpaceConcept`, `ElementSpaceConcept`, `QuadratureSpaceConcept`,
`FunctionalSpace`, `LinearElementSpaceConcept` and
`LinearQuadratureSpaceConcept`.

Let `SpaceType` be a type.

- `SpaceType` satisfies the `SpaceConcept` concept if:
  - it exposes an integral type named `size_type` and a `size` method
    (returning a `size_type`). This method shall return the total number
    of elements in the discretization space.
- `SpaceType` satisfies the `ElementSpaceConcept` concept if:
  - it exposes an `element_index_type` type.
- `SpaceType` satisfies the `LinearElementSpaceConcept` concept if:
  - it satisfies the `ElementSpaceConcept`,
  - it exposes a boolean value named `linear_element_indexing` which
    must be `true`,
  - elements of the space are indexed by integers ranging from from 0 to
    `size()-1`. The `element_index_type` must be identical to
    `size_type`.
- `SpaceType` satisfies the `QuadratureSpaceConcept` concept if:
  - it satisfies the `SpaceConcept`,
  - it exposes a method named `getNumberOfCells` which returns a
    `size_type`,
  - it exposes an `cell_index_type` type,
  - it exposes an `quadrature_point_index_type` type. In `MGIS/Function`,
    this type must be an integer type,
  - it exposes a method named `getNumberOfQuadraturePoints` which returns
    the number of quadrature points of a given cell.
- `SpaceType` satisfies the `LinearQuadratureSpaceConcept` concept if:
  - it satisfies the `QuadratureSpaceConcept`,
  - it exposes a boolean value named `cell_element_indexing` which
    must be `true`,
  - cells of the space are indexed by integers ranging from from 0 to
    `getNumberOfCells()-1`. The `cell_index_type` must be identical to
    `size_type`,
  - it exposes a method `getQuadraturePointOffset` which must be a
    bijection between all valid cell and quadrature points indices and
    the range `[0:size()-1]`.

A functional space is either an element space and/or a quadrature space.

# Functions

## `C++` Concepts defined `MGIS/Function`

In order to adapt to the structures used by as many solvers as possible,
spaces in `mgis` are defined using a set of `C++` concepts:
`FunctionConcept`, `ElementFunctionConcept`, `QuadratureFunctionConcept`.

## Element-major or function-major storage {#sec:mgis:functions:storage_policy}

The documentation of the `PETsC` library discusses various patterns used
in the litterature to store the values of a function
[@petscsection_2025]. This discussion distinguishes so-called
"point-major" or "field-major" storages. In `MGIS/Function`, we will use
the terms "element-major" and "function-major" for consistency.

In a field-major pattern, the values a single-valuated field are usually
stored as follows:

~~~~
+-----------++-----------++---------++-----------+
| Element 1 || Element 2 || ....... || Element N |
+-----------++-----------++---------++-----------+
~~~~

Multi-component values can be stored:

- in an interleaved manner (all the components associated with one
  element are stored contiguously):

  ~~~~
  | <------------- Element 1 ------------> | .... | <------------- Element N ------------> |
  +-------------++---------++--------------+------+-------------++---------++--------------+
  | Component 1 || ....... || Component Nc | .... | Component 1 || ....... || Component Nc |
  +-------------++---------++--------------+------+-------------++---------++--------------+
  ~~~~

- in a non-interleaved manner (all the values associated with a given
  component of the function are stored continuously):

  ~~~~
  | <------------------ Component 1 -------------> | .... | <------------- Component Nc -----------------> |
  +-----------++-----------++---------++-----------+------+-----------++-----------++---------++-----------+
  | Element 1 || Element 2 || ....... || Element N |......| Element 1 || Element 2 || ....... || Element N |
  +-----------++-----------++---------++-----------+------+-----------++-----------++---------++-----------+
  ~~~~

Non interleaved storage allows coalescent memory accesses, which may be
beneficial on some architectures, including GPUs.

In a element-major pattern, values usually are stored as follows:

~~~~
| <------------------ Element 1 ------------------> | .... | <------------------- Element N -----------------> |
+------------++------------++---------++------------+      +------------++------------++---------++------------+
| Function 1 || Function 2 || ....... || Function F | .... | Function 1 || Function 2 || ....... || Function F |
+------------++------------++---------++------------+      +------------++------------++---------++------------+
~~~~

This element major pattern is used by the `MaterialStateManager` class
to store gradients, thermodynamic forces and internal state variables.

## The `Function` class

The `Function` class is an implementation of the field-major function
based on a linear functional space. Multiple components are stored in an
interleaved manner. The `Function` class has two template parameters:

- the functional space,
- the number of components which defaults to `mgis::dynamic_extent`,
  meaning that the number of components can be choosen at runtime.

## The `FunctionView` class

The `FunctionView` is a lightweight class that allows to make a view on
top of a contiguous memory which acts as a function over a given linear
functional space as follows:

~~~~
|--------------------------------------------------------------------------------------------------------------------------|
<-                         Raw data                                                                                       ->
|--------------------------------------------------------------------------------------------------------------------------|
<- Data of the first element          --><- Data of the second element        -->....<- Data of the Nth element          -->                    
|---------------|xxxxxxxxxxxxxxxxxxxxxxx||---------------|xxxxxxxxxxxxxxxxxxxxxx|....|---------------|xxxxxxxxxxxxxxxxxxxxx|
<-function data->                        <-function data->                           <-function data->                      
^               ^                       ^                ^                      ^                    ^                     ^
|               |                       |                |                      |                    |                     |
            data_size                   |            data_size                  |                data_size                 |
                                    data_stride                              data_stride                               data_stride
~~~~

This figure shows that:

- The components of the function are interleaved, i.e. stored
  contiguously. The number of components is called the data size on the
  figure.
- The values of the function are not contiguous in memory, i.e. the
  values associated with the first component of two consecutive elements
  are separated by a constant stride, named `data_stride`. The data
  between `data_size` and `data_stride`, marked by an `x` in the figure,
  are not accessible by the view.

By definition of a view, a function view does not handle the lifetime of
the memory on which it is built. The user is responsible for:

- ensuring that the function view is destroyed before the underlying
  memory (or a least not used after this memory is released),
- releasing the underlying memory.

This data structure is very close to an element-major storage (see
Section @sec:mgis:functions:storage_policy), as the one used by the
`MaterialStateManager` class. Indeed, the primary intent of the
`FunctionView` class is to manipulate as each gradient, thermodynamic
force or internal state variable as if it were stored in an interleaved
element-major function.

The `FunctionView` class has two template parameters:

- the functional space,
- a data layout containing both the number of components (data size) and
  the data stride. If the data size or the date stride is equal
  `mgis::dynamic_extent`, which is the default, its value must be given
  at runtime.

The `FunctionView` class matches the `FunctionEvaluator` concept (See
Section @sec:mgis:function:overview and [this page](evaluators.html) for
details).

## Tensorial functions

### The `TensorView` class

The `TensorView` class allows to make views which returns tensorial
objects from functions returning data in contiguous memory.

A `TensorView` is generally created by combining a `FunctionView` which
a modifier such as `as_stensor` (See [this page](evaluators.html) for
details).

### The `CoalescedMemoryAccessTensorView` class

Coalescent memory access refers to an access pattern where each
components of a function is accessed through its dedicated memory
location (see the non-interleaved storage pattern described previously).

The `CoalescedMemoryAccessTensorView` class allows to make a tensorial
function view from scalar function view on each components, as follows:

~~~~{.cxx}
real values[8] = {1, 10,  // components XX
                  2, 20,  // components YY
                  3, 30,  // components ZZ
                  4, 40}; // components XY
const auto ne = size_type { 2 };
auto space = BasicLinearSpace{ne};
auto c_xx = ScalarFunctionView{space, std::span{values, ne}};
auto c_yy = ScalarFunctionView{space, std::span{values + ne, ne}};
auto c_zz = ScalarFunctionView{space, std::span{values + 2 * ne, ne}};
auto c_xy = ScalarFunctionView{space, std::span{values + 3 * ne, ne}};
auto f = CoalescedMemoryAccessTensorView<BasicLinearSpace,
                                         tfel::math::stensor<2, real>>{
            std::array{c_xx, c_yy, c_zz, c_xy}};
~~~~