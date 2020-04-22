% Introduction to the FEniCS interface to MGIS
% Bleyer Jeremy
% April 5th 2020

\newcommand{\bsig}{\boldsymbol{\sigma}}
\newcommand{\bg}{\mathbf{g}}
\newcommand{\dx}{\,\text{dx}}
\newcommand{\dS}{\,\text{dS}}
\newcommand{\TT}{\mathbb{T}}

These pages are intended to describe the `FEniCS` interface to `MGIS`
available in the `mgis.fenics` Python module. This module has been
developed in the spirit of providing sa upport to the implementation of
generalized nonlinear behaviours i.e. including multiple gradients and
dual flux variables. It can also be used as a simple interface to
standard mechanical behaviours (see the plasticity demos). We briefly
describe [here](#a-brief-overview-of-the-mgis.fenics-module) the general
concepts underlying the module implementation.

# Documented demos

The provided documented demos have been designed to progressively
illustrate the use of the interface and the versatility of the approach
when implementing complex generalized behaviours, both on the `MFront`
and `FEniCS` sides. We recommend browsing the demos in the following
order:

* [Stationnary non-linear heat
  transfer](mgis_fenics_nonlinear_heat_transfer.html)
* [Stationnary non-linear heat transfer: 3D problem and performance
  comparisons](mgis_fenics_nonlinear_heat_transfer_3D.html)
* [Small-strain von Mises
  elastoplasticity](mgis_fenics_small_strain_elastoplasticity.html)
* [Finite-strain elastoplasticity within the logarithmic strain
  framework](mgis_fenics_finite_strain_elastoplasticity.html)
* [Multiphase model for fiber-reinforced
  materials](mgis_fenics_multiphase_model.html)
* [Phase-field approach to brittle
  fracture](mgis_fenics_phase_field.html)
* [Transient non-linear heat
  equation](mgis_fenics_transient_nonlinear_heat_equation.html)

# A brief overview of the `mgis.fenics` module

This module has been developed based on the concepts exposed in the [Elasto-plastic analysis implemented using the MFront code generator](https://comet-fenics.readthedocs.io/en/latest/demo/plasticity_mfront/plasticity_mfront.py.html) demo published on [Numerical Tours of Computational Mechanics using FEniCS](https://comet-fenics.readthedocs.io). We suggest reading first this demo as an introduction to this module implementation concepts. In particular, the module relies on the notion of `Quadrature` function spaces to express the constitutive relation at the quadrature points level. Their definition will however be handled directly inside the [`MFrontNonlinearProblem`](#the-mfrontnonlinearproblem-class) class. The module will also rely on the `NewtonSolver` available in FEniCS instead of programming the Newton method manually. We also provide access to non-linear optimisation solvers (see the [`MFrontOptimisationProblem` class](#the-mfrontoptimisationproblem-class)). Finite-strain behaviours are now supported as well as parallel computations.

## The `MFrontNonlinearMaterial` class

This class handles the loading of a MFront behaviour through MGIS. In particular, it will contain the following important attributes:

* `behaviour`: an instance of MGIS `Behaviour` class which handles all the information about a specific MFront behaviour. It is created by the `load` function which takes the path to a library, the name of a behaviour and a modelling hypothesis.

  Note that before loading the behaviour, it is checked if the behaviour is a finite-strain one or not. In the former case, specific finite-strain options are used when calling `load`. Such options specify that the stress measure will be post-processed by MGIS from Cauchy to First Piola Kirchhoff (PK1) stress and that the tangent operator will be given by $\dfrac{\partial \text{PK1}}{\partial F}$ (`DPK1_DF`).

* `data_manager`: an instance of MGIS `MaterialDataManager` class which handles a bunch of integration points. It is instantiated using `behaviour`and the number of integration points

A set of helper functions enables to retrieve information about MFront object names and sizes of the form
```python
self.get_{1}_{2}
```
where `{1}` can be any of `material_property`, `external_state_variable`, `internal_state_variable`, `gradient`, `flux` or `tangent_block` and `{2}` can be any of `names` or `sizes`.


## The `MFrontNonlinearProblem` class

This class handles the definition and resolution of a nonlinear problem associated with a `MFrontNonlinearMaterial`. Its main attributes are:

* `u`: the unknown mechanical field $u$ (a `dolfin.Function`)

* `material`: the associated `MFrontNonlinearMaterial`

* `state_variables`: a dictionary of internal and external state variables. External state variables are represented as `dolfin` objects whereas internal state variables are represented as [Quadrature functions](#flux-and-gradient-objects).

* `gradients`: a dictionary of `Gradient` objects

* `fluxes`: a dictionary of `Flux` objects

* `residual`: the nonlinear residual $F(u)$

* `tangent_form`: the tangent bilinear form associated with the residual

* `solver`: the non-linear solver (default is `NewtonSolver`)

Its main methods are:

* `register_gradient`: registers a MFront gradient with a UFL expression (see [the registration concept](#the-registration-concept)).

* `register_external_state_variable`: registers a MFront external state variable with a UFL expression (see [the registration concept](#the-registration-concept)).

* `set_loading`: defines the external forces linear form $L(v)$ in the default residual expression 

* `initialize`: initializes the functions associated with gradients, fluxes, external and internal state variables objects and the corresponding tangent blocks. All gradients and external state variables must have been registered first. [Automatic registration](#the-registration-concept) is performed at the beginning of the method call. This method has to be called once before calling `solve`.

* `update_constitutive_law`: performs the consitutive law update. Currently, this method is called for all quadrature points before assembling the residual or tangent form. Ideally, this should be done locally during the assembly (see [the module current limitations](#current-limitations)).

* `solve`: solves the associated non-linear problem using `solver`

* `get_flux`: returns the function associated with a flux

* `get_state_variable`: returns the function associated with an internal state variable

### Residual and tangent bilinear form
By default, the nonlinear residual is assumed to take the following form: Find $u\in V$ such that:

$$
\begin{equation}
 \sum_{i=1}^p \int_{\Omega}\bsig_i(u)\cdot\delta\bg_i(v) \dx - L(v) = 0 \quad \forall v\in V 
\tag{residual}
\label{residual}
\end{equation}
$$

where the $\bsig_i(u)$ are a set of **fluxes** as defined in the MFront behaviour using `@Flux` and $\bg_i$ are the corresponding set of **gradients**, declared using `@Gradient`. $\delta \bg_i(v)$ denotes the directional derivative of $\bg_i$ along direction $v$, a `TestFunction` of the function space $V$, and $L(v)$ is a linear form which can be expressed using standard UFL operators. From the mechanical point of view, this residual expresses the balance between internal and external forces.

This generic form is suitable for most quasi-static mechanical problems but it does not necessarily encompass all possible situations. In particular, evolution equations such as transient heat transfer cannot be expressed in this form. However, this is not a limitation since the residual can also be redefined explicitly by the user (see [Transient non-linear heat equation](demos/transient_nonlinear_heat_equation.html)).

Either for the default or a user-defined one, the tangent operator associated with the residual is computed automatically using either the UFL symbolic derivation `ufl.derivative` for simple terms (e.g. involving the unknown field $u$) or using tangent operator blocks defined in the MFront behaviour for the nonlinear fluxes and internal state variables. More precisely, the variation of each flux is given by:

$$
\dfrac{\partial \bsig_i}{\partial u} = \sum_{j\in \text{blocks}(i)} \TT^{\bsig_i}_{\bg_j}\cdot \delta \bg_j(u)
$$
where the flux $\bsig_i$ is assumed to depend on the gradients $\bg_j$ for $j\in \text{blocks}(i)$ (which at least contains $i$ itself in general) and $\TT^{\bsig_i}_{\bg_j} = \dfrac{\partial \bsig_i}{\partial \bg_j}$ is the tangent operator associated with the corresponding block. Variations of internal state variables are computed in the same manner.

In the default case, the tangent bilinear form therefore reads as:

$$
a_\text{tangent}(u, v) = \sum_{i=1}^p \sum_{j\in\text{blocks}(i)}\int_{\Omega}\delta\bg_j(u)\cdot \TT^{\bsig_i}_{\bg_j} \cdot\delta\bg_i(v) \dx 
$$

### `Flux` and `Gradient` objects

Two helper classes have been defined to handle flux and gradient objects:

* the `Gradient` class provides a representation of MFront gradient objects. Its main purpose is to provide the corresponding UFL expression, linking MFront and FEniCS concepts. It also handles:

    - the reshaping from UFL tensorial representation to MFront vectorial conventions
    - the symbolic expression of the gradient variation (directional derivative)
    - the representation as a Quadrature function

    This class is intended for internal use only. Gradient objects must be declared by the user using the [registration concept](#the-registration-concept).

* the `Flux` class provides a representation of MFront flux objects. Its main purpose is to store the flux values in the form of a Quadrature function, as well as the tangent block structure of the flux, each of them also represented by a Quadrature function.

* the `InternalStateVariable` class is the same as the `Flux` class but for internal state variables. Both classes derive from the abstract `QuadratureFunction` class.

### The registration concept

Inspecting the default case $\eqref{residual}$, MFront provides access to the flux names, shapes and values when performing the constitutive update and also to the corresponding gradient. The definition of the tangent operator blocks inside the MFront behaviour also gives access to the block structure $\text{blocks}(i)$ for each flux. The remaining information which must be provided from the FEniCS side are the unknown field $u$ and its discretization space $V$, the chosen integration measure $\dx$ (through the `quadrature_degree` keyword) and, finally, the expression of each declared gradients $\bg_i$ in terms of the unknown field $u$. This step is what we call *registration* of each gradient which will be discussed in depth in the demos. Let us just mention that the gradients can registered using the `register_gradient` method of the `MFrontNonlinearProblem` class or via an automatic procedure if the gradient name matches predefined common gradient objects e.g. `"Strain"`, `"TemperatureGradient"`, etc. 

## The `MFrontOptimisationProblem` class

This class is a sibing to the `MFrontNonlinearProblem` class. It enables to solve nonlinear optimisation problems, especially bound-constrained problems using the default `PETScTAOSolver` of the form:

$$
\min_{b_l \leq u \leq b_u}  f(u)
$$
where $b_l$ (resp. $b_u$) denotes a lower (resp. upper) bound on the optimisation variable $u\in V$.

By default, the objective function $f(u)$ corresponds to the material total energy (stored + dissipated) computed from the `get_total_energy()` method. The optimisation problem requires the definition of the gradient $F(u)$ which, by default, corresponds to $\eqref{residual}$ and its jacobian which is computed as discussed before.

## Current limitations

The module has been developed using FEniCS version 2019.1.0. An important reimplementation will be planned once the [`dolfinx` project](https://github.com/FEniCS/dolfinx) will officially release a stable version. In particular, it will aim at fixing the following current limitations:

* the constitutive update is performed before any assembly procedure of the tangent and residual forms. This adds an extra cost of looping over quadrature points and, more importantly, an important memory cost since all tangent blocks at all quadrature points must be saved. The `dolfinx` project should offer the possibity of using custom assemblers in which constitutive integration should be possible at the local assembly level

* multi-materials are not completely supported yet. More precisely, spatially varying material properties can be defined but it is not possible to define two different constitutive behaviours on two disting parts of the mesh. This feature has not been supported since it is not possible to define functions on sub-meshes yet. This should also be available soon in the next developments.

* memory transfers between FEniCS and MGIS objects have not been optimized

