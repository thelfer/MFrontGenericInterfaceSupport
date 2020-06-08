% Stationnary non-linear heat transfer
> **Source files:**
>
> * Jupyter notebook: [mgis_fenics_nonlinear_heat_transfer.ipynb](https://gitlab.enpc.fr/navier-fenics/mgis-fenics-demos/raw/master/demos/nonlinear_heat_transfer/mgis_fenics_nonlinear_heat_transfer.ipynb)
> * Python file: [mgis_fenics_nonlinear_heat_transfer.py](https://gitlab.enpc.fr/navier-fenics/mgis-fenics-demos/raw/master/demos/nonlinear_heat_transfer/mgis_fenics_nonlinear_heat_transfer.py)
> * MFront behaviour file: [StationaryHeatTransfer.mfront](https://gitlab.enpc.fr/navier-fenics/mgis-fenics-demos/raw/master/demos/nonlinear_heat_transfer/StationaryHeatTransfer.mfront)

# Description of the non-linear constitutive heat transfer law

The thermal material is described by the following non linear Fourier
Law:

$$
\mathbf{j}=-k\left(T\right)\,\mathbf{\nabla} T
$$

where $\mathbf{j}$ is the heat flux and $\mathbf{\nabla} T$ is the
temperature gradient.

## Expression of the thermal conductivity

The thermal conductivity is assumed to be given by:

$$
k\left(T\right)={\displaystyle \frac{\displaystyle 1}{\displaystyle A+B\,T}}
$$

This expression accounts for the phononic contribution to the thermal
conductivity.

## Derivatives

As discussed below, the consistent linearisation of the heat transfer
equilibrium requires to compute:

-   the derivative
    ${\displaystyle \frac{\displaystyle \partial \mathbf{j}}{\displaystyle \partial \mathbf{\nabla} T}}$
    of the heat flux with respect to the temperature gradient.
    ${\displaystyle \frac{\displaystyle \partial \mathbf{j}}{\displaystyle \partial \mathbf{\nabla} T}}$
    is given by: $$
      {\displaystyle \frac{\displaystyle \partial \mathbf{j}}{\displaystyle \partial \mathbf{\nabla} T}}=-k\left(T\right)\,\matrix{I}
    $$
-   the derivative
    ${\displaystyle \frac{\displaystyle \partial \mathbf{j}}{\displaystyle \partial T}}$
    of the heat flux with respect to the temperature.
    ${\displaystyle \frac{\displaystyle \partial \mathbf{j}}{\displaystyle \partial T}}$
    is given by: $$
      {\displaystyle \frac{\displaystyle \partial \mathbf{j}}{\displaystyle \partial T}}=-{\displaystyle \frac{\displaystyle \partial k\left(T\right)}{\displaystyle \partial T}}\,\mathbf{\nabla} T=B\,k^{2}\,\mathbf{\nabla} T
    $$

# `MFront`’ implementation

## Choice of the the domain specific language

Every `MFront` file is handled by a domain specific language (DSL), which
aims at providing the most suitable abstraction for a particular choice
of behaviour and integration algorithm. See `mfront mfront --list-dsl`
for a list of the available DSLs.

The name of DSL’s handling generic behaviours ends with
`GenericBehaviour`. The first part of a DSL’s name is related to the
integration algorithm used.

In the case of this non linear transfer behaviour, the heat flux is
explicitly computed from the temperature and the temperature gradient.
The `DefaultGenericBehaviour` is the most suitable choice:

``` cpp
@DSL DefaultGenericBehaviour;
```

## Some metadata

The following lines define the name of the behaviour, the name of the
author and the date of its writing:

``` cpp
@Behaviour StationaryHeatTransfer;
@Author Thomas Helfer;
@Date 15/02/2019;
```

## Gradients and fluxes

Generic behaviours relate pairs of gradients and fluxes. Gradients and
fluxes are declared independently but the first declared gradient is
assumed to be conjugated with the first declared fluxes and so on…

The temperature gradient is declared as follows (note that Unicode characters are supported):

``` cpp
@Gradient TemperatureGradient ∇T;
∇T.setGlossaryName("TemperatureGradient");
```

Note that we associated to `∇T` the glossary name `TemperatureGradient`.
This is helpful for the calling code.

After this declaration, the following variables will be defined:

-   The temperature gradient `∇T` at the beginning of the time step.
-   The increment of the temperature gradient `Δ∇T` over the time step.

The heat flux is then declared as follows:

``` cpp
@Flux HeatFlux j;
j.setGlossaryName("HeatFlux");
```

In the following code blocks, `j` will be the heat flux at the end of
the time step.

## Tangent operator blocks

By default, the derivatives of the gradients with respect to the fluxes
are declared. Thus the variable `∂j∕∂Δ∇T` is automatically declared.

However, as discussed in the next section, the consistent linearisation
of the thermal equilibrium requires to return the derivate of the heat
flux with respect to the increment of the temperature (or equivalently
with respect to the temperature at the end of the time step).

``` cpp
@AdditionalTangentOperatorBlock ∂j∕∂ΔT;
```

## Parameters

The `A` and `B` coefficients that appears in the definition of the
thermal conductivity are declared as parameters:

``` cpp
@Parameter real A = 0.0375;
@Parameter real B = 2.165e-4;
```

Parameters are stored globally and can be modified from the calling
solver or from `python` in the case of the coupling with `FEniCS`
discussed below.

## Local variable

A local variable is accessible in each code blocks.

Here, we declare the thermal conductivity `k` as a local variable in
order to be able to compute its value during the behaviour integration
and to reuse this value when computing the tangent operator.

``` cpp
@LocalVariable thermalconductivity k;
```

## Integration of the behaviour

The behaviour integration is straightforward: one starts to compute the
temperature at the end of the time step, then we compute the thermal
conductivity (at the end of the time step) and the heat flux using the
temperature gradient (at the end of the time step).

``` cpp
@Integrator{
  // temperature at the end of the time step
  const auto T_ = T + ΔT;
  // thermal conductivity
  k = 1 / (A + B ⋅ T_);
  // heat flux
  j = -k ⋅ (∇T + Δ∇T);
} // end of @Integrator
```

## Tangent operator

The computation of the tangent operator blocks is equally simple:

``` cpp
@TangentOperator {
  ∂j∕∂Δ∇T = -k ⋅ tmatrix<N, N, real>::Id();
  ∂j∕∂ΔT  =  B ⋅ k ⋅ k ⋅ (∇T + Δ∇T);
} // end of @TangentOperator 
```
# `FEniCS` implementation

We consider a rectanglar domain with imposed temperatures `Tl` (resp. `Tr`) on the left (resp. right) boundaries. We want to solve for the temperature field `T` inside the domain using a $P^1$-interpolation. We initialize the temperature at value `Tl` throughout the domain.


```python
%matplotlib notebook
import matplotlib.pyplot as plt
from dolfin import *
import mgis.fenics as mf

length = 30e-3
width = 5.4e-3
mesh = RectangleMesh(Point(0., 0.), Point(length, width), 100, 10)

V = FunctionSpace(mesh, "CG", 1)
T = Function(V, name="Temperature")

def left(x, on_boundary):
    return near(x[0], 0) and on_boundary
def right(x, on_boundary):
    return near(x[0], length) and on_boundary

Tl = 300
Tr = 800
T.interpolate(Constant(Tl))

bc = [DirichletBC(V, Constant(Tl), left),
      DirichletBC(V, Constant(Tr), right)]
```

## Loading the material behaviour

We use the `MFrontNonlinearMaterial` class for describing the material behaviour. The first argument corresponds to the path where material librairies have been compiled, the second correspond to the name of the behaviour (declared with `@Behaviour`). Finally, the modelling hypothesis is specified (default behaviour is `"3d"`).


```python
material = mf.MFrontNonlinearMaterial("./src/libBehaviour.so",
                                      "StationaryHeatTransfer",
                                      hypothesis="plane_strain")
```

The `MFront` behaviour declares the field `"TemperatureGradient"` as a Gradient variable, with its associated Flux called `"HeatFlux"`. We can check that the `material` object retrieves `MFront`'s gradient and flux names, as well as the different tangent operator blocks which have been defined, namely `dj_ddgT` and `dj_ddT` in the present case:


```python
print(material.get_gradient_names())
print(material.get_flux_names())
print(["d{}_d{}".format(*t) for t in material.get_tangent_block_names()])
```

    ['TemperatureGradient']
    ['HeatFlux']
    ['dHeatFlux_dTemperatureGradient', 'dHeatFlux_dTemperature']


## Non-linear problem definition

When defining the non-linear problem, we will specify the boundary conditions and the requested quadrature degree which will control the number of quadrature points used in each cell to compute the non-linear constitutive law. Here, we specify a quadrature of degree 2 (i.e. 3 Gauss points for a triangular element). 


```python
problem = mf.MFrontNonlinearProblem(T, material, quadrature_degree=2, bcs=bc)
```

## Variable registration

The `MFront` behaviour implicitly declares the temperature as an external state variable called `"Temperature"`. We must therefore associate this external state variable to a known mechanical field. This can be achieved explicitly using the `register_external_state_variable` method. In the present case, this can be done automatically since the name of the unknown temperature field matches the [TFEL Glossary](http://tfel.sourceforge.net/glossary.html) name `"Temperature"` which is a predefined external state variable. In this case, the following message will be printed:
```
Automatic registration of 'Temperature' as an external state variable.
```
For problems in which the temperature only acts as a parameter (no jacobian blocks with respect to the temperature), the temperature can be automatically registered as a constant value ($293.15 \text{ K}$ by default) or to any other (`dolfin.Constant`, `float` or `dolfin.Function`) value using the `register_external_state_variable` method.

In the `FEniCS` interface, we instantiate the main mechanical unknown, here the temperature field `T` which has to be named `"Temperature"` in order to match `MFront`'s predefined name. Using another name than this will later result in an error saying:
```
ValueError: 'Temperature' could not be associated with a registered gradient or a known state variable.
```
Finally, we need to associate to `MFront` gradient object the corresponding UFL expression as a function of the unknown field `T`. To do so, we use the `register_gradient` method linking `MFront` `"TemperatureGradient"` object to the UFL expression `grad(T)`.


```python
problem.register_gradient("TemperatureGradient", grad(T))
```

Similarly to the case of external state variables, common gradient expressions for some [TFEL Glossary](http://tfel.sourceforge.net/glossary.html) names have been already predefined which avoid calling explicitly the `register_gradient` method. Predefined expressions can be obtained from:


```python
mf.list_predefined_gradients()
```

    'TFEL gradient name'   (Available hypotheses)
    ---------------------------------------------
    'Strain'               ('Tridimensional', 'PlaneStrain', 'Axisymmetrical')
    'TemperatureGradient'  ('Tridimensional', 'PlaneStrain', 'Axisymmetrical')
    'DeformationGradient'  ('Tridimensional', 'PlaneStrain', 'Axisymmetrical')
    


We can see that the name `"Temperature Gradient"` is in fact a predefined gradient. Omitting calling the `register_gradient` method will in this case print the following message upon calling `solve`:
```
Automatic registration of 'TemperatureGradient' as grad(Temperature).
```
meaning that a predefined gradient name has been found and registered as the UFL expression $\nabla T$. 
> Note that automatic registration is not supported when using mixed function spaces.

## Problem resolution
No external loading has been specified so that the residual form is automatically defined as:
\begin{equation}
F(\widehat{T}) = \int_\Omega \mathbf{j}\cdot \nabla \widehat{T} \text{dx} = 0 \quad \forall \widehat{T}
\end{equation}

From the two tangent operator blocks `dj_ddgT` and `dj_ddT`, it will automatically be deduced that the heat flux $\mathbf{j}$ is a function of both the temperature gradient $\mathbf{g}=\nabla T$ and the temperature itself i.e. $\mathbf{j}=\mathbf{j}(\mathbf{g}, T)$. The following tangent bilinear form will therefore be used  when solving the above non-linear problem:

\begin{equation}
a_\text{tangent}(\widehat{T},T^*) = \int_{\Omega} \nabla \widehat{T}\cdot\left(\dfrac{\partial \mathbf{j}}{\partial \mathbf{g}}\cdot \nabla T^*+\dfrac{\partial \mathbf{j}}{\partial T}\cdot T^*\right) \text{dx}
\end{equation}

We finally solve the non-linear problem using a default Newton non-linear solver. The `solve` method returns the number of Newton iterations (4 in the present case) and converged status .


```python
problem.solve(T.vector())
```

    Automatic registration of 'Temperature' as an external state variable.
    
    (4, True)



We finally check that the thermal conductivity coefficient $k$, computed from the ratio between the horizontal heat flux and temperature gradient matches the temperature-dependent expressions implemented in the `MFront` behaviour.


```python
j = problem.fluxes["HeatFlux"].function
g = problem.gradients["TemperatureGradient"].function
k_gauss = -j.vector().get_local()[::2]/g.vector().get_local()[::2]
T_gauss = problem.state_variables["external"]["Temperature"].function.vector().get_local()
A = material.get_parameter("A");
B = material.get_parameter("B");
k_ref = 1/(A + B*T_gauss)
plt.plot(T_gauss, k_gauss, 'o', label="FE")
plt.plot(T_gauss, k_ref, '.', label="ref")
plt.xlabel(r"Temperature $T\: (K)$")
plt.ylabel(r"Thermal conductivity $k\: (W.m^{-1}.K^{-1})$")
plt.legend()
plt.show()
```


    <IPython.core.display.Javascript object>



<img src="data:image/png;base64,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" width="640">

