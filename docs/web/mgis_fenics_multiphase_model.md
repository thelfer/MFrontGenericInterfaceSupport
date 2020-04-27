% Multiphase model for fiber-reinforced materials
% Jérémy Bleyer
% 17/04/2020

\newcommand{\bsig}{\boldsymbol{\sigma}}
\newcommand{\beps}{\boldsymbol{\varepsilon}}
\newcommand{\bkappa}{\boldsymbol{\kappa}}
\newcommand{\DD}{\mathbb{D}}
\newcommand{\CC}{\mathbb{C}}
\newcommand{\bI}{\boldsymbol{I}}
\newcommand{\bV}{\boldsymbol{V}}
\newcommand{\bU}{\boldsymbol{U}}
\newcommand{\bt}{\boldsymbol{n}}
\newcommand{\bn}{\boldsymbol{t}}
\newcommand{\T}{{}^\text{T}} 
\newcommand{\div}{\operatorname{div}}
\newcommand{\jump}[1]{[\![#1]\!]}
\newcommand{\avg}[1]{\left\langle#1\right\rangle}

This demo explores the full capability of `MFront`'s recent extension to
handling generalized behaviours, namely involving different pairs of
fluxes and gradients, as well as `FEniCS` versatility to solve generic
PDEs.

This will be illustrated on a generalized continuum model called
*multiphase model* in the context of fiber-reinforced materials [@bleyer2018multiphase].

> **Source files:**
>
> * Jupyter notebook: [mgis_fenics_multiphase_model.ipynb](https://gitlab.enpc.fr/navier-fenics/mgis-fenics-demos/raw/master/demos/multiphase_model/mgis_fenics_multiphase_model.ipynb)
> * Python file: [mgis_fenics_multiphase_model.py](https://gitlab.enpc.fr/navier-fenics/mgis-fenics-demos/raw/master/demos/multiphase_model/mgis_fenics_multiphase_model.py)
> * MFront behaviour file: [MultiphaseModel.mfront](https://gitlab.enpc.fr/navier-fenics/mgis-fenics-demos/raw/master/demos/multiphase_model/MultiphaseModel.mfront)


## A quick primer on the multiphase model

![](img/multiphase_kinematics.svg ""){width=75%}

The multiphase model is a higher-grade (i.e. with enhanced kinematics)
generalized continuum which represents a biphasic material (the main
application being fiber-reinforced materials) by a superposition of two
phases (say *matrix* and *fiber* phases), each possessing their own
kinematics and in mutual interaction. Each phase kinematics is described
by a displacement field $\bU^1$ and $\bU^2$ for the matrix and fiber
phase respectively.

In the present case, each phase is a standard Cauchy continuum but other
variants, including fiber bending effects for instance, exist, see
@debuhan2017elastic for more details.

* **Generalized strains**: infinitesimal strain in each phase
  $\beps^j=\nabla^s \bU^j$ and relative displacement between both phases
  $\bV = \bU^2-\bU^1$

* **Generalized stresses**: partial stress $\bsig^j$ in both phases and
  interaction force $\bI$ (exerted by phase 2 over phase 1)

* **Equilibrium equations**:

\begin{align}
\div \bsig^1 +\rho_1F + \bI &=  0 \label{multiphase-eq1}\\
\div \bsig^2 + \rho_2F - \bI &= 0 \label{multiphase-eq2}
\end{align}

* **Traction boundary conditions**:

\begin{align}
\bsig^1\cdot\bn &= \bt^1 \text{ on }\partial \Omega_{T} \\
\bsig^2\cdot\bn &= \bt^2 \text{ on } \partial \Omega_{T} 
\end{align}

* **Deformation work density**:

\begin{equation}
w_{\text{def}} = \bsig^1:\beps^1 + \bsig^2:\beps^2+ \bI\cdot\bV
\end{equation}

* General form of the **constitutive equations** based on a convex free
  energy density $\psi$:

\begin{align}
\bsig^1 &= \dfrac{\partial \psi}{\partial \beps^1}(\beps^1,\beps^2,\bV) \\
\bsig^2 &= \dfrac{\partial \psi}{\partial \beps^2}(\beps^1,\beps^2,\bV)\\
\bI &= \dfrac{\partial \psi}{\partial \bV}(\beps^1,\beps^2,\bV)
\end{align}

which will be particularized for the present demo to the following
*elastic behaviour*:

\begin{align*}
\bsig^1 &= \DD^{11}:\beps^1+\DD^{12}:\beps^2\\
\bsig^2 &= (\DD^{12})\T:\beps^1+\DD^{22}:\beps^2\\
\bI &= \bkappa \cdot \bV
\end{align*}
in which $\DD^{ij}$ are partial stiffness tensors and $\kappa$ can be
seen as an interaction stiffness between both phases.

## `MFront` implementation

The `MFront` implementation expands upon the detailed `MFront`
implementation of the [stationnary heat transfer
demo](nonlinear_heat_transfer.html). Again, the
`DefaultGenericBehaviour` is used here and specify that the following
implementation will only handle the 2D plane strain case.

``` cpp
@DSL DefaultGenericBehaviour;
@Behaviour MultiphaseModel;
@Author Jeremy Bleyer;
@Date 04/04/2020;

@ModellingHypotheses {PlaneStrain};
```

## Flux/Gradient dual pairs
The three pairs of generalized flux/gradient of the multiphase model are
then defined, namely $(\bsig^1,\beps^1)$, $(\bsig^2,\beps^2)$ and
$(\bI,\bV)$:

```cpp
@Gradient StrainStensor e₁;
e₁.setEntryName("MatrixStrain");
@Flux StressStensor σ₁;
σ₁.setEntryName("MatrixStress");

@Gradient StrainStensor e₂;
e₂.setEntryName("FiberStrain");
@Flux StressStensor σ₂;
σ₂.setEntryName("FiberStress");

@Gradient TVector V;
V.setEntryName("RelativeDisplacement");
@Flux TVector I;
I.setEntryName("InteractionForce");
```

We now declare the various tangent operator blocks which will be needed
to express the generalized tangent operator. Of the nine possible
blocks, only 5 are needed here since the partial stresses do not depend
on the relative displacement and, similarly, the interaction force does
not depend on the phase strains:

```cpp
@TangentOperatorBlocks{∂σ₁∕∂Δe₁,∂σ₁∕∂Δe₂,∂σ₂∕∂Δe₁,∂σ₂∕∂Δe₂,∂I∕∂ΔV};
```

### Defining material properties

We consider the case of a 2D plane strain bilayered material made of
isotropic materials for both the matrix and fiber phases. In addition to
the four material constants need for both constituents, the volume
fraction of both phases (here we ask for the volume fraction $\rho$ of
the fiber phase) and the size $s$ of the material unit cell (i.e. the
spacing between two consecutive matrix layers for instance) are the two
other parameters characterizing the multiphase elastic model.

```cpp
@MaterialProperty stress Y1;
Y1.setEntryName("MatrixYoungModulus");
@MaterialProperty real nu1;
nu1.setEntryName("MatrixPoissonRatio");
@MaterialProperty stress Y2;
Y2.setEntryName("FiberYoungModulus");
@MaterialProperty real nu2;
nu2.setEntryName("FiberPoissonRatio");
@MaterialProperty real ρ;
ρ.setEntryName("FiberVolumeFraction");
@MaterialProperty real s;
s.setEntryName("Size");
```

## Generalized elastic behaviour

It has been shown in @bleyer2018multiphase that, for materials made of
two constituents, the partial stiffness tensors $\DD^{ij}$ can be
expressed as functions of the material individual stiffness
$\CC^1,\CC^2$ and the macroscopic stiffness $\CC^{hom}$ of the composite
(obtained from classical homogenization theory) as follows:

\begin{align*}
\DD^{11} &= \phi_1\CC^1-\CC^1:\jump{\CC}^{-1}:\Delta\CC:\jump{\CC}^{-1}:\CC^1\\ 
\DD^{22} &= \phi_2\CC^2-\CC^2:\jump{\CC}^{-1}:\Delta\CC:\jump{\CC}^{-1}:\CC^2\\ 
\DD^{12} &= \CC^1:\jump{\CC}^{-1}:\Delta\CC:\jump{\CC}^{-1}:\CC^2
\end{align*}

with $\phi_1=1-\rho$, $\phi_2 = \rho$ are the phases volume fractions,
$\jump{\CC}=\CC^2-\CC^1$ and $\Delta\CC =
\phi_1\CC^1+\phi_2\CC^2-\CC^{hom}$. They satisfy the following property
$\DD^{11}+\DD^{22}+\DD^{12}+(\DD^{21})\T=\CC^{hom}$.

In the present case of a 2D bilayered material made of isotropic
constituents, the macroscopic stiffness is given by @debuhan2017elastic:

\begin{equation}
\CC^{hom} = \begin{bmatrix}
\avg{E_\text{oe}} - \avg{\lambda^2/E_\text{oe}}+\avg{\lambda/E_\text{oe}}^2/\avg{1/E_\text{oe}} & \avg{\lambda/E_\text{oe}}/\avg{1/E_\text{oe}} & \avg{\lambda/E_\text{oe}}/\avg{1/E_\text{oe}} & 0 \\
\avg{\lambda/E_\text{oe}}/\avg{1/E_\text{oe}} & 1/\avg{1/E_\text{oe}} & \avg{\lambda} & 0 \\
\avg{\lambda/E_\text{oe}}/\avg{1/E_\text{oe}} & \avg{\lambda} & \avg{E_\text{oe}} & 0 \\
0 & 0 & 0 & 2/\avg{1/\mu}
\end{bmatrix}_{(xx,yy,zz,xy)}
\end{equation}

where $\avg{\star}=\phi_1\star_1+\phi_2\star_2$ denotes the averaging
operator and $E_\text{oe}=\lambda+2\mu$ is the oedometric modulus.

As regards the interaction stiffness tensor, it is obtained from the
resolution of an auxiliary homogenization problem formulated on the
heterogeneous material unit cell. In the present case, one obtains
@bleyer2018multiphase:

\begin{equation}
\bkappa = \begin{bmatrix}
\dfrac{12}{\avg{\mu}s^2} & 0 \\
0 & \dfrac{12}{\avg{E_\text{oe}}s^2}
\end{bmatrix}_{(x,y)}
\end{equation}

These relations are all implemented in the behaviour `@Integrator` which
also defines the required tangent blocks. `MFront` linear algebra on
fourth-order tensor makes the implementation very easy and Unicode
support provides much more readable code.

```cpp
@ProvidesTangentOperator;
@Integrator {
  // remove useless warnings, as we always compute the tangent operator
  static_cast<void>(computeTangentOperator_);

  const auto λ₁ = computeLambda(Y1,nu1);
  const auto μ₁ = computeMu(Y1,nu1);
  const auto λ₂ = computeLambda(Y2,nu2);
  const auto μ₂ = computeMu(Y2,nu2);
  const auto Eₒₑ¹ = λ₁+2*μ₁;
  const auto Eₒₑ² = λ₂+2*μ₂;
  const auto Eₒₑ = (1-ρ)*Eₒₑ¹ + ρ*Eₒₑ²;
  const auto iEₒₑ = (1-ρ)/Eₒₑ¹ + ρ/Eₒₑ²;
  const auto λ = (1-ρ)*λ₁ + ρ*λ₂;
  const auto λiEₒₑ = (1-ρ)*λ₁/Eₒₑ¹ + ρ*λ₂/Eₒₑ²;
  const auto λ2iEₒₑ = (1-ρ)*λ₁*λ₁/Eₒₑ¹ + ρ*λ₂*λ₂/Eₒₑ²;
  const auto iμ = (1-ρ)/μ₁ + ρ/μ₂;
  const auto C₁₁₁₁ʰᵒᵐ = Eₒₑ - λ2iEₒₑ + λiEₒₑ*λiEₒₑ/iEₒₑ;
  const auto C₁₁₂₂ʰᵒᵐ = λiEₒₑ/iEₒₑ;
  const Stensor4 Cʰᵒᵐ = {
      C₁₁₁₁ʰᵒᵐ , C₁₁₂₂ʰᵒᵐ, C₁₁₂₂ʰᵒᵐ, 0., //
      C₁₁₂₂ʰᵒᵐ, 1/iEₒₑ, λ, 0., //
      C₁₁₂₂ʰᵒᵐ, λ, Eₒₑ, 0., //
      0., 0., 0., 2/iμ
  };
  const auto C¹ = λ₁ ⋅ (I₂ ⊗ I₂) + 2 ⋅ μ₁ ⋅ I₄;
  const auto C² = λ₂ ⋅ (I₂⊗I₂) + 2 ⋅ μ₂ ⋅ I₄;
  const auto iΔC = invert(C²-C¹);
  const auto Cᵃᵛᵍ = (1-ρ)⋅C¹ + ρ⋅C²;
  const auto H = iΔC ⋅ (Cᵃᵛᵍ-Cʰᵒᵐ) ⋅ iΔC;
  const auto D¹¹ = (1-ρ)⋅C¹-(C¹ ⋅ H ⋅ C¹);
  const auto D²² = ρ⋅C²-(C² ⋅ H ⋅ C²);
  const auto D¹² = C¹ ⋅ H ⋅ C²;
  ∂σ₁∕∂Δe₁ = D¹¹;
  ∂σ₁∕∂Δe₂ = D¹²;
  ∂σ₂∕∂Δe₁ = transpose(D¹²);
  ∂σ₂∕∂Δe₂ = D²²;
  σ₁ = ∂σ₁∕∂Δe₁ ⋅ (e₁ + Δe₁) + ∂σ₁∕∂Δe₂ ⋅ (e₂ + Δe₂);
  σ₂ = ∂σ₂∕∂Δe₁ ⋅ (e₁ + Δe₁) + ∂σ₂∕∂Δe₂ ⋅ (e₂ + Δe₂);

  const auto κʰ = 12/iμ/s/s; // horizontal interaction stiffness
  const auto κᵛ = 12/iEₒₑ/s/s; // vertical interaction stiffness
  const tmatrix<N, N, real> κ = {
      κʰ, 0., 
      0., κᵛ
  };
  ∂I∕∂ΔV = κ;
  I = κ ⋅ (V + ΔV);
}
```

# `FEniCS` implementation

## Geometry and generalized boundary conditions

We consider a rectangular domain with a standard $P_1$ Lagrange
interpolation for both displacement fields $\bU^1,\bU^2$. The total
mixed function space `V` is therefore made of two vectorial $P_1$
Lagrange elementary spaces. The same boundary conditions will be applied
for both phases: fixed vertical displacement on the bottom boundary,
imposed vertical displacement on the top boundary and fixed horizontal
displacement on the left boundary. The considered problem corresponds to
the pure compression of a multilayered block, an example treated in the
previously mentioned references.

In such a problem, both phase strains are uniform along the vertical
direction solution but vary along the horizontal direction due to a
boundary layer effect appearing near the right free-edge. Such a feature
cannot be captured by a classical Cauchy continuum but is well
reproduced by the multiphase model.

```python
%matplotlib notebook
import matplotlib.pyplot as plt
from dolfin import *
import mgis.fenics as mf
import numpy as np
from ufl import diag

width = 0.5
height = 1.
mesh = RectangleMesh(Point(0., 0.), Point(width, height), 100, 100)

Ve = VectorElement("CG", mesh.ufl_cell(), 1)
V = FunctionSpace(mesh, MixedElement([Ve, Ve]))
u = Function(V, name="Displacements")
(u1, u2) = split(u)

def bottom(x, on_boundary):
    return near(x[1], 0) and on_boundary
def top(x, on_boundary):
    return near(x[1], height) and on_boundary
def left(x, on_boundary):
    return near(x[0], 0) and on_boundary


bc = [DirichletBC(V.sub(0).sub(0), Constant(0), left),
      DirichletBC(V.sub(1).sub(0), Constant(0), left),
      DirichletBC(V.sub(0).sub(1), Constant(0), bottom),
      DirichletBC(V.sub(1).sub(1), Constant(0), bottom),
      DirichletBC(V.sub(0).sub(1), Constant(-1), top),
      DirichletBC(V.sub(1).sub(1), Constant(-1), top)]

facets = MeshFunction("size_t", mesh, 1)
ds = Measure("ds", subdomain_data=facets)
```

## Generalized gradients registration

After having defined the `MFrontNonlinearMaterial` and
`MFrontNonlinearProblem` instances, one must register the UFL expression
of the gradients defined in the `MFront` file since, in this case,
automatic registration is not available for this specific model. The
matrix (resp. fiber) strains $\beps^1$ (resp. $\beps^2$) is simply given
by `sym(grad(u1))` (resp. `sym(grad(u2))`) whereas the relative
displacement $\bV$ is the vector `u2-u1`:


```python
mat_prop = {"MatrixYoungModulus": 10.,
            "MatrixPoissonRatio": 0.45,
            "FiberYoungModulus": 10000,
            "FiberPoissonRatio": 0.3,
            "FiberVolumeFraction": 0.01,
            "Size": 1/20.}

material = mf.MFrontNonlinearMaterial("./src/libBehaviour.so",
                                      "MultiphaseModel",
                                      hypothesis="plane_strain",
                                      material_properties=mat_prop)

problem = mf.MFrontNonlinearProblem(u, material, quadrature_degree=0, bcs=bc)
problem.register_gradient("MatrixStrain", sym(grad(u1)))
problem.register_gradient("FiberStrain", sym(grad(u2)))
problem.register_gradient("RelativeDisplacement", u2-u1)

problem.solve(u.vector())
```

    Automatic registration of 'Temperature' as a Constant value = 293.15.
    
    (1, True)



## Comparison with an analytical solution

We compare the horizontal displacements in both phases with respect to
known analytical solutions for this problem in the case of very stiff
inclusions $E_2 \gg E_1$ in small volume fraction $\eta \ll 1$ (see
@debuhan2017elastic). Note that the more complete solution for the
general case can be found in @bleyer2018multiphase.

```python
x = np.linspace(0, width, 21)
plt.figure()
plt.plot(x, np.array([u(xi, height/2)[0] for xi in x]), "ob", label=r"$u_x^\textrm{matrix}$ (FE)")
plt.plot(x, np.array([u(xi, height/2)[2] for xi in x]), "or", label=r"$u_x^\textrm{fiber}$ (FE)")

E1 = mat_prop["MatrixYoungModulus"]
nu1 = mat_prop["MatrixPoissonRatio"]
E2 = mat_prop["FiberYoungModulus"]
nu2 = mat_prop["FiberPoissonRatio"]
eta = mat_prop["FiberVolumeFraction"]
s = mat_prop["Size"]
mu1 = E1/2/(1+nu1)
mu2 = E2/2/(1+nu2)
lmbda1 = E1*nu1/(1+nu1)/(1-2*nu1)
lmbda2 = E2*nu2/(1+nu2)/(1-2*nu2)
kappa = 12*mu1*mu2/((1-eta)*mu2+eta*mu1)/s**2 # horizontal interaction stiffness
Eo1 = lmbda1+2*mu1 # oedometric modulus for phase 1
Eo2 = lmbda2+2*mu2 # oedometric modulus for phase 1
Er = eta*E2/(1-nu2**2) 
ahom = Eo1+Er
l = sqrt(Eo1*Er/ahom/kappa) # material characteristic length scale
um = lmbda1/ahom*(x+l*(Er/Eo1)*np.sinh(x/l)/np.cosh(width/l))
ur = lmbda1/ahom*(x-l*np.sinh(x/l)/np.cosh(width/l))
plt.plot(x, um, '--b', label=r"$u_x^\textrm{matrix}$ (ref.)")
plt.plot(x, ur, '--r', label=r"$u_x^\textrm{fiber}$ (ref.)")
plt.legend()
plt.show()
```

<img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAwAAAAJACAYAAAA6rgFWAAAgAElEQVR4nOzde5QcV30v+i0/sQ24DMGagAnQwElMAuZUQpIbVswJ6pBgLbgX3JfkkOSYgN08DjHnBJ8OAZ8IRgHRYREnDsTuQG6MwwnYnYFBPGTSyGNkjF7TeoyERpamHcmxNRYauz3vmZru/t4/Srunu6v6UbsfVbXr+1mrlq2uftTMb2pmf6v2Q4CIiIiIiCJD+H0AREREREQ0OAwAREREREQRwgBARERERBQhDABERERERBHCAEBEREREFCEMAEREREREEcIAQEREREQUIQwAREREREQRwgBARERERBQhDABERERERBHCAEBEREREFCEMAEREREREEcIAQEREREQUIQwAREREREQRwgBARERERBQhDABERERERBHCAEBEREREFCEMAEREREREEcIAQEREREQUIQwAREREREQRwgBARERERBQhDABERERERBHCAEBEREREFCEMAEREREREEcIAQEREREQUIQwAREREREQRwgBARERERBQhDABERERERBHCAEBEREREFCEMAEREREREEcIAQEREREQUIQwAREREREQRwgBARERERBQhDABERERERBHCAEBEREREFCEMAEREREREEcIAQEREREQUIQwAREREREQRwgBARERERBQhDABERERERBHCAEBEREREFCEMAEREREREEcIAQEREREQUIQwAREREREQRwgBARERERBQhDABERERERBHCAEBEREREFCEMAEREREREEcIAQEREREQUIQwAREREREQRwgBARERERBQhDABERERERBHCAEBEREREFCEMAEREREREEcIAQEREREQUIQwAREREREQRwgBARERERBQhDABERERERBHCAEBEREREFCGhDgC5XA7JZBLpdLr6X1X5fB6JRMLTa9LpNLLZrPJnEhERERENWmgDQDabhWmadY+l02nPjXgAyGQyEEIgHo93/JpisQghBDKZjOfPIyIiIiLyS2gDgGEYyOVyro93elU+kUggHo8jk8nAMAxPASCZTDIAEBEREVHohDIAZLNZCOF+6IlEQukugJcAkMvlkE6nGQCIiIiIKHRCGQCSySQMw3Ddl0qlmu5rxUsASKVSKBQKDABEREREFDqhDACxWMzR/1+S/fnz+byn9+w0AGQyGRQKBQYA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H/sAeN5XqVT6tq/T79sDDzwAwzCq+06ePImZmZmOv9+lUglLS0uuj+tYJ9M0sWnTpoHXqZN9hmEgkUh4ft3U1BQKhQIsy6o7D8+dO4eJiQmsrKxgbW2tbt/s7CwWFxerX1/jvoWFhernNe6bn5+v1qlx39zcXLUWjftqGzRzc3OOffJ7Oj8/79gnv/6FhQXHPln/xcVFx761tTUAwNLSkmOfdb61sby83HTfysqKY1+z33uzs7NYWVkBAEctZmdnsby83HTf0tISAHRcp3PnzuHQoUN45plnWKcA10nlfJK1PXfunFZ1+sxn7Mb/z/5sGXv2rB9TWOukcj796Ec/YgBQwABwXr8DwJYtWyCEcN3uvPNOjI6OVjfZyJ+amqp7fNeuXThx4gQA+xdzsVis2+QJaFmWY588AdfW1hz75AlYKpUc++QJWC6XHfvkCVipVBz7nn322erX3rivth/2s88+69gnf8HOzs469slfsHNzc4598hfs/Py8Y5/8BbuwsODYJ3/Bfvvb38Z1112HYrGIm266Cdddd131F+Wf/umfOl4nf1EuLy879slflFGp06lTp3DllVfijjvu6HudFhcXHftknZaWluoev/XWW6uzA3mt06OPPopCoYDp6em683B0dBTj4+MAgJmZGce+3bt3V78vjft27dpVPc7GfTt37qzWt3Hfjh07qrXYvn173b7t27dX9+3YscPxWvn93rlzp2Of/Pp37drl2Cd/bnbv3u3YNzMzAwAYHx937JuengYATExMOPadPn0aADA5OenY1+z33ujoKCYnJwEAp0+fduybmJgAANaJdYpknSoV4KMfPYtM5vuRrdPdd9/NAKBA2wAQtC5Are4APPLIIx1dYeEdgP5eWd6/fz9uvPFG3Hfffdi3bx9uueUW3HXXXUin03j66ac9fb9LpRJWVlZcH9e1Tvv27YMQAvv27QvEHYD77rsPhmHg5PlR7r26A/D000/j5MmTsCwrcFfCwnDFEgjuleWnn34ax48frwZj1imYdVI5n2Rtn3766dDW6Zln9K+Tyvn0yCOPMAAoCGUAMAyj7aw8vR4E3OvPAzgGQGelUnTGAOiKYwCih7XVV9hru3s3cO21zWftizKOAVATygBgmmbT1UflIk9eG8qtAkA/Pg9gANAZA0D4MQBED2urrzDX9oEHgMsvt/v6v+IVwJkzfh9RgFgW9n3xiwwACkIZAOTiQ25aLdrVSqsA0I/PAxgAdMYAEH4MANHD2uorrLX92tfsqT1rZ/v5y7/0+6gCwLKA4WFg40YcPT+ekgHAm1AGALm6qNsf51gs5jrfd7sFgFoFAJXP6wQDgL4YAMKPASB6WFt9hbG2X/jC+iJecvv4x+2VfyPNsoAbbqiucsYAoCaUAQCwV2ZNpVJ1j8lVSRsb+6Zptl0FVAjRcjEjL5/XKRkAOp0/nQEgPMrlMhYXF6uDwih8mp1vy8vLmJiYqA6WI32wtvoKU20rFWDLFueKvnfc4feRBcTwcN03hgFATWgDQLFYrFupNZ/PwzRN18G4yWQSsVjM0VBPpVKIx+MwDKM6JWcsFkM8HnfM6uPl8zolA0CnP7QMAESDw/ONiAatVAI+9KH6hv+FFwL33uv3kQWEZQEbN9bdGmEAUBPaACDlcjmk02lks1nlK/F+fZ4MAAcPHuzo+WyQhEe5XIZlWbwDEGLNzjfLsjA9PV2dgo/0wdrqKwy1XV0Ffu/36hv/l10GfOc7fh9ZBwa1/PDYmOPWCAOAmtAHgDDjGAB9cQxA+HEMQPSwtvoKem3n54Hf/u36tq1hAD/6kd9H1kbNYNy6gx8ash/vdRAYGWEA6BEGAB8xAOiLASD8GACih7XVV9Br++1v17drf/ZngfOL8gZXw2Dcui9A/nvz5u5DwDPPAPv2Af/yL/YoaAaAnmAA8BEDQPDl83lkMhmk02nXwdrFYhG5XA7ZbLY6PgRgANABA0D0sLb6CkNt77jDbtO+8pXAY4/5fTQdaBiM23TburX1+1QqwE9/Cvz4x/Zgh7/4C+Dd7wZ+9VeBF7yg/r1e9SqOAegRBgAfMQAEV6FQgGmaSKVSKBQKKBaLyGazjpmgCoUCEolEdQC5xAAQfgwA0cPa6isstf3iF4GnnvL7KDrgMhjXdduwwe4OtLpqj3Bu9E//BDz/+Z0FCSGAiy5yTJHEAKCGAcBHMgCMj4939HwGgMEoFoswDKPuij6wvuqzGznTlFQul7GwsMBBwCHW7HxbWlrC+Pg4lpaWfDgq6ifWVl9Bq+3qqt9H0CWXwbgtt+c8x+7n1MilT3/L7aUvBU6dsrsWCa4D0A0GAB9xGtBgymazrg39VCrVdNXnxgBA4cfzjYi60WxinH/7N+BlLwOOHfPx4LrlteEumixkcPhw/XMuuACIxYC3vMWeD/Wv/xrYvh34yU+A2vBmWXbXoqEhBgBFDAA+kgHg0KFDHT2fDZLBUGnMN76mUqlgbW0Nlcgv2Rhezc63tbU1zMzMYG1tzYejon5ibfU1yNq2mhjnXe8CLr7Y/vc11wCPP973w+mPH/zAewD47//d+T6Li8Df/R2wYwdw4oT3WyOWhcNf/jIDgAIGAB9xDEAw9SIAcAxA+HEMQPSwtvoaVG1bTYzjtv31X/f1cPrHsuz++J00/C+/HPjCFwCXiTR6Yc+ePQwAChgAfMQAEEwMAAQwAEQRa6uvQdW204lxhAA+97m+Hkr3ymXg4YeBVApwu3Pylrf0ZhagLjEAqGEA8BEDQLDkcjnE43EYhgEhBOLxeHUzTRPifD9Dt+lAawNAOp3Gtm3bcOutt+KWW25pWbN0Oo1UKoV0Oo1kMlk38DiXy8E0TRiGgXg8DsAeh5BIJJBIJHr81VMjBoDoYW31NYjadjoxjhDAlVf2b7HcrlblLZWAH/4Q+PCH7cUI5AHv3Ol87qOPAi95SfPZf4TozToAbTAAqGEA8BEDQDA1uwOQyWRaBgDDMJBMJgGs3wE4ceIEDMNALpdzvCYejyOdTtc9JqceleR0pKZpIp1Oo1gsVmcjKhaL3X6p1AIDQPSwtvoaRG29TowzNtbjA1BdlbdUAh56yO6jPzTkfrAf+EDzzzw/GNfxmVu39r3xDzAAqGIA8JEMAPv37+/o+QwAg9EsAORyuZYBoLZRXiqVMD8/j1Kp5Dp7UDqdbvkZtXWWx5PNZgGguiaBqlZfhx+y2SwMwwjcz3az821xcRG7d+/G4uKiD0dF/cTa6msQtfU6Mc7ISA8/3OuqvKUS8OCDwAc/6AwMjZthALfd1v7zVe86dGn//v0MAAoYAHzEaUCDSTUANJsiNJ/PQwhRd7XfMIym3XiEEHUN/MZw0Y1CoeB4/34oFovV74nsTtWKfG6Q7mrwfCMiL3y9A+B1Vd53vKP1817wAuC977Vn5wn4ogVe21JkYwDwkfyhPXz4cEfPZ4NkMHoRACqVCkqlEiqVCorFIoQQ1Qa/bITLLkCNWzwer+syJANAL8RiMZim2ZP3avc5sktUKpWqdo1q95ogjW1odr6VSiXMzs5ygLeGWFt9DaK2cgxAuza4XBy3ZxfJva7Ka1n2ksON+1/4QuDmm4Hvf3+gV/C7dfjwYQYABQwAPuIYgGDqRQBonAVICFFteMv3qe3r3+54mt1d8EIucNbvrj9yrITXn1V5fEH5GecYgOhhbfUVtFmAejoxjsqth+lpOxD8zM8AySSQy4Wq0V+LYwDUMAD4iAEgmPoVAORV8MY7AqrH45VpmgNZrTiRSCjfsagdSO03BoDoYW31Nch1ADZvHvDEOKqDD/bvd5/eM2QYANQwAPiIASCYeh0AZJef2ik+DcNo2RWn9jN6EQDkOIRO7zp0Ix6PKweAXo536BYDQPSwtvrqR23X1oBPfAL40pfqHx/4xDj33uvj4AP/MQCoYQDwEQNAMPVqFiAZAFKplOP9WnWTSafTPQ8A6XR6IN1/gO4CgPy+uE2bOmgMANHD2uqr17U9cwZ405vs9vSllwKHDjmf09eJcSoVu6++nPmnk63ngw+CgQFADQOAj2QA2LdvX0fPZwAYjGZ97tsFgEQiUb3CXiqVMDc3hwceeACGYbi+JpFIwDTNuqvd+XzesTZALwYBt2uU196lkGsNxGKxursWgB0k5KJopmnWzSYkG++Nm5e7DoO8U9FOs/NtYWEBu3btwsLCgg9HRf3E2uqrl7X9wQ+Aq6+ub1u/7nX2wrkDMT4OvOY13q7692XwQTDs27ePAUABA4CPOA1o/9XOrNMon8/XNfTz+Xzdir9y8S0AdSsE167MKyWTSRQKhWrjWc58k0wmW3ZnyWaz1eCQTqfrGtyNxxOLxRCPx5W6xxiG0XIgsQwAjQ382ivxsm9/IpGofk9ruzYVCgVkMhnEYrHq45lMxvNdh06mDR0Enm9EVKtUArZscU6286pXAQcPDvBAnnoKuOSS+oO45BLgmmvqBxv0ffBBMHAaUDUMAD6SP7QTExMdPZ8NEm/kyrmA3QBuvJrttkBXr1QqFZTLZVQqlb68v1e1sxC5kQFABo3GRru8ut+4foDb3ZJuugABdq0GMVi5nWbnW7lcxtLSEsoDu9xHg8La6qvb2k5PA29+s/OC+rveBfStx1ilApw65b7vppvsA9i4EfjUp+xQEIBVef0wMTHBAKCAAcBHHAPQX7LBL7uVNPYr7+e8842zAPlJzjrU6qp6bQBw+xmTdx+ava72e9ttAJBrCPiNYwCih7XVVze1ffBBZ5v6kkvsqfT7co1nddUe2PvLvww873nuCePoUfs5KyvOfT6uyusHjgFQwwDgIwaA/pHdcQD3K/2yUdx4V6BXghQAZCO91fSatYuTNZLfq1Zb7fex2wAguyD5jQEgelhbfanUtly2L5xfcEF94z8Ws7vht+W1IX72rH01vzFt/O3fdnzMUcQAoMb/v7IRxgAwGLFYzNH47feiU2ENAG7PkXdQkskk8vm861Y7LoEBgMKKtdWXSm0/9jFnl593vhN49tk2L7Qse0WwxmWBh4bsxxuDwMGDwHve4+zXL7c3vMH7FxwhDABq/P8rG2EMAP3XrPtPr1bXbSZIAaCThcdkAGg2+067AFGrF2MA2AWI/MDa6kultk88YS+UKwRw8cX2hfi2XX4sa31qzlaDcZeXgW98A7j++uYz9vzCLwB///cAZ6VqiQFADQOAj/o2DejZs3b/wE63kyfd3+fkSW/vc/as8z2Wl7t/jy7I+e8b1fb/r+0u1CulUgmzs7OBCABA+5l12gUAL/3yuw0AciCy35qdb/Pz89i5cyfm5+d9OCrqJ9ZWX6q1feAB4BWvAPbu7fAFw8OdTcdpGM33vfWt9gdzMHpHOA2oGgYAH/VtGtAtW7zNC/ya17i/j9d5hrdscfsiu3+PLrgtwtXY/79x3n0dtZtZp10AkF2m3O4CZLPZup/LXgQATgNKRIM0M2OPvXXT7HEHy7K7/TRe+e9ku/xy4EMfAiYne/Y1RQWnAVXDAOAj+UN75MiRjp7PAOCdbLjWSqVSdd2CGhu9cj5+Obd/Op1GKpXyNGC4UqlUtyDodCGwVgtwyQXJTNOsrnUg++vXThva6rPkegrN7rh0chyD0ux8q1QqKJVKgakt9Q5rqyfLAh58sIJstoQHH6w4uuD/6Ef2FPof+UiXHzQ25r3h/7KXAZ/7HPDMM11+eHQdOXKEAUABA4CP+jYGgAGgTiqVqi5elUqlUCwWqw18+W+p9m5AMpmszp0vX9+pII0BANa7QjVblKvThnc2m6028GU3qsafyVYBQAaGZgFArjfQOGbDDxwDED2srV7ajcVdWQH+6q+ACy9c3zcy0sUHjox4+3v3v/4XsLbWs683qjgGQA0DgI/6FgA4BkBZbQPZNE3laUKDFgDkYOigd3eS4aHXYzJUMABED2urj07G4l59tbNN/u53d/GhO3d6CwBjY736ciONAUANA4CPOAtQsHUzTWjQAgBgB5pWqwEHgZfZhvqNASB6WFt9dDoWV24XXgik04rjbksl4OtfB37xFzv7sA0b7NsQmi/QNSgMAGoYAHzEABBcuVyubtYbrzMFBTEAyPEQzboB+U12/wnKzzgDQPSwtnrwOhb3JS8BHn5Y4YPW1oCvftWertNL2hDCXmGMeoIBQA0DgI9kANjb4fxiDAD9lcvlqlfIa/v/A86Bwu2USiU8++yzgQoAgH0XIAgz7LgxDKPlWgWD1ux8m5ubw44dOzA3N+fDUVE/sbZ68DoW95vf9PgBlgX80z8Br361+xs+5znNr/wLYa8DwKv/PbN3714GAAUMAD7q2zSgpCSfzyORSCCbqG7LtQAAIABJREFUzSKfzyOZTCKTySCdTgeiT3ovFItFGIahPLahX1KpVMvZgfzA840onLyOxfU88PeTn3R/o1e8AvjSl+yFu7Zutbv51O4fGrIfZ+O/pzgNqBoGAB8xAJAf5IDgoHQFymazMAwjcD/bPN+IwsnrHQDPY3HPnAEuvXT9DV71KvuOQGPD3rLsNx8Zsf/Lhn9fMACoYQDwEccA6CuIYwDIm1ZjALZv385+4hpibfXQ6RiAtmNxl5aAffvc9334w3bf/69+lVN5+oxjANQwAPiIAUBfDADhx0HA0cPahp9cw612FqCLhIX/LMbxWnEI14sxXCSs1mNxFxaAz3/eTgdXXQW4jQmZn7dn/yHfMQCoYQDwEQOAvhgAwo8BIHpY23B78EHgda8DTp2yr+q//a0WbhfDmBb1K4GdEUO4XQzj7W+16q/+z8/bc4G+6EX1twq2bfPta6L2GADUMAD4iAFAXwwA4ccAED2sbTjNzwMf+tB6ez0eByqrFspvtVcCK4v6vkDy3+Ubzs/GMzsLfPrTwAtf6N5X6Hd+x+8vkVpgAFDDAOAjBgB9MQCEH8cARA9rGz47dwIvf7nLzD6/f39nI4Df/GbAMNz3/dqvAd/97nq/IgokBgA1DAA+4ixARMHF840ouObmgA9+0Nlmv/xy4M47SihfPdT5SmCN2xvfCPzbv7HhHxKcBUgNA4CPGACIgovnG1Ew/eAHwMte5my3X389MDUF7/OAyu2//Bd7IAEb/qHCAKCGAcBHXAlYX0FdCZg6x5WAo4e1Dba5OeD973e/6v93fweUy+ef6HUlsOuuA3bt8vVrI3VcCVgNA4CPOAZAXxwDEH4cBBw9rG1w7dkD/NzPOdvub3oT4DhN+74SGAUJxwCoYQDwEQOAvhgAwm9qaooBIGJY2+CamrKv9Ms2+xVXAF/4Qs1V/1rPPgs897ntG/5tVwKjMGAAUBPqAJDL5ZBMJpFOp6v/7ed7yOfKLZVKoVgsKh+/SgA4efKk8ufR4DAAhN/JkycZACKGtQ22O+9c76rvei1sddVOBUNDnV/9d10JjMKEAUBNaANANpuFaZp1j6XTaSQSib68RyaTQSaTqXssl8s5Xu+F1wDw+OOPY3JyEmXXSx4UJAwA4VYulzE5OYnHH3/csY+NRH2xtsEwO2u35RuVy0A263LVv1QC7rnHfT7QZlf+hQA2b+bVfw0wAKgJbQAwDAO5XM718Ww229P3KBaLTYNFOp1WuvMArAeAI0eOdPT8mZkZHDt2jAPUQqBSqVQ3Cp+5uTkcO3YMMzMzjn2VSgWlUom11RBr21+WZXe3Hxmx/+vW9n7gAeClLwW2bOngDWUiuPba5vP4//EfO+8IDA3ZV/7Z+NfCkSNHGAAUhDIAZLNZCOF+6IlEoqO7AF7eQ3b9cZPP5xGPxzs4aievU1etrq5icnISJ06cwNLSEv9IEfVYpVLB4uIiTpw4gcnJSay6XYYkIk8sCxgeBjZudLbDh4ft/c8+C7zvfev7LroIOHiwyRtWKsCOHYBpujf8X/c6YPv29ek8O0keFFqcBlRNKANAMpmEYRiu+1KpVNN9qu+RzWYRi8Vcn5vL5Tx1O6olf2j37dvX8WtmZ2dx7NgxHDt2DCdOnKgOVOQWrG1qagrHjx9nfUK0TU1N4cSJE9Xzq1k3kPn5eezcuRPz8/NK5z0FF2vbe5YF3HBDfc+bxp44b3gDcM01znb8hz/c5E3PnAEuucT5gle/Gvja11xHBrO2+tq3bx8DgIJQBoBYLNa0730mk4EQAvl8vmfvUSgUIISAaZqOQb+JRMK1G1EnvI4BkJaXlzE9PY3HHnvM90YTN/ft+PHj2LVrF44fP+77sXDrfHvssccwPT2N5eXlpucf+4nri7XtveFh7+txPfe5wN13t1mP60/+ZP0FL30p8OUvA2trTZ/O2uqLYwDUhDIAGIbRtNuN7NrTrlHu9T0SiQSEEHXjBjKZDFKplOJXoR4AKPj4x0ZfrK2+WNvesiy720/jlf9WWzwOnDpV8yY/+QmwuOh88+lp4BWvAP72b4EWgV1ibfXFAKAmlAFACNG28d5uILDKe8gQIO8GdDrYuBkGAH3xj42+WFt9sba91bge10XCwhvFw/g58e+ujf+PfrTmqn+hAPzRH9npodlEGx5mWWNt9cUAoEbbANA4ZWcv3qNYLMIwjLoQUCh0tjDX2bNncfTo0bptdHQUQgiMjY1hdna2usmBh6urq3WPz87OYmVlBQBgWZZjn+y24LZvaWkJALC2tubYt3j+6kqpVHLsW1hYAGBPi9i4T/alrFQqjn21MxU17qv9BTw3N+fYJwc3z8/PO/bJKVAXFhYc++SUm4uLi459a+dvDS8tLTn2WecHhC0vLzfdt7Ky4tjXqk7nzp3D6OgoZmZmWKcA10nlfDpz5ky1IcE6BbdOKueTrO1TTz3FOvWgTv/8z4vVhv/tYhjHxM/jxeIJ56ycogwhgPvvL2FuchKr730vKhddtP6EF7wApaef7up8krU9c+YM66RZO+Khhx5iAFDAANDhe+RyOcRiMeRyOaTT6WoIMAyjoxCwZcuW6msatzvvvBOjo6PVbWpqCoC9Emnt46Ojo5icnAQAnD592rFvYmICADA9Pe3YNz4+DsCeSrRx3+7duwGsXyGp3Xbt2gXA/sXUuG/nzp0A7BO+cd+OHTuqX/v27dvr9m3fvr26b8eOHY7Xyl+UO3fudOyTv4B27drl2Cd/IezevduxT07nOD4+7tg3PT0NAJiYmHDsO336NABgcnLSsa9VnY4dO4alpSWcOnWKdQpwnVTPp71796JcLrNOAa+T6vn0wx/+kHXqQZ22bn0YFwkL3xH2KOCy2IA/FPdW2/XPE7P4kngftovN2CjO4MT/8ycouQ3uveQSLGazPJ9Cej71u0533303A4ACbQNAL7sAFQoFR0M/n88jFotV7wS00+oOwJ49eyKb3HmFhXVinVgn1knPOs3MzGLb5Z+qa8w/Ja6GIZ7B74gdOC1eWn18RVzqaPhXLrwQazfdBJw+zTrxfGpap7179zIAKAhlADAMo+0MPp0MAu70PeLxeNPFvuLxeEezDrlRmQaUwmFhYQG7du2q/uIjfbC2+mJte2fnTuB729cwf8VGlEX9KOApEUOl3Wjg3/994NFHe3Y8rK2+OA2omlAGANM0m87Ln0qlIIRo2y3Hy3u0ez8vqw/X4iBgfc3OcsCZrlhbfbG23XvqKeAP/9Buw//sC1cwK57nbQ7Qt70NOHSo58fF2uqLg4DVhDIAyD74blot8KX6HkIIx/z/teLxuNJaAAwA+uIfG32xtvpibdWVy8BddwGGUd+ev1X8TeeN/898pm/Hx9rqiwFATSgDQKFQaHpVPhaLIZlMOh5vbMB7eY92DfxmdxLaYQDQF//Y6Iu11Rdrq+bgQeDXfs3Znn/e5Wv4gvhQ5wFgbKxvx8ja6osBQE0oAwBgz8nfuAhXPp93vVpvmqbr452+hxzw6yaTybSdcagZBgB98Y+NvlhbfbG23szNAf/zfwIXXOBsy7/rXcBPv70HcJvVxzEP6AZgaMheOaxPWFt9MQCoCW0AKBaLME2z2vjO5/MwTdP1Sn0ymUQsFnMEAC/vkc1mq8/N5/PI5XJIJpNNBwd3QgaAw4cPK78HBZOcCUHO/ED6YG31xdp2plIBRkaAl7zE2ZaPxYCH/nHKHsTrpe//1q19PWbWVl+HDx9mAFAQ2gAgyXn5s9lsy376vXqPXnyeJAMAf2iJiCgMZmaAzZud7feLLwa2/c+zWPvAh4HaRbwar/S7/Xvz5r5e/Se9sS2lJvQBIMzkD+3+/fv9PhTqscXFRezevbs6NzLpg7XVF2vbnmUBv/RL9e34G35zDj/98CeB5z7X2ei/7DLgYx8DPvEJu5tP7b6hIfvK/wAa/6ytvvbv388AoIABwEccA6Av9jfVF2urL9a2Mz/6kd1+f/HPrGLvf/sCKldf7Wz4X3gh8P73A08+uf5Cy7IH+o6M2P8d4FV/1lZfHAOghgHARwwA+uIfG32xtvqKUm07aYv/9KfAqVPur//eX+xG6RWvdO/qc+ONwPHj/Tx8z6JU26hhAFDDAOAjBgB98Y+NvlhbfUWhtpYFDA8DGzc6e+MMD9v7y2XgS18CrroK+K3fsgf9Opw+DVx6af2bXH89sHv3wL+mTkShtlHFAKCGAcBHDAD64h8bfbG2+tK9tpYF3HBD6/G4v/mbwK//ev2+f/7nJm/40Y/aT3jta4HvfrdJUggG3WsbZQwAahgAfCQDwKE+LHtO/lpbW8PMzAzW1tb8PhTqMdZWX7rXdnh4vVF/kbDwJjGGd4gRvEmM4UJhufbmiYkpPHTtB4CFBecbzswAX/kKEIKpNXWvbZQdOnSIAUABA4CPOHUVERENgmXZ3X4uFhZuF8OYFut9gEbF2/Fi8URdw/9F4iz+6XkfRvnC81N6fvrTfn8JRK7YllLDAOAj+UM7Pj7u96FQjy0tLWF8fBxLS0t+Hwr1GGurL51rOzZmX/X/jrD7AJXFBpwSP4e3iW/VNfyfK+bwF2ILli5qmNLz+c8Hzp3z+8tQpnNto258fJwBQAEDgI84BkBf7G+qL9ZWXzrXdmQEuF3YfYBWxCXYJv4Ml4uF9YW8xCo+JL6Ap8VV7lN6JpP2tEAhpXNto45jANQwAPiIAUBf/GOjL9ZWXzrX9qGchWmxEWWxAWfFi3ClKNqDf0UZ7xJfR0G83H1Kz3e+E5ic9Pvwu6ZzbaOOAUANA4CPGAD0xT82+mJt9aVzbddyY3UN+78Rt+LN4gfYL37ZteFfee3rAjulpwqdaxt1DABqGAB8xACgL/6x0Rdrqy8da7uwAMzNwe4DVNPAt8SFeFxc437VXwjgX//V70PvKR1rSzYGADUMAD6SAeDgwYN+Hwr1mGVZmJ6ehjXApe5pMFhbfelU20oFuO8+4JprgP/xP2CPAm7W2HfbxsZ8/gp6S6faUr2DBw8yAChgAPARp64iIqJeO3LEXsFXtuWvuuBZPPq9KWDjRlQaVwBr7PqzYYO9LDAbyhQSbEupYQDwkfyhzefzfh8K9djy8jImJiawvLzs96FQj7G2+gp7bYtF4NZb7Ul75Mw+t4q/wTnxQhx7+e/WrwTWatu61e8vpefCXltqLp/PMwAoYADwEccA6Iv9TfXF2uorrLUtl4Evfxl40YtkG76C/1fchykRq2/Y79gBbN5s/3/jnQD5782btbz6H9baUnscA6CGAcBHDAD64h8bfbG2+gpjbffuBd7whvV2/G+KH2KP+FX3K/u33mo37rdutbv51O4bGrIf17DxD4SzttQZBgA1DAA+YgDQF//Y6Iu11VeYanv2LPDe9663339BHMOoeLt7w/91rwMeeMAeGSxZlj3Qd2TE/q+mDX8pTLUlbxgA1DAA+IgBQF/8Y6Mv1lZfYartj350/sK9OIO7RRIlcYGz4X/NNcA99wClkt+H67sw1Za8YQBQwwDgI04Dqi/LsnD69GlOOach1lZfYavt59/0LcyLK5wN/+c/H9i2DVha8vsQAyNstaXOcRpQNQwAPuLUVURE0dZJT5wzZ+p770hP5Z/AyoWXrTf8L74Y+MhHgHPn+n3YRIHBtpQaBgAfyR/aAwcO+H0o1GMrKyuYnJzEysqK34dCPcba6muQtbUse1bOjRudY3GHh+39y8v2uNzLLgO++tUmb/SJT9gvfNe7gKmpvh93WPG81deBAwcYABQwAPiIYwD0xf6m+mJt9TWo2loWcMMNcq5+C28SY3iHGMGbxBguFhaEAH7lV4DY+Vk8f03sxncvfQfmzsy7HTTAvyFt8bzVF8cAqGEA8BEDgL74x0ZfrK2+BlXb4WHgImHhdjGMaVF/C+Bh8Rt4pTgJIYBXipO4XySq+44kPtnX49IZz1t9MQCoYQDwEQOAvvjHRl+srb4GUVvLAl5ytYXvCPsWQFnYC3DNieciJT6Li8Uqfkb8FH8r/gSWuKi+f9AVVwDT0307Np3xvNUXA4AaBgAfMQDoi39s9MXa6msQtR0bA24Xw9VGfUUIfFW8Gy8WT+AysYiPic/gWfF858w+V1wBfOpTwMJC345NZzxv9cUAoIYBwEecBlRfq6urmJqawurqqt+HQj3G2uprELX9xn0WpsVGlMUG5MV/xhvFw7hAlPAe8f/hP8RLHA3/8gUXAh/8IK/8d4nnrb44DagaBgAfceoqIqJoOXjHWLVx/0XxAfyO2IHD4rXuK/gKgWMf+4rfh0wUaGxLqWEA8BHvAOiLV5v0xdrqaxC1XbtvpNq4XxMX4FHx6qaNfwhhP5+6xvNWX7wDoIYBwEccA6Av9jfVF2urr37VtlIBDh06/4+xsboG/iHxupYBAGNjPT2WqOJ5qy+OAVDDAOAjBgB98Y+NvlhbffWjthMTwJvfDDx/wxyOPXAasCxUrt5Ynf2n2VYWG1DZOOS+NDB5xvNWXwwAahgAfMQAoC/+sdEXa6uvXtb23Dl77O4lGyx8QPw9nhJXI3/Vm1EpV+yFAFpd9Zfb1q09+KoI4HmrMwYANQwAPmIA0Bf/2OiLtdVXL2prWcAddwDGlRW8TXwLk+Ln6xr1577yXftJmzfb04BuqL8TUP335s28+t9DPG/1xQCghgHARzIAHDhwwO9DoR5bWVnB5OQkVlZW/D4U6jHWVl/d1vZ73wN+4ReAXxH78JC43v2qfjJpP9my7Cv8Q0P1+4eG7MfZ+O8pnrf6OnDgAAOAAgYAH3HqKiKi8JucBN76VuBl4t/xf8R/dW/4X3st8O1v2yOCa1mWPdB3ZMT+Lxv+RJ6wLaWGAcBHnAZUX5Zl4fTp07D4x1w7rK2+vNa2UgH+9E+Bn7nwGfyVuA0r4hJnw3/jRuDuu4G1tT4fPbXC81ZfnAZUDQOAjzgGQF/sb6ov1lZfKrX9q9/5AWbECxwN/8pllwH/+38Dc3N9PGLqFM9bfXEMgBoGAB8xAOiLf2z0xdrqybKA73xnAX/2Z3vxne8sdNwT5+zhacyLK+oH8b73vcATT/T3gMkTnrf6YgBQwwDgIwYAffGPjb5YW71Ylj0r50uutvAmMYZ3iBG8SYzhmo0WhofXu+QXCsAPfuD+Hkd/7/y0nr/7u/bE/xQ4PG/1xQCghgHARwwA+uIfG32xtvqwLOBtv2vhdjGMabGxrgvPGTGE28UwbvhtC6kUcO3FJ/Evz/ljzD7+rPONFhaAXG7wXwB1jOetvhgA1DAA+EgGgHw+7/ehUI8tLy9jYmICy8vLfh8K9Rhrq4+/3GLhO+KG6sq7tQFgTVyAe8R/Q0ycxB3iI1gVFwNC4MFf/3O/D5sU8LzVVz6fZwBQwADgI05dRUTkD8sCPnuF+4q8j4j/C28Uu3Cb+CsUxZV1+1YvfA7w+ON+Hz4Rnce2lBoGAB9xGlB9WZaF6elpTjmnIdZWDw/lLEyLjXVX/v9DvAR/IO7FfxX/B/8uXuac2efCC4EPfhA4d87vwyePeN7qi9OAqmEA8BHHAOiL/U31xdrq4YefGqs27BfE5dgituAtYgf2i192vStw7ud/w17xi0KJ562+OAZADQOAjxgA9MU/NvpibfVw5FMjgBCYEjFcL8YwKt7u2vCXdwiOfGrE70OmLvC81RcDgBoGAB8xAOiLf2z0xdrqYS03BgiBktiACfFLjob/mriw/t+5Mb8PmbrA81ZfDABqQh0Acrkckskk0ul09b/9fo98Po9EIlHdUqmU6uEzAGiMf2z0xdpqwrIwf4U9BuCAuK5pw78sNmD+uUPoeGUwCiSet/piAFAT2gCQzWZhmmbdY+l0GolEom/vkclkYJomCoVC3XuohgAZAMbHx5VeT8G1tLSE8fFxLC0t+X0o1GOsbTjNzgJ//mdlHPzB+gDe0pb1WYDOiRei4tIFCEKg9MmtPh459QLPW32Nj48zACgIbQAwDAM5l4VXDMNANpvt+Xvk83kYhoFisVj3eCwWc4SITnHqKiKi/iqVgC99CXjHVWMYFyYOPf83USlX7J2WhfINm89f6Xfv+1++YTOv/hMFGNtSakIZALLZLIRwP3TZNafX7xGLxZBMJh3PbfZ4J+QP7aFDh5ReT8G1traGmZkZrK2t+X0o1GOsbXg8+CDwf//8JL4l3lbXuP9x6pvrT7IsYOtWVDYO1U/5uXEI2LqVjX9N8LzV16FDhxgAFIQyACSTSRiG4bovlUo13af6Hv9/e/ceH8dZnwv8tWznnmYT4tiBlMsCAVJCyiZAWiABS21pQjkpEadwCr1n6QU4lNaipaUGqaRsKagf9wI6tFwSWiCKQDUhbtkoK8uO17qs7omSIDlREltaW8lal9VltLvP+UOejVYzq8urWb27v32+n89+aHak0biPRvM+s/POhMNhKKVcPy3YDM4BkIvXm8rFbIvfk08Cv/WeOP5Z/bHjmn4ohambf8n5TZaFmQceQPunP42ZBx7gwF8Y7rdycQ6AnpIsAKtddtPQ0AClFGKxmGfrCAaDeT8t2AwWALl4sJGL2RavRAKo+fgs/qribkyqSx0Df+uiy5D5+y8Bc3Ou389s5WK2crEA6CnJAuDz+VBVVeW6zL60Z62z9RtZh9/vzxaAUCiUfW3mDkAAC4BkPNjIxWyLz+Ii8C//lMYfXXIPRtXPOifxVuxA6k8+AUxMrLoeZisXs5WLBUBPSRYApdSag/e1JgJvZB1KKfh8PjQ0NORMAm5oaIDf73dMDF4vFgC5eLCRi9kWn3/648fQpQKud/CZv/39S9cErQOzlYvZysUCoEdsAWhoaPBsHcsLwEp+v39dnwTE43EMDg7mvJqbm6GUwuHDhzE5OZl9LSwsAAAWFhZy3p+cnMT8/DwAwLIsx7K5cx9ruy2zb322uLjoWJZMJgEAqVTKsWxmZgYAkE6nHcump6cBAJlMxrFsamoq+29fuWz5H+CpqSnHskxm6Q4d09PTjmXpdBoAMDMz41iWSqUAAMlk0rHMnvg1OzvrWGadu9Z3bm4u77L5+XnHstVySiQSiEajOHv2LHMq4px09qd4PI6jR48imUwypwLkNDOzgAcemMG99ybxwAMzmJhYO6fTj53BWfUzOQP/6etuwvxDD21of4rH4zhy5Agmzn1SwJwKvz9NTm7N3z0723g8zpyKOCed/amtrY0FQAMLwDrWoZSCUsr1TP965wfs378/u56VrwMHDqC5uTn7Gh4eBgAMDw/nvN/c3IyhoSEAwOjoqGNZf38/AGBsbMyxzH7WwMTEhGNZNBoF8OIZkuWvtrY2AEt/mFYua2lpAbC0w69cdujQoey//eDBgznLDh48mF126NAhx/fafyhbWlocy+w/QG1tbY5l9h+EaDTqWGYf0Lu6uhzLxsbGAAD9/f2OZaOjowCAoaEhxzLmxJyYk3c5DQ0No7YW2POSWdyqIvh11YRbVQRXXjaFj3/8NCxrKaemJvec/qfy84BSeOGyq9HxZ3+O5h/+kDlxf2JOZZDT1772NRYADWILgNeXAOUb5NfU1Kxr0vFqnwAcO3asbJu71DMs9ja4rZM5FU9OOvtTIpFAPB5HKpViTh7lNDExidt/eR5/rWoxpnbnnMk/pfbgr1Ut3vsrFr78t9P44hV/h2f6n3bkNDMxhZkvfgmT8bj2/pRIJHDq1Kns/7+Zk5wzy3a2iUSCORVxTjr7UzQaZQHQUJIFwOfzrXkHn/VMAl7vOnw+X95bhoZCoXUVDjecAyDX5CSvN5WK2Xrvb/dbeEDdlvMAruwEXrUNB9Xt+Lj6x+wE38iNnyrIdjBbuZitXJwDoKckC0AgEIDf73ddZp+RHxkZ8WwdVVVVa34CoPOMABYAuXiwkYvZesuygC9eXOs6gbdPXY9PqS85JvguqJ0YP7b633gdzFYuZisXC4CekiwA9ll3N6s94Et3HfbXrjYHQOdOQCwAcvFgIxez9VZr2MKY2p1z5j+uduFv1H4cVO91LQZnX/9W4Ny1yl5itnIxW7lYAPSUZAEYGRnJe5bf7/cjGAw63l85QN/IOuyvdTvLHwgE8l5KtBYWALl4sJGL2Xrr8Ocj2YH9nDofB9Sf4Ovq912f4DupLsHgnZ8Fzl2H7TVmKxezlYsFQE9JFgAAqK6udtx+MxaLuZ6NDwQCru9vZB3BYNAxadj+2rUmAOdjF4COjg6t76fiNTMzg7a2tuzkJ5KD2Xpr4PNNyCiF+9X78UW1z/UJvrPqAsyq8wGlMPD5poJtC7OVi9nK1dHRwQKgoWQLQCKRQCAQyN6qMxaLIRAIuJ6lDwaDrg/s2sg6gKW5AMFgELFYLPsQMJ1r/212AeAvLRGVq8VwBGm1DS3qXY6B/6La7ri//2I4YnqTiaiIcCylp2QLgC0cDiMUCqGxsVH7ibwbWYc9+N/Mz7PZv7T9BbiWlcxKp9OYnZ3N3haO5GC2HrMsTF+8Gz3qTTkD/dPqJTn/nVbbMH3JnqVZwwXCbOVitnL19/ezAGgo+QJQyjgHQC5ebyoXs92cqSngxE9TOe+l9i/dBahH3YAz6iVIu0z8hVJIfa6uoNvGbOVitnJxDoAeFgCDWADk4sFGLmarJ5UC7q0/g3+76GPoufSdyKSXTeS1LKRvu/3cmX7lOPMPpZaWF/DsP8BsJWO2crEA6GEBMIgFQC4ebORithv38INz+PKev0dCXZYd2D/yie/lfpFlAXV1yOzek1MAMrv3AHV1BR/8A8xWMmYrFwuAHhYAg1gA5OLBRi5mu35Dj2UQevN3cUK90nE5z9ilrwXcrse2LCASAZqalv53Cwb+NmYrF7OViwVADwuAQbwNqFzT09NoaWnB9PS06U0hjzHbtU1MAF+58yiOq7e5Xst/6i3vQ2bocdOb6cBs5WK2cvE2oHpYAAzirauISJKFBeDfPzOM5p13ug78x68JwPpJxPRmEpEgHEvpYQEwyP6vOGvNAAAgAElEQVSlHRgYML0p5LFMJoNUKoVMgZ5YSuYwW3d9bWfxjcs+iQW10zHwf+GSazD1r/e6X/JTRJitXMxWroGBARYADSwABnEOgFy83lSucsrWSlroqY8guq8JPfURWMn81+O/MPw8Xth2ec7AP7n9Epz+5BeA2dkt3Gp95ZRtuWG2cnEOgB4WAINYAOTiwUaucsjWSlqI7K1FvGJ37iU8FXsQ2VubtwgcvuMr2Sf4PvvePwTGx7d4yzenHLItV8xWLhYAPSwABrEAyMWDjVzSs7WSFjp23ZZzD/7l9+SfURchdMFnMd57yvG9C9MLePLddyHV/6iBLd886dmWM2YrFwuAHhYAg1gA5OLBRi7p2Ub21rpO4E2rbWhSd+C/1HsBpfDjiz9gelM9Jz3bcsZs5WIB0MMCYJBdANrb201vCnlsamoKhw4dwtTUlOlNIY9JztZKWohX7Hac+T+qfgH3qA9jXp2Xfc9SO/D4D2TdwEBytuWO2crV3t7OAqCBBcAg3rqKiIpJT30kZ+D/U+XH19XvYUJd4fhEIKkuRPR9d5veZCIqcxxL6WEBMIi/tERUTKL7mpZu26kuwzfU7+BJ9WrHwD+lKjCsXoWMUktfT0RkEMdSelgADOIcALkmJydx8OBBXm8qkORsO/6hFd9V/xuPqJtd5wGMqFdhftl9/nvqI6Y32VOSsy13zFYuzgHQwwJgEAuAXJxwJpfkbL/3pi+4DvyfVS9FQl2WMyF4vGLPqs8FKEWSsy13zFYuFgA9LAAGsQDIxYONXJKzHfxGe87A/3l1OZ5TV7uWgkhlnenN9ZzkbMsds5WLBUAPC4BBLABy8WAjl5RsMxn394++/INIqgsxrPxIu9wKFEqh/arbxZ39B+RkS07MVi4WAD0sAAaxAMjF603lKvVsZ6Yz+I8P/QgP7fog0otpx/LJoZOYe/wpRCrrMF6xx/kk4Mo6kYN/oPSzpfyYrVwsAHpYAAzizHUi2iqpFNC8vxtt5+3NDuiPfPSeVb/HSlroqY8guq8JPfURsQN/IipdHEvpYQEwiL+0RLQVjnzvOfzX5b/teMDXye3XwDqbNL15RETaOJbSwwJgEJ8ELBefOilXKWX7RGwa37v2s0iqCx2TeF/YuQtPffqrwOKi6c0sGqWULW0Ms5WLTwLWwwJgEOcAyMUJZ3KVQrYT8RT+893/D2Nqt2PgP7ftAgy9/zPInC3e7TelFLIlPcxWLs4B0MMCYBALgFw82MhVzNnOzwP33/XfeLTija637nz0xg9j/slR05tZtIo5W9ocZisXC4AeFgCDWADk4sFGLhPZrncy7tmTMziz7UrHwP+nL70FiXDnlm1vqeJ+KxezlYsFQA8LgEEsAHLxYCPXVmZrJS1E9tYiXrHbeTvOvbWuRaD1N/4l+3XPXXwtnvuX5vw3/acc3G/lYrZysQDoYQEwyC4AAwMDpjeFPJbJZJBKpZDhwEucrcrWSlro2HVbzgO4lj+Q61n1UrRccaejBCzOLeKxK9+BJz52ALB4286N4H4rF7OVa2BggAVAAwuAQbx1FRHlE9lb63od/5S6GPer9+Okuhpt6h2IVNY5v5mDHCIqExxL6WEBMMj+pe3o6DC9KeSx6elptLS0YHp62vSmkMe2IlsraSFesTvnzH9KVeABdRv61Rtz3utQN/IBXR7hfisXs5Wro6ODBUADC4BBnAMgF683lWsrsu2pj+Sc9T+ifhGt6hbXTwS61JuXvp42jfutXMxWLs4B0MMCYBALgFw82Mi1FdlG9zUBSmFIXYsfqduxqLY7Bv4n1dUYUa8ElFr6eto07rdyMVu5WAD0sAAYxAIgFw82cm1FtpHPtuAH6g6cVT/jGPifVT+DAXUdMsve4ycA3uB+KxezlYsFQA8LgEEsAHLxYCNXobONfLYFT297hWPgb6kd6FY3YEHtyLkb0HjFHs4B8Aj3W7mYrVwsAHpYAAyyC0B/f7/pTSGPpdNpzM7OIp1Om94U8lihsx36bo/jtp996nokXD4NgFLudwEiLdxv5WK2cvX397MAaGABMIi3riIiN22v+Z2lOQDbXo9n1UtdnwMApdB+1e08+09EZY1jKT0sAAbxNqByzczMoK2tDTMzM6Y3hTzmVbYnup7Hj99xN1ILKceysa7n8PhnvwNreh6RyjqMV+xxPgm4so6Df49xv5WL2crF24DqYQEwiHMA5OL1pnJtNtsXxubRfOuX8YLyAUqh7be+vub3WEkLPfURRPc1oac+woF/gXC/lYvZysU5AHpYAAxiAZCLBxu5dLO1FjL48e824kSF33E2f3qMDycqBtxv5WK2crEA6GEBMIgFQC4ebOTaaLaZDHD474+j64K3u07i7Xv5e5F49GSBt5rWg/utXMxWLhYAPSwABrEAyMWDjVwbyXbgR0/hoV0fdB34D1/68zjx9Ye2YItpvbjfysVs5WIB0MMCYJBdAPr6+kxvCnkslUphcnISqZRzgieVtvVk+9yjZ/HAz9VgTp3vGPiP73gp+j/1TYC/G0WH+61czFauvr4+FgANLAAG8dZVRKVlPZNxU4sZPHbemxwD/6S6CJ3v+zwWz/IuJEREXuFYSg8LgEH2L21nZ6fpTSGPJZNJRKNRJJNJ05tCHrCSFiJ7axGv2O28HefeWkcROPxb/5Zzz/7OG34f00+eMrT1tF7cb+VitnJ1dnayAGhgATCIcwDk4vWmclhJCx27bst5ANfywb2ltjseyJVaSOGJC65H/9W/hLH/4SV+pYL7rVzMVi7OAdDDAmAQC4BcPNjIEdlb6zqJ90n1ajyk9uKweiegFCKVdTnfN/vc80u3AKKSwf1WLmYrFwuAnpIuAOFwGMFgEKFQKPu/W7mOUCiExsbGDf9MGwuAXDzYyGAlLcQrduec+T+tXoIfq/cgqS7MfgrwmLoW4xV7+ICuEsf9Vi5mKxcLgJ6SLQCNjY0IBAI574VCIVRXV2/JOhKJBJRSaGhoWPfPW4kFQC4ebGToqY9kB/5z6jwcUr+McXWV49OAThUAlFr6eipZ3G/lYrZysQDoKdkC4PP5EA6HXd9f71n5zawjGAx6VgB6e3u110HFaXFxERMTE1hcXDS9KbQJ0X1NyCiFh9WtGFKvc70UqEPdhIS6DFAK0X1NpjeZNoH7rVzMVq7e3l4WAA0lWQAaGxuhlPumV1dXr+sM/mbWEQ6HEQqFPCsA/KUlKk4//M3v45h6m+vAf0Bdh2fUNTnv8RMAIqKtxbGUnpIsAMFgED6fz3VZTU1N3mVeraOmpgYjIyOeFYCuri7tdVBxmp2dRVdXF2ZnZ01vCml4pnMc4VfdhZSqcAz8n1Yvx6C6znE3IM4BKH3cb+VitnJ1dXWxAGgoyQLg9/sd1+7bGhoaoJRCLBYryDoaGhowMjLiaQHgHAB5eL1paXvyh4OOwf/z6nJ0qJuQdvk0wO0uQFR6uN/KxWzl4hwAPSVZAHw+H6qqqlyX2Zf2uF3bv9l1JBKJ7F2CWABoNTzYlL7Drw8CSmFenYef7HwPZtX5rs8BgFKO5wBQaeJ+KxezlYsFQE9JFgCl1JqD97Um8eqso6amJvt/swDQaniwKQ2ZDHDsS0ddB++nB8bR/rqP4Ez7yNKTgCvrMF6xx/kk4Mo6Dv6F4H4rF7OViwVAj9gCsNbAfKPriMViOYVgowUgHo9jcHAw59Xc3AylFCKRCCYnJ7OvhYUFAMDCwkLO+5OTk5ifnwcAWJblWDY3N5d3mX3d4+LiomOZ/Wj0VCrlWDYzMwMASKfTjmXT09MAgEwm41g2NTWV/bevXLb8D/DU1JRjWebcw5Omp6cdy9LpNABgZmbGsSyVSgFYeuT7ymX2nR9mZ2cdyyxrafA2NzeXd9n8/Lxj2Wo5nTlzBs3NzZiYmGBORZpT+z0xHLni15Yu3/nAP60rp4mxCbR+rhn33PFFHLv7QcxOOreTOXm/P23V371Tp06hubkZ4+PjzKmIc9LZn+xsT506xZyKOCed/am1tZUFQAMLwDrXEQwGc75mowVg//79UEq5vg4cOIDm5ubsa3h4GAAwPDyc835zczOGhoYAAKOjo45l/f39AICxsTHHMnui8cTEhGNZNBoF8OIZkuWvtrY2AEt/mFYua2lpAbC0w69cdujQoey//eDBgznLDh48mF126NAhx/fafyhbWlocy+w/QG1tbY5l9h+EaDTqWDYxMQFgabLQymVjY2MAgP7+fsey0dFRAMDQ0JBj2Wo5DQ4OYmxsDCdOnGBORZbTPV/5Pppe9juw1I7smfwz267EfV///rpzOnr0KCzLYk5btD9t9d+9w4cPM6cSyIn7E3MCgK997WssABrEFgAvLwGyJ/4u5+UnAMePHy/b5s4zLMxpq3KaeG4SD9x6d/ae/ctfZ7ddhlj9fzOnIsiJ+xNzYk7MaSM5tbe3swBoKMkC4PP51ryDz3omAa9nHYlEIufaf5uXcwDWumMRlZ65uTn09/dn/5iSOanFDB76g+9itOKVjoG/pXbg2Fs+jumnzqx7fcxWLmYrF7OVKxaLsQBoKMkCEAgE4Pf7XZfV1NRAKeU4Y6+7jsbGRgQCAVRVVeW8AoEAlFLw+/2oqqpyLQlr4SRguSYnOeGsGBz/yiPou9D9QV5d1/wvnGp9YsPrZLZyMVu5mK1cnASspyQLgP0UXjerPeDLy3XYjZN3ASI3PNgUlpW00FMfQXRfE3rqI6534Qlf9wnXgf/jlwTweEOr9s9mtnIxW7mYrVwsAHpKsgDYl9+4neX3+/2OCbvA0j38N7uO5VgAaDU82BSGlbQQ2VuLeMVu5+0499bmFIEjd30752tO7bgGXZ+8F5lUelPbwGzlYrZyMVu5WAD0lGQBAIDq6mrHZTf2oHzlYN++XGfl+xtZx0rhcJgFgPLiwcZ7VtJCx67bch7AtdoDudKLaTx2YQDT6hI8cvsXsHB21pPtYLZyMVu5mK1cLAB6SrYAJBIJBAKB7AA8FoshEAi4Tv4NBoPw+/2OQf1G1mGz5wT4fL7sbTyXr2Mj7ALQ09Oz4e+l4mZZFkZHR7N3aaDNi+ytdVzOk1YKreqdeEi9O/tepLIu+z0nfjSIxOPjnm4Hs5WL2crFbOXq6elhAdBQsgXAFg6HEQqF0NjYuOZZ+0KuQ4ddAPhLS7Q6K2khXrE758x/t7oBMfXm7H8PqjcgrbZhvGIPn8xLRFQmOJbSU/IFoJTZv7Td3d2mN4U8Nj8/j6Ghoez9lmlzeuoj2YH+CfUKtKpbHJ8G9Ko3IXPu/+6pjxRsW5itXMxWLmYrV3d3NwuABhYAgzgHQC5eb+qt6L4mPK98+ImqxKy6wDH4j6tdOK7eki0A0X1NBdsWZisXs5WL2crFOQB6WAAMYgGQiwcb7ywkF3H/mz6PuNrlGPgn1YU4qn4Rc+q8nPcL+QkAs5WL2crFbOViAdDDAmAQC4BcPNhsXiadwdHP/BhP7rzOZfLvNhxTb8OEutzxfqHnADBbuZitXMxWLhYAPSwABrEAyMWDzeY9c+RpLKrtrtf6j6hXuj7ka+VdgAqB2crFbOVitnKxAOhhATCItwGVa2FhAcPDw1hYWDC9KSXt8PV/8uLk3/Nfhwcu/sC6nwNQKMxWLmYrF7OVi7cB1cMCYBBvXUW0ZPLUDBamnQfmiaHTeHrna3DsN/8ZqTlr6UnAlXUYr9jjfBJwZR1v/0lEVGY4ltLDAmAQPwGQi2eb1mdxPoWHP/zvGKu4Gq131Lt+TWYx5XjPSlroqY8guq8JPfWRLR34M1u5mK1czFYufgKghwXAIM4BkIvXm64ukwGifxvG0Pk3ZM/iv7Dtcrww/LzpTVsTs5WL2crFbOXiHAA9LAAGsQDIxYNNfkNNjyL6kttdJ/AeudP9U4BiwmzlYrZyMVu5WAD0sAAYxAIgFw82TmN9cTz8+j90vbPP6M5Xo+sv70cmnTG9mWtitnIxW7mYrVwsAHpYAAxiAZCrnA42a12PP3NmFg9V3o1Jdalj4P/Ctstx9M6vYDFZOtflllO25YbZysVs5WIB0MMCYJBdALq7u01vCnlsfn4eQ0NDmJ+fN70pBWMlLUT21iJesdt5R569tbCSFh75i4N4dvvLHQP/BbUTR276JCafKv5r/lcqh2zLFbOVi9nK1d3dzQKggQXAIN66ikqVlbTQseu2Ne/J3/aH33EM/tuveT+ejfzU9D+BiIgE4FhKDwuAQbwNqFyWZWF0dBSWJfO+9JG9tXmfxLv89fC7azFwydsApfDoJW/Bo19rM73pmyY923LGbOVitnLxNqB6WAAM4hwAuSRfb2olLcQrduec+T+jrsBh9U7HJwHjFXvQ/7VHcPwT/4FMKm160z0hOdtyx2zlYrZycQ6AHhYAg1gA5JJ8sOmpj2QH+XPqPDyk3o0XlA9QCn3qjY5PAXrqI6Y32VOSsy13zFYuZisXC4AeFgCDWADkknywie5rQkYptKm34yn1ipzB/oD6OaRXFIDovibTm+wpydmWO2YrF7OViwVADwuAQSwAckk+2PzoA99Gj7rB9Zr/TnUjJtUl/ASAShKzlYvZysUCoIcFwCC7AMRiMdObQh6bm5tDf38/5ubmTG+KZ54+/DSO/OyHXAf+j6vXoldd7zoHYOVzAUqdxGxpCbOVi9nKFYvFWAA0sAAYxFtXUSl44emzaHnLpzGnzncM/MfVVXhE3ey47Md+RSrrTG8+EREJxrGUHhYAg3gbULksy8LY2FjJ33Ju9oU5PFdxjWNgP6Muwo923oFZdf6qzwGQdvYfkJMtOTFbuZitXLwNqB4WAIM4B0AuSdebRm78VM7g/ui1v4vTvSeXngRcWYfxij3OJwFX1okc/AOysqVczFYuZisX5wDoYQEwiAVArlI92GQyzvcSJ17A89uuQPdLKjF8v/PTKitpoac+gui+JvTUR8QO/G2lmi2tjdnKxWzlYgHQwwJgEAuAXKV2sHmu4yRa/b+LyO1fcl0+1j7q3g7KUKllS+vHbOVitnKxAOhhATCIBUCuUjnYTJ2aRss7/gYz6iJAKZxVP4Mzj502vVlFrVSypY1jtnIxW7lYAPSwABhkF4Curi7Tm0Iem52dRVdXF2ZnZ01viqvF+RQiH/k3jFVc7Zjg23rDx01vXlEr9mxJH7OVi9nK1dXVxQKggQXAIN66ikzo+MJP8Pj5b3K9becj/t/EqeOjpjeRiIhoXTiW0sMCYJD9S9vb22t6U8hji4uLmJiYwOLi4pb9zLUm4z7xg0G0X/mrrgP/vsveicfv7diybS1lJrKlrcFs5WK2cvX29rIAaGABMIhzAOTayutNraSFyN5axCt2O2/HubcWz7U/i9bXfxQpVeEY+D+18zXo+MsfIJPmBN/14rXEcjFbuZitXJwDoIcFwCAWALm26mBjJS107Lot5wFcKx/I9YMLP+QY+D+/7Qq03fmPsGYWCrp9EnEgIRezlYvZysUCoIcFwCAWALm26mAT2VvreknP8ldGKXRtfyugFBbUThy+6VM4+9QLBd0uyTiQkIvZysVs5WIB0MMCYBALgFxbcbCxkhbiFbtzzvw/rl7rKABptQ2PqJtx7GV34pnIcMG2p1xwICEXs5WL2crFAqCHBcAguwB0dnaa3hTyWDKZRDQaRTKZLNjP6KmPZAf5I+qVOKZuBpRCj7rB9ZOAnvpIwbalnGxFtmQGs5WL2crV2dnJAqCBBcAg3rqKNiO6rwkT6nJE1K2w1I7sQP8x9XqkXQpAdF+T6U0mIiLyFMdSelgADLJ/afv6+kxvCnkslUphcnISqVSqIOufOzuP5uv+Agl1mWOgf1ZdimHl5ycABVLobMkcZisXs5Wrr6+PBUADC4BBnAMgV6GuN82kMzj6yfvw9A7nAH9RbUebegcm1OWOOQDjFXsczwUgPbyWWC5mKxezlYtzAPSwABjEAiBXIQ42fQ1R9F3yi67X93eom/CUekXeOwFFKus8245yx4GEXMxWLmYrFwuAHhYAg1gA5PLyYJNJZ9D26t9xHdgPXfjz+K+LP7jqcwDar7qdZ/89xIGEXMxWLmYrFwuAHhYAg1gA5PL6YBO5+S9yBvdjFS/F0bu+ibSVWnoScGUdxiv2OJ8EXFnHwb/HOJCQi9nKxWzlYgHQwwJgkF0AOjo6TG8KeWxmZgZtbW2YmZnxZH2Tz07i9LZdmFYXo7WyFsnTzvVaSQs99RFE9zWhpz7CgX+BeJ0tFQ9mKxezlaujo4MFQAMLgEG8dRUtl0lnEP3MQRz+yNddl/fURxDvPbXFW0VERFS8OJbSwwJgkP1L29/fb3pTyGPpdBqzs7NIp9Pr+vpHv9ON2GXvBpTCtLoY8b6xAm8h6dpotlQ6mK1czFau/v5+FgANLAAGcQ6AXOu93vRk+7No8/+2YwLv4dfftUVbShvFa4nlYrZyMVu5OAdADwuAQSwAcq11sJk6NY2H3/7XSKoLHXf2OV1xFY78tvtlQGQeBxJyMVu5mK1cLAB6SroAhMNhBINBhEKh7P8Wch0NDQ0IBoMIBAIIBAJaP285FgCZrKSFY3c/iHvu+CKO3f1gzmTcxfkUWj/ydcQrdjsG/rPqArS+/TOYPskDVDHjQEIuZisXs5WLBUBPyRaAxsZGBAKBnPdCoRCqq6sLso6amho0NjZm/3tkZAQ+nw9+v3+DW/4iFgBZrKSFyN5ax+B+vGIPIntrEd3/IJ44/3rX+/kf9X8Yp9qfMf1PoHXgQEIuZisXs5WLBUBPyRYAn8+HcDjs+v7ygboX64jFYqipqXF8XTgchlIKwWBwA1v+It4GVA4raaFj1215H8jVp97oOvDvvewWDN3baXrzaQOmp6fR0tKC6elp05tCHmO2cjFbuXgbUD0lWQAaGxuhlPumV1dXr+tTgI2sIxgMIpFIuH6tz+fLu5618NZVckT21roO8O1XRil0qTdn//vEztei/S9/iEw6Y3rTiYiIShbHUnpKsgAEg0H4fD7XZTU1NXmX6a4jEAjA5/O5loBAIAClFEZGRta59S+yf2kHBgY2/L1UPKykhXjF7uyZ/3m107UEPKZeh9PqSkTuqIc1s2B6s0lTJpNBKpVCJsPyJg2zlYvZyjUwMMACoKEkC4Df73dcu29raGiAUgqxWMyzdVRVVeUd5HtRADgHoLT11EfOXeqjcES9HSfV1Tln+5e/5tR5S19PJYvXEsvFbOVitnJxDoCekiwAPp8PVVVVrsvsS3vcru3XXUcikchbKLy4BIgFoLRF9zWhW92AR9UbsgP9J9VrkFoxF8B+Rfc1md5k2gQOJORitnIxW7lYAPSUZAFQSq05eF9rIrAX67AnAbtNEF4PFoDSN3LoCbRd+qvuA331Ntf3+QlAaeNAQi5mKxezlYsFQI/YAtDQ0FDwdQQCAfj9/rwThJeLx+MYHBzMeTU3N0MphUgkgsnJyexrYWHp+vCFhYWc9ycnJzE/Pw8AsCzLsWxubi7vstnZWQDA4uKiY1kymQQApFIpx7KZmRkAS49RX7nMvptCJpNxLJuamsr+21cuW/4HeGpqyrHMvkZzenrascx+jPvMzIxjWSqVAgAkk0nHssXFRQDA7OysY5llLd2nf25uLu+y+fn5nPef7noKD1//MVhqh2OA/4y6BsfUzciseD+ttmG8Yg+mnp9iTluUUyH2p1OnTmUHEsypeHPS2Z/sbMfHx5lTEeeksz/Z2Z46dYo5FXFOOvtTa2srC4AGFgDNdYRCoQ1d+79//34opVxfBw4cQHNzc/Y1PDwMABgeHs55v7m5GUNDQwCA0dFRx7L+/n4AwNjYmGNZV1cXAGBiYsKxLBqNAnjxDMnyV1tbG4ClP0wrl7W0tABY2uFXLjt06FD2337w4MGcZQcPHswuO3TokON77T+ULS0tjmX2H6C2tjbHMvsPQjQadSybmJgAAHR1dTmWjY2NAQD6+/sdy0ZHRwEAQ0NDaG5uxv3/0YTv3vhnSKjLHAP/s+pStKpb8k4ChlK4702fYE5bkFOh96eDBw9iamqKORV5Trr7UyQSYU4lkBP3J+YEAF/96ldZADSILQCFvATIfgjYWvMMllvtE4Djx4+XbXMvpTMsrf/3uxjd/irHoH5RbcfDb/goWi6/I+9zAKAUju+6DRNjE8yJZ8KYE3NiTsyJOXmUU3t7OwuAhpIsAD6fb807+KxnErDuOgKBwIYG//nw3rWlJXLL3zgG/+17fg0nHlw6m2IlLUQq6zBescf5JODKOlhJy/C/gIiISBaOpfSUZAGwr713U1NTs65Lc3TXUV1dve4nDa+Fk4BLy/TYdHZwP3Thm9HzlYddv85KWjh294O499dDOHb3gxz4CzM5OYmDBw/mnIEiGZitXMxWLk4C1lOSBcC+/t7Nag/42uw6QqGQ67yAhoYGPgdAkMRTCXR+4X9clx370/tw9K5vIb2YXnUdk5O844RUzFYuZisXs5WLBUBPSRaAkZGRvGfo/X4/gsGg4/2Vd+rZ6DoaGxsRCoVct6e6unojm5/FAlBcrKSFyJ0HMLHtJUiqC3Gq41ntdfFgIxezlYvZysVs5WIB0FOSBQBYGnSvvP9+LBaDUsox2Lef1rvy/fWuw570W1VVlfOyLyPKdynRWlgAikMmnUH0L5oxsvPanGv3j7z6t7XXyYONXMxWLmYrF7OViwVAT8kWgEQigUAgkL0kJxaL5Z2cGwwGXe/Xv9512AUi3yvf3YTWwgJQWFbSQk99BNF9Teipj7hei//od2Lovuxdrrfs7Nz1HizO6l2/z+tN5WK2cjFbuZitXCwAekq2ANjC4TBCoRAaGxvX9UCuQq1DB2euF4aVtBDZW4t4xW7n3Xj21sJKWjjZ/iyO+H/LcctOKIUnz/85dP3tobV/EBERERnFsZSeki8ApYy/tN6zkhY6dt2W9378U+piHNzx60iqCx0D/9MVV6Htww1IzS+a/mcQERHROjtFq9oAABSqSURBVHAspYcFwCD7l7a9vd30pogR2Vub9ym8cXUl4mqX4/1ZdQEib/8rTJ2cWvsHrNPU1BQOHTqU8yATkoHZysVs5WK2cvFBYHpYAAziHABvWUkL8Yrdrpf1ZK/rVzfm/PdR/0dwqv0Zz7eFE87kYrZyMVu5mK1cnAOghwXAIBYAb/XUR/IO/O3XT5UfKVWBHnUDDv3mPQXbFh5s5GK2cjFbuZitXCwAelgADGIB8FZ0X1P2Up/D6p1oVze5loAn1GuQUWrp6wuEBxu5mK1czFYuZisXC4AeFgCDWAC8dbzuf/CwuhVT6hJAKQwrPyy1Pe+nAT31kYJtCw82cjFbuZitXMxWLhYAPSwABtkFYGBgwPSmlLT0Yhptd30bJyuucQzyW9UtjvfSahvGK/a4PhfAK5lMBqlUCplMpmA/g8xgtnIxW7mYrVwDAwMsABpYAAziras2r/srETx2YcD1DP8xdTOeUS9zXRaprDO96URERLRJHEvpYQEwyP6l7ejoML0pJWfkwcdxfPf7XAf3veoGDKjrXJ8DAKXQftXtBT37DwDT09NoaWnB9PR0QX8ObT1mKxezlYvZytXR0cECoIEFwCDOAdi4M0Nn0Hr9x2CpHY6B/+iOV+HYJ7+PhekFRCrrMF6xx/kk4Mq6gg/+AV5vKhmzlYvZysVs5eIcAD0sAAaxAGxc5Ff+zjHwP6suQ+t7v4T5yfmcr7WSFnrqI4jua0JPfWRLBv42HmzkYrZyMVu5mK1cLAB6WAAMYgHYuOSZJE5uX5rsa6kdaL3h43j+iTOmN8uBBxu5mK1czFYuZisXC4AeFgCDWABWd7LjOdf3j3z0HrTveR9OHHp8i7do/XiwkYvZysVs5WK2crEA6GEBMMguAP39/aY3paiMRkZw7JoPYFZdgGePPOVYnkkX/23c0uk0ZmdnkU6nTW8KeYzZysVs5WK2cvX397MAaGABMIi3rsqVeCqByI1/hnl1Xvb6/qOv+JDpzSIiIqIixbGUHhYAg8rpNqCrTchdmLEQef8BTGx7iWOC76S6FKcH4wa3XM/MzAza2towMzNjelPIY8xWLmYrF7OVi7cB1cMCYFA5zAGwkhYie2sRr9jtuCXnw+/6PI7+eRNGdl7rGPinVAXaXn8X4n1jpv8JWni9qVzMVi5mKxezlYtzAPSwABgkvQBYSQsdu27LeQiX/RpQb0C3usH1QV6dV/4KnvzBgOnN3xQebORitnIxW7mYrVwsAHpYAAySXgAie2sdg3tLbccR9XbXgf+T578RXV/4b9Ob7QkebORitnIxW7mYrVwsAHpYAAySXACspIV4xW7HmX8ohePqrTn/HVe70Pp/vorU/KLpzfYMDzZyMVu5mK1czFYuFgA9LAAG2QWgr6/P9KZ4rqc+4nqWH0rhhHoFLLUDSXUhIuoWTKuLlr5ekFQqhcnJSaRSKdObQh5jtnIxW7mYrVx9fX0sABpYAAySfOuqY3/ehA51I7rUm11LwCPqF3BKvTgxOLqvyfQmExERUYmRPJYqJBYAg+xf2s7OTtOb4qkn7u9H9KJ3A0rhafVyLKgdeT8NsF/SPgFIJpOIRqNIJpOmN4U8xmzlYrZyMVu5Ojs7WQA0sAAYJG0OwHjvGA5f+wdIqYqcwX1E3ZJ34J9W2zBesSfnuQAS8HpTuZitXMxWLmYrF+cA6GEBMEhKAZg5ncTDe+swrS52DPAn1BVoy3PXn2xBqKwz/U/wHA82cjFbuZitXMxWLhYAPSwABpV6AUgvptH2B9/GqYqXOQb18+o8PHzTn6P1ijtcnwNg/3f7VbeLO/sP8GAjGbOVi9nKxWzlYgHQwwJgUCkXgNiXI3j0woD7BN+f/Q08c/gEgHNPAq6sw3jFHseTgCOVdSIH/wAPNpIxW7mYrVzMVi4WAD0sAAbZBaC3t9f0pmzIsT+733Xg33/JL2Dg61HX77GSFnrqI4jua0JPfUTswN+2uLiIiYkJLC7KebYBLWG2cjFbuZitXL29vSwAGlgADCrVW1fNnZ3H6I5XZQf+oztehWN/eh8y6YzpTSMiIqIyUqpjKdNYAAyyf2m7urpMb0pe6cW06/uPfPL7SGzzIfLef8D85PwWb1Xxm52dRVdXF2ZnZ01vCnmM2crFbOVitnJ1dXWxAGhgATComOcAZNIZHP3YdzGy81o8/dBPXZcnTrxgYMtKA683lYvZysVs5WK2cnEOgB4WAIOKtQD0/etRDFz8thef0vvS95vepJLDg41czFYuZisXs5WLBUAPC4BBJgrAapNxn24ZxrGXVbtO8H30Wx1bto0S8GAjF7OVi9nKxWzlYgHQwwJg0FYWACtpIbK3FvGK3Y7bcR66+XN4+Oc/iQW10zHwP7n9Ghz9w3vzzgUgdzzYyMVs5WK2cjFbuVgA9LAAGGQXgJ6enoL+HCtpoWPXbY4Hci2oHYioW/C8utwx8J9SlyDyS1/A7POcMKXDsiyMjY3BsmTf7rQcMVu5mK1czFaunp4eFgANLAAGbdWtqyJ7ax0D/OPqrXhKvcLxfkpV4PAbPorTA+MF3SYiIiKizeJtQPWwABhk/9LGYrGC/QwraSFesTvnzD+UwjF1s2spGPpeYT+NKBdzc3Po7+/H3Nyc6U0hjzFbuZitXMxWrlgsxgKggQXAoK2YA9BTH3Gd1DuqrsG8Og9QCk+o16JLvRlQaunradN4valczFYuZisXs5WLcwD0sAAYtBUFIPLR7+Jxda1rCYioW9Gm3oHUsk8HovuaCrYt5YQHG7mYrVzMVi5mKxcLgB4WAIMKWQAW5xZx+ENfRXzbVXhWvQxz5872r/XiJwDe4MFGLmYrF7OVi9nKxQKghwXAoEIUgEw6g47P/Rg/Pf+6nIH9w+rWVQf+abUN4xV7cp4LQPp4sJGL2crFbOVitnKxAOhhATDI69uAPv79XnRdUeU6wI+qtyKzxtn/SGWdJ9tBS7ecGx0d5S3nBGK2cjFbuZitXLwNqB4WAIO8unXVWOwk2l77e447/UApjJz3Ohz79A9x/ErncwCW/3f7Vbfz7D8RERGVFN4GVA8LgEH2L213d7fW98/EZxB51+cwoy5yDPzPbLsSrR/45+yg3kpaiFTWYbxij+NJwJHKOg7+PTY/P4+hoSHMz8+b3hTyGLOVi9nKxWzl6u7uZgHQwAJg0GbmAPz0vx7FqYqXOgb+8+o8RN5ag7OjZ12/z0pa6KmPILqvCT31EQ78C4TXm8rFbOVitnIxW7k4B0APC4BBmykAC9MLeGrna3IG/4+8/IN49shT3m8obRgPNnIxW7mYrVzMVi4WAD0lXQDC4TCCwSBCoVD2fwu5Di9+3nKbvQtQtOYHgFLou/QXMfBvhXuWAG0cDzZyMVu5mK1czFYuFgA9JVsAGhsbEQgEct4LhUKorq4uyDq8+HkrracAnHnsNFqv/xOM/HjIsSyTziAWCiOTzmhvAxUGDzZyMVu5mK1czFYuFgA9JVsAfD4fwuGw6/uNjY2er8OLn7fSarcBnUvMIfKeL2JSXQooheO7f03rZ5AZCwsLGB4exsLCgulNIY8xW7mYrVzMVi7eBlRPSRaAxsZGKOW+6dXV1es6K7+RdXjx89zYBeD+mm9kJ+Nm0hk88rH/xLPbX+GY4Nv95Ye1fg4RERGRRLwNqJ6SLADBYBA+n891WU1NTd5luuvw4ue5yf7SnrsdZ+O1f4mBi9/q/iCvq38dT7cMa/0c2no82yQXs5WL2crFbOXiJwB6SrIA+P1+x/X4toaGBiilEIvFPFuHFz/PjV0AfqL24Jh6m+vA/9GLbkLvgcMbXjeZxetN5WK2cjFbuZitXJwDoKckC4DP50NVVZXrMvtyHbfr9XXX4cXPc5OdA6B2OAb+z23/WRz9o+8gvZje8HrJPB5s5GK2cjFbuZitXCwAekqyACil1hyQrzUxdyPr8OLnuVl+CZA98J9Ul+K/1S9j8ln+kSplPNjIxWzlYrZyMVu5WAD0iC0ADQ0Nnq3Di58Xj8cxODiY87rvvvuglEKzUuhT2/AtFUCb8mHw3MTgnp4eHD9+POfV3d2NwcFB12WxWCzvsq6uLgwODqK3t9exrLOzE4ODg+jr63Ms6+jowODgIPr7+/MuGxgYcCxrb2/P/jtXLjt+/Hh2WXt7u2PZwMAABgcH0dHR4VjW39+fd1lfXx8GBwfR2dnpWNbb24vBwUF0dXU5lvX09GBwcBCxWCzvsu7u7rzL3P7/feTIERw4cABHjx5lTkWck87+FIlEcODAgez2MqfizElnf7KzbW1tZU5FnJPO/mRnG4lEmFMR56SzP33zm9/UvhS7nLEArGMdXvy8/fv3QynFF1988cUXX3zxxZfHr29961urjsMol9gCUGyXALl9AvCd73wHSincd999jmV8lfarubkZSik0Nzcb3xa+mC1fzLbcX8xW7su+muLYsWOrjsMoV0kWAJ/Pt+ZdedYzCXi96/Di57kZHOS9a6VitnIxW7mYrVzMVi5mq6ckC0AgEIDf73ddVlNTA6UURkZGPFuHFz/PDX9p5WK2cjFbuZitXMxWLmarpyQLQCgUglLum77aQ7t01+HFz3PDX1q5mK1czFYuZisXs5WL2eopyQIwMjKS96y73+9HMBh0vJ9IJLTXofPz1oO/tHIxW7mYrVzMVi5mKxez1VOSBQAAqqurUVNTk/NeLBaDUsox2A8EAq7vb2QdG/na9YrH49i/fz/i8bjW91PxYrZyMVu5mK1czFYuZqunZAtAIpFAIBDI3n4zFoshEAi4TsYNBoPw+/2OgfpG1rGRryUiIiIiKlYlWwBs4XAYoVAIjY2N2mfiN7IOL34eEREREZEpJV8AiIiIiIho/VgAiIiIiIjKCAuAx8LhMILBIEKhUPZ/TayDvOdlLrFYDNXV1R5uHW2WF/k2NDQgGAwiEAggEAhw3y0SXmYbDAZRVVXluCkEmVGI46V9mS+Ztdls7e9ZfgdH+9gbi8W83tySwwLgocbGRscTg0Oh0IYGel6sg7znZS7206Orqqq82jzaJC/yrampyRk0jIyMwOfz5X2IIG0NL7INBoOOAeFqD4ikrVGI42UikYBSKnvDDzLDi2yrq6uhlMq+fD4flFI8MXMOC4CHfD6f612BfD7fus8meLEO8p4XuVRXV6OqqgoNDQ3w+XwsAEVks/nGYjHXM8LhcBhKKe1nhdDmbTbbkZER1321sbGRgwnDCnG8DAaDLABFwKtjrn0bePuYyzP/L2IB8Ih9MHBTXV29rtbqxTrIe4XIhQWgeHiRbzAYzHtXMPusE209L/8u53sODPdjMwrxd9m+yx8LgFleZcsTL6vjUckjwWAQPp/PdVlNTU3eZV6vg7xXiFxYAIqHF/kGAgH4fD7XEmCfgXJ7kjgVlhfZ2pdyrRwQ2gWAJ2bMKMTf5ZqaGoyMjLAAGOZVtiwAq2MB8Ijf73dcr2azr/le66MnL9ZB3itELiwAxcOLfKuqqvIO8lkAzCnk31T7+zlQNMPrbBsaGjAyMsICUAS8ypYFYHUsAB5ZbUBnf5y11lODvVgHea8QubAAFA8v8k0kEnkPSLwEyJxC/k3lJGCzvMw2kUhk53KwAJjnVbb2pZk1NTU5L1rCo5JHVrsW1P6FXWviihfrIO8VIhcWgOJRyP3OngTMg44Zhcq2oaEBfr+fn+oY5GW2y/dPFgDzvMq2urra8be3pqaGx95zWAA8sp5f2LX+oHixDvJeIXJhASgehdzv7LPE+SYIU2F5ma19JrGqqgp+v5+XYxrmVbaxWMxx+14ea83yMtt862e+LACeYQGQiwVAtkLtd/bdRHiW2JxCZZtIJOD3+zkB2CCvsl15nTgLgHmFHgvx8r0lLAAe4SVAcvESINkKka995xjO2TGrkH9T7bsAcaKhGV5ka0/8XY4FwLxCj4XsmzaUO/5/wCM+n2/NWevrmQS82XWQ9wqRCwtA8ShEvoFAgPtqESj031R7gjc/5dl6m83WvqRrJRYA8wq939pPCC73y/hYADyy2kdKNTU16zpIeLEO8l4hcmEBKB5e51tdXc1P6oqEF9k2NDTkncRt3+KVeW+9zWbb2NiIQCCAqqqqnJedqd/vR1VVFSfwG+DFfmsP8t3YnwCU+9wsFgCP2Nf7ulntoRZer4O8V4hcWACKh5f5hkIh1zOHbpcaUOF5ka1SKu+Aw+/385NZQwp1vLQv7eInAOZ4ke1qJcJ+cGO5YwHwiP2xYb6DhNt1oivbp846qPC8yHYlFoDi4VW+jY2N2XuJr8TJomZ4kW0gEMh7FtguB+V+JtGEQvxdBlgAioEX2dbU1LjmnUgkeGvmc1gAPOR2z1n7j4nbQcXt/Y2sg7aOF9kup5TKe40jbb3N5mtP+nW7nMDv9/OOEwZtNtt8xc6+Fjlf6aPC8/rvMvDisztYAMzabLaJRMK1KNTU1PDs/zksAB5KJBIIBALZPxyxWCzvZMBgMOh6f/CNrIO2jhfZ2vcQtycOLr/OlNcQm7XZfO0DUL4XP+0xx4t9t6GhAdXV1WhoaEAsFkNDQwN8Ph8H/4Z5ka3NnhOw/O/z8nXT1vIi23A4jOrqaoTDYYTDYQSDQQQCAZ5MPYcFoADC4TBCoRAaGxu1f9G8WAd5j7nIxnzl2my2iUQi+2kAfz+KC/dbubzab+3yTi9iASAiIiIiKiMsAEREREREZYQFgIiIiIiojLAAEBERERGVERYAIiIiIqIywgJARERERFRGWACIiIiIiMoICwARERERURlhASAiIiIiKiMsAEREREREZYQFgIiIiIiojLAAEBERERGVERYAIiIiIqIywgJARERERFRGWACIiIiIiMoICwARERERURlhASAiIiIiKiMsAEREREREZYQFgIiIiIiojLAAEBERERGVERYAIiIiIqIywgJARERERFRGWACIiIiIiMoICwARERERURlhASAiIiIiKiMsAEREREREZYQFgIiIiIiojLAAEBERERGVERYAIiIiIqIywgJARERERFRGWACIiIiIiMoICwARERERURlhASAiIiIiKiMsAEREREREZeT/A2lwMzr3iVpDAAAAAElFTkSuQmCC" width="640">
