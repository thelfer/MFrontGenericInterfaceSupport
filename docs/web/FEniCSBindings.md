% `FEniCS` and `MFront` for complex non linear solid mechanics simulation
% Thomas Helfer, Jérémy Bleyer
% 15/01/2019

This document is meant to show how `FEniCS` and `MFront` can be combined
to build almost arbitrary complex non linear solid mechanics simulation.
As `MFront` currently allows the implementation of small and finite
strain mechanical behaviours, the current scope of this coupling is
limited to standard first order theories[^1]. Generalised behaviours,
usable in monolithic multiphysics simulations and/or higher order
mechanical theories, will be available in the next version[^2].

`FEniCS` will handle:

- The computation of the gradients (the strain in the infinitesimal
  strain assumption, the deformation gradient in finite strain).
- The assembly of the internal and external forces.
- The assembly of the stiffness matrix

The `MFront` behaviour will handle, at each integration point:

- The computation of the thermodynamic forces (the Cauchy stress here).
- The update of the internal state variables.
- The computation of the consistent tangent operator.

The interface with `FEniCS` will be handled by the
`MFrontGenericInterfaceSupport` project[^3]. This project aims at
proving tools (functions, classes, bindings, etc...) to handle
behaviours written using `MFront` generic interface. Written in `C++`,
it provides bindings for various languages: `C`, `python`, `fortran`.
This library also provides the so-called `FEniCS` bindings described
below.

This document considers two approaches:

- The first one will rely on the `python` bindings. This approach
  requires the user to write its own version of the Newton solver.
- The second one, strongly inspired by the `fenics-solid-mechanics`
  project written by Kristian B. Ølgaard and Garth N. Wells[^4], will
  allow the usage of the standard Newton-Raphson procedure delivered by
  `dolfin`. This approach is the basis of the `FEniCS` bindings.

Each approach has its strengths and weaknesses:

- The first approach is very flexible can be easily used to describe a
  large number of situations. It is however probably be less efficient
  numerically. In particular, the consistent tangent operator is stored
  at each integration point, which can be a memory bottleneck for large
  simulations.
- The second approach is currently limited to \(3D\) simulations written
  in `C++` and small strain behaviours[^5]. The computation of the
  gradients and the construction of inner forces and the consistent
  tangent operator are handled at element level: this explain why only a
  limited range of behaviours are currently supported. This approach
  also requires to go deep in the very gory details of `FEniCS`
  internals which are far from trivial and, to the authors opinion,
  currently insufficiently described for new comers in the source code
  and the `doxygen` documentation.

# Coupling at the `python` level

## Rapid description of useful classes and functions introduced by the `MFrontGenericInterface` project

In this paper, a limited number of classes and functions of the
`MFrontGenericInterface` project will be considered:

- the `Behaviour` class handles all the information about a specific
  `MFront` behaviour. It is created by the `load` function which takes
  the path to a library, the name of a behaviour and a modelling
  hypothesis.
- the `MaterialDataManager` class handles a bunch of integration points.
  It is instantiated using an instance of the `Behaviour` class and the
  number of integration points[^6]. The `MaterialDataManager` contains
  various interesting members:
    - `s0`: data structure of the `MaterialStateManager` type which
      stands for the material state at the beginning of the time step.
    - `s1`: data structure of the `MaterialStateManager` type which
      stands for the material state at the end of the time step.
    - `K`: a `numpy` array containing the consistent tangent operator at
      each integration points.
- the `MaterialStateManager` class describe the state of a material. The
  following members will be useful in the following:
    - `gradients`: a numpy array containing the value of the gradients
      at each integration points. The number of components of the
      gradients at each integration points is given by the
      `gradients_stride` member.
    - `thermodynamic_forces`:a numpy array containing the value of the
      thermodynamic forces at each integration points. The number of
      components of the thermodynamic forces at each integration points
      is given by the `thermodynamic_forces_stride` member.
    - `internal_state_variables`: a numpy array containing the value of the
      internal state variables at each integration points. The number of
      internal state variables at each integration points is given by the
      `internal_state_variables_stride` member.
- the `setMaterialProperty` and `setExternalStateVariable` functions can
  be used to set the value a material property or a state variable
  respectively.
- the `update` function updates an instance of the
  `MaterialStateManager` by copying the state `s1` at the end of the
  time step in the state `s0` at the beginning of the time step.
- the `revert` function reverts an instance of the
  `MaterialStateManager` by copying the state `s0` at the beginning of
  the time step in the state `s1` at the end of the time step.
- the `integrate` function triggers the behaviour integration at each
  integration points. Various overloads of this function exist. We will
  use a version taking as argument a `MaterialStateManager`, the time
  increment and a range of integration points.

Those classes and functions are available in the `mgis.behaviour`
`python` module.

### Initialisation of a single material

Assuming that number of integration points is known and stored in the
`ng` variable, a possible initialisation of a behaviour and a material
data manager can be:

~~~~{.python}
# modelling hypothesis
h = mgis_bv.Hypothesis.PlaneStrain
# loading the behaviour        
b = mgis_bv.load('src/libBehaviour.so','IsotropicLinearHardeningPlasticity',h)
# material data manager
m = mgis_bv.MaterialDataManager(b,ngauss)
for s in [m.s0, m.s1]:
    mgis_bv.setMaterialProperty(s, "YoungModulus", 150e9)
    mgis_bv.setMaterialProperty(s, "PoissonRatio", 0.3)
    mgis_bv.setMaterialProperty(s, "HardeningSlope", 150e6)
    mgis_bv.setMaterialProperty(s, "YieldStrength", 200e6)
    mgis_bv.setExternalStateVariable(s, "Temperature", 293.15)
~~~~

### An example of a Newton-Raphson loop for a single material

The following code shows how a simple Newton-Raphson can be setup in the
case of a single material:

~~~~{.python}
for (i, t) in enumerate(load_steps):
    loading.t = t
    A, Res = assemble_system(a_Newton, res, bc)
    nRes0 = Res.norm("l2")
    nRes = nRes0
    u1.assign(u)
    print("Increment:", str(i+1))
    niter = 0
    while nRes/nRes0 > tol and niter < Nitermax:
        solve(A, du.vector(), Res, "mumps")
        # current estimate of the displacement at the end of the time step 
        u1.assign(u1+du)
        # project the strain on the Quadrature space with `MFront` conventions
        local_project(eps_MFront(u1), Weps, deps)
        # assigning the values in the material data manager
        m.s1.gradients[:, :] = deps.vector().get_local().reshape((m.n, m.gradients_stride))
        # behaviour integration
        it = mgis_bv.IntegrationType.IntegrationWithConsitentTangentOperator
        mgis_bv.integrate(m, it, 0, 0, m.n);
        # getting the thermodynamic forces and the consistent tangent operator for FEniCS
        sig.vector().set_local(m.s1.thermodynamic_forces.flatten())
        sig.vector().apply("insert")
        Ct.vector().set_local(m.K.flatten())
        Ct.vector().apply("insert")
        # solve Newton system
        A, Res = assemble_system(a_Newton, res, bc)
        nRes = Res.norm("l2")
        print("    Residual:", nRes)
        niter += 1
    # update for new increment
    u.assign(u1)
    mgis_bv.update(m)
~~~~

`m` is an instance of the `MaterialDataManger` class. The key function
here is `local_project` which computes the current estimate of the
strain at the end of the time step (stored using `MFront` conventions).

This code also relies on the previous definition of the consistent
tangent operator and the stress tensor, which is standard `FEniCS` code:

~~~~
a_Newton = inner(eps_MFront(v), dot(Ct, eps_MFront(u_)))*dxm
res = -inner(eps(u_), as_3D_tensor(sig))*dxm + F_ext(
~~~~

This code suffers from the following drawbacks:

- Integration failures are not handled
- Useless copies from `FEniCS` to `MGIS` and back:
- Parallelization:

# Coupling at the `C++` level

## `UFL` file

~~~~{.python}
from math import sqrt

scheme= "default"
degree = 3
dx = Measure("dx")
dx = dx(degree=degree, scheme=scheme)

elementA = VectorElement("Lagrange", tetrahedron, 1)
elementS = VectorElement("Quadrature", tetrahedron, degree, dim=6, quad_scheme=scheme)
elementT = VectorElement("Quadrature", tetrahedron, degree, dim=36, quad_scheme=scheme)

v = TestFunction(elementA)
u = TrialFunction(elementA)

t = Coefficient(elementT)
s = Coefficient(elementS)

f = Coefficient(elementA)
h = Coefficient(elementA)
disp = Coefficient(elementA)

# eps_xx, eps_yy, eps_zz, sqrt(2)*eps_xy, sqrt(2)*eps_xz, sqrt(2)*eps_yz
def eps(u):
    return as_vector([u[i].dx(i) for i in range(3)] +
                     [(u[i].dx(j) + u[j].dx(i))/sqrt(2) for i, j in [(0, 1), (0, 2), (1, 2)]])

def sigma(s):
    return as_matrix([[s[0], s[3]/sqrt(2), s[4]/sqrt(2)],
                      [s[3]/sqrt(2), s[1], s[5]/sqrt(2)],
                      [s[4]/sqrt(2), s[5]/sqrt(2), s[2]]])

def tangent(t):
  return as_matrix([[t[i*6 + j] for j in range(6)] for i in range(6)])

# Bilinear and linear forms
a = inner(eps(v), dot(tangent(t), eps(u)) )*dx
L = inner(grad(v), sigma(s))*dx - dot(v, f)*dx - dot(v, h)*ds(1)

# Functionals to compute displacement and load
M_d = disp[0]*ds(1)
M_f = h[0]*ds(1)

forms = [a, L, M_d, M_f]
~~~~

~~~~{.cpp}
// Extract elements for stress and tangent
std::shared_ptr<const dolfin::FiniteElement> element_t;
std::shared_ptr<const dolfin::FiniteElement> element_s;
{
  MGISSmallStrainFormulation3D::BilinearForm::CoefficientSpace_t Vt(mesh);
  element_t = Vt.element();
}

MGISSmallStrainFormulation3D::LinearForm::CoefficientSpace_s Vs(mesh);
element_s = Vs.element();
~~~~

~~~~{.cpp}
auto b = mgis::behaviour::load(library, "Elasticity",
                               mgis::behaviour::Hypothesis::TRIDIMENSIONAL);
auto m = mgis::fenics::NonLinearMaterial(u,element_t,element_s,b);
const auto yg = 150e9;
const auto nu = 0.3;
setMaterialProperty(m.s0, "YoungModulus", yg);
setMaterialProperty(m.s1, "YoungModulus", yg);
setMaterialProperty(m.s0, "PoissonRatio", nu);
setMaterialProperty(m.s1, "PoissonRatio", nu);
setExternalStateVariable(m.s0,"Temperature", 293.15);
setExternalStateVariable(m.s1,"Temperature", 293.15);
}
~~~~

~~~~{.cpp}
// Create forms and attach functions
auto a = std::make_shared<MGISSmallStrainFormulation3D::BilinearForm>(V, V);
a->t =  m.getTangentOperatorFunction();
a->ds = load_marker;
auto L = std::make_shared<MGISSmallStrainFormulation3D::LinearForm>(V);
L->f = f;
L->h = boundary_load;
L->s = m.getThermodynamicForcesFunction();
~~~~

~~~~{.cpp}
// Solver loop
mgis::size_type step = 0;
mgis::size_type steps = 10;
double dt = 0.1;
while (step < steps) {
  auto extract = [](const std::vector<double> v) {
    return std::array<double, 6>{v[0], v[1], v[2], v[3], v[4], v[5]};
  };
  m.setTimeIncrement(dt);
  t += dt;
  ++step;
  std::cout << "step begin: " << step << std::endl;
  std::cout << "time: " << t << std::endl;
  // solve the non-linear problem
  nonlinear_solver.solve(nonlinear_problem, *u->vector());
  // update state variables
  mgis::behaviour::update(m);
  // post-processings
  const double u_avg = assemble(*M_d);
  ux.push_back(u_avg);
  const double force = assemble(*M_f);
  fx.push_back(force);
  e.push_back(extract(m.s1.gradients));
  s.push_back(extract(m.s1.thermodynamic_forces));
}
~~~~


# References

[^1]: Cohesive zone models are not considered here, mostly because the
  authors don't know if it is even possible to introduce them in
  `FEniCS`.
[^2]: Indeed, being able to test this new feature is one of the main reason why an interface to
`FEniCS` has been considered.
[^3]: <https://github.com/thelfer/MFrontGenericInterfaceSupport>
[^4]: <https://bitbucket.org/fenics-apps/fenics-solid-mechanics>
[^5]: \(2D\) simulations is theoretically supported, but the authors
  acknowledges that this is not currently functional.
[^6]: Note that an instance of `MaterialDataManager` keeps a reference
  to the behaviour which has been used for its initialization: the user
  must ensure that this behaviour outlives the instance of the
  `MaterialDataManager`, otherwise memory corruption may occur.
