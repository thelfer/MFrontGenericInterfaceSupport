// Copyright (C) 2013-2017 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.

// This program tests the elastic and plastic load-displacement
// responses for a unit cube in uniaxial tension with an Von Mises
// (J2) model with linear strain hardening. The program will throw an
// error is any problems are detected. It uses P2 elements.

#include <memory>
#include <cstdlib>
#include <dolfin/common/Array.h>
#include <dolfin/fem/assemble.h>
#include <dolfin/fem/DirichletBC.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/SystemAssembler.h>
#include <dolfin/function/Function.h>
#include <dolfin/log/Logger.h>
#include <dolfin/nls/NewtonSolver.h>
#include <dolfin/function/Constant.h>
#include <dolfin/function/Expression.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/generation/UnitCubeMesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/SubDomain.h>

#include "Plas3D.h"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/FEniCS/NonLinearMaterial.hxx"
#include "MGIS/FEniCS/NonLinearMechanicalProblem.hxx"

// Displacement right end
struct Load : public dolfin::Expression {
  Load(const double& t) : dolfin::Expression(3), t(t) {}
  void eval(dolfin::Array<double>& values,
            const dolfin::Array<double>& x) const {
    values[0] = 1.0e5 * t;
    values[1] = 0.0;
    values[2] = 0.0;
  }
  const double& t;
};

// Right boundary (x=1)
struct Right : public dolfin::SubDomain
{
  bool inside(const dolfin::Array<double>& x, bool on_boundary) const {
    return std::abs(x[0] - 1.0) < DOLFIN_EPS;
  }
};

// Left boundary (x=1)
struct Left : public dolfin::SubDomain {
  bool inside(const dolfin::Array<double>& x, bool on_boundary) const
  { return x[0] < DOLFIN_EPS; }
};

// Point x = (0, 0, 0)
struct CornerPoint : public dolfin::SubDomain {
  bool inside(const dolfin::Array<double>& x, bool on_boundary) const
  { return (x[0] < DOLFIN_EPS) && (x[1] < DOLFIN_EPS) && (x[2] < DOLFIN_EPS); }
};


int main(){
  // getting the path to the test library
  auto library = std::getenv("MGIS_TEST_BEHAVIOURS_LIBRARY");
  if(library==nullptr){
    std::exit(EXIT_FAILURE);
  }
  // create mesh
  // auto mesh = std::make_shared<dolfin::UnitCubeMesh>(4, 4, 4);
  auto mesh = std::make_shared<dolfin::UnitCubeMesh>(1, 1, 1);
  // Time parameter
  double t = 0.0;
  // Source term, RHS
  auto f = std::make_shared<dolfin::Constant>(0.0, 0.0, 0.0);
  // function space
  auto V = std::make_shared<Plas3D::FunctionSpace>(mesh);

  // Extract elements for stress and tangent
  std::shared_ptr<const dolfin::FiniteElement> element_t;
  std::shared_ptr<const dolfin::FiniteElement> element_s;
  {
    Plas3D::BilinearForm::CoefficientSpace_t Vt(mesh);
    element_t = Vt.element();
  }

  Plas3D::LinearForm::CoefficientSpace_s Vs(mesh);
  element_s = Vs.element();

  // Create boundary conditions (use SubSpace to apply simply
  // supported BCs)
  auto zero = std::make_shared<dolfin::Constant>(0.0);
  auto zero_vector = std::make_shared<dolfin::Constant>(0.0, 0.0, 0.0);
  auto boundary_load = std::make_shared<Load>(t);
  auto Vx = V->sub(0);

  auto left = std::make_shared<Left>();
  auto corner = std::make_shared<CornerPoint>();
  auto left_bc = std::make_shared<dolfin::DirichletBC>(Vx, zero, left);
  auto corner_bc = std::make_shared<dolfin::DirichletBC>(V, zero_vector, corner,
                                                         "pointwise");
  std::vector<std::shared_ptr<const dolfin::DirichletBC>> bcs;
  bcs.push_back(left_bc);
  bcs.push_back(corner_bc);

  // Mark loading boundary
  Right right;
  auto load_marker = std::make_shared<dolfin::MeshFunction<std::size_t>>(mesh,mesh->topology().dim()-1, 0);
  right.mark(*load_marker, 1);

  // Solution function
  auto u = std::make_shared<dolfin::Function>(V);

  auto b = mgis::behaviour::load(library, "Plasticity",
                                 mgis::behaviour::Hypothesis::TRIDIMENSIONAL);
  auto m = mgis::fenics::NonLinearMaterial(u,element_t,element_s,b);
  setExternalStateVariable(m.s0,"Temperature", 293.15);
  setExternalStateVariable(m.s1,"Temperature", 293.15);
  
  // // Create forms and attach functions
  auto a = std::make_shared<Plas3D::BilinearForm>(V, V);
  a->t =  m.getTangentOperatorFunction();
  a->ds = load_marker;
  auto L = std::make_shared<Plas3D::LinearForm>(V);
  L->f = f;
  L->h = boundary_load;
  L->s = m.getThermodynamicForcesFunction();

  // Create PlasticityProblem
  mgis::fenics::NonLinearMechanicalProblem nonlinear_problem(a, L, u, bcs);

  // // Displacement and load integration functionals
  auto M_d = std::make_shared<Plas3D::Form_M_d>(mesh, u);
  auto M_f = std::make_shared<Plas3D::Form_M_f>(mesh, boundary_load);
  M_d->ds = load_marker;
  M_f->ds = load_marker;

  // Create nonlinear solver and set parameters
  dolfin::NewtonSolver nonlinear_solver;
  nonlinear_solver.parameters["convergence_criterion"] = "incremental";
  nonlinear_solver.parameters["maximum_iterations"] = 50;
  nonlinear_solver.parameters["relative_tolerance"] = 1.0e-6;
  nonlinear_solver.parameters["absolute_tolerance"] = 1.0e-15;

  // Structures to hold load-disp data
  std::vector<double> disp, load;
  disp.push_back(0.0);
  load.push_back(0.0);

  // Solver loop
  mgis::size_type step = 0;
  mgis::size_type steps = 10;
  double dt = 0.001;
  while (step < steps) {
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
    disp.push_back(u_avg);
    const double force = assemble(*M_f);
    load.push_back(force);
  }

  for (mgis::size_type i = 0; i != disp.size(); ++i) {
    std::cout << disp[i] << " " << load[i] << '\n';
  }

  return 0;
}
