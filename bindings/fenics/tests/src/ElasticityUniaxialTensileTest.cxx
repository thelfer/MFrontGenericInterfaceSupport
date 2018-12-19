/*!
 * \file   bindings/fenics/tests/ElasticityUniaxialTensileTest.cxx
 * \brief  This program tests the elastic and plastic load-displacement
 * responses for a unit cube in uniaxial tension with an Von Mises
 * (J2) plastic behaviour with linear strain hardening.
 * \author Thomas Helfer
 * \date   14/12/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 *
 * \note This file contains material strongly inspired by the
 * fenics-solid-materials by Kristian Oelgaard and Garth N. Wells.
 * See <https://bitbucket.org/fenics-apps/fenics-solid-mechanics> for
 * details.
 */

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

#include "MGISSmallStrainFormulation3D.h"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/FEniCS/NonLinearMaterial.hxx"
#include "MGIS/FEniCS/NonLinearMechanicalProblem.hxx"
#include "MGIS/FEniCS/FEniCSTestingUtilities.hxx"

// force at the right end
struct Load : public dolfin::Expression {
  Load(const double& t) : dolfin::Expression(3), t(t) {}
  void eval(dolfin::Array<double>& values,
	    const dolfin::Array<double>& x) const {
    values[0] = 100e6 * t;
    values[1] = values[2] = 0.0;
  }
private:
  const double& t;
};


int main(){
  // getting the path to the test library
  auto library = std::getenv("MGIS_TEST_BEHAVIOURS_LIBRARY");
  if(library==nullptr){
    std::exit(EXIT_FAILURE);
  }
  // create mesh and boundaries
  auto mesh = std::make_shared<dolfin::UnitCubeMesh>(4, 4, 4);
  auto boundaries = mgis::fenics::getUnitCubeBoundaries();
  // Time parameter
  double t = 0.0;
  // Source term, RHS
  auto f = std::make_shared<dolfin::Constant>(0.0, 0.0, 0.0);
  // function space
  auto V = std::make_shared<MGISSmallStrainFormulation3D::FunctionSpace>(mesh);

  // Extract elements for stress and tangent
  std::shared_ptr<const dolfin::FiniteElement> element_t;
  std::shared_ptr<const dolfin::FiniteElement> element_s;
  {
    MGISSmallStrainFormulation3D::BilinearForm::CoefficientSpace_t Vt(mesh);
    element_t = Vt.element();
  }

  MGISSmallStrainFormulation3D::LinearForm::CoefficientSpace_s Vs(mesh);
  element_s = Vs.element();

  // Create boundary conditions (use SubSpace to apply simply
  // supported BCs)
  auto zero = std::make_shared<dolfin::Constant>(0.0);
  auto boundary_load = std::make_shared<Load>(t);

  std::vector<std::shared_ptr<const dolfin::DirichletBC>> bcs;
  bcs.push_back(std::make_shared<dolfin::DirichletBC>(V->sub(0), zero, boundaries["sx1"]));
  bcs.push_back(std::make_shared<dolfin::DirichletBC>(V->sub(1), zero, boundaries["sy1"]));
  bcs.push_back(std::make_shared<dolfin::DirichletBC>(V->sub(2), zero, boundaries["sz1"]));

  // Mark loading boundary
  auto load_marker = std::make_shared<dolfin::MeshFunction<std::size_t>>(mesh,mesh->topology().dim()-1, 0);
  boundaries["sx2"]->mark(*load_marker, 1);

  // Solution function
  auto u = std::make_shared<dolfin::Function>(V);

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
  
  // // Create forms and attach functions
  auto a = std::make_shared<MGISSmallStrainFormulation3D::BilinearForm>(V, V);
  a->t =  m.getTangentOperatorFunction();
  a->ds = load_marker;
  auto L = std::make_shared<MGISSmallStrainFormulation3D::LinearForm>(V);
  L->f = f;
  L->h = boundary_load;
  L->s = m.getThermodynamicForcesFunction();

  // create non linear material problem
  mgis::fenics::NonLinearMechanicalProblem nonlinear_problem(a, L, u, bcs);

  // // Displacement and load integration functionals
  auto M_d = std::make_shared<MGISSmallStrainFormulation3D::Form_M_d>(mesh, u);
  auto M_f = std::make_shared<MGISSmallStrainFormulation3D::Form_M_f>(mesh, boundary_load);
  M_d->ds = load_marker;
  M_f->ds = load_marker;

  // Create nonlinear solver and set parameters
  dolfin::NewtonSolver nonlinear_solver;
  nonlinear_solver.parameters["convergence_criterion"] = "incremental";
  nonlinear_solver.parameters["maximum_iterations"] = 50;
  nonlinear_solver.parameters["relative_tolerance"] = 1.0e-6;
  nonlinear_solver.parameters["absolute_tolerance"] = 1.0e-15;

  // post-processings data
  std::vector<double> ux, fx, sxx, exx, eyy;
  ux.push_back(0.0);
  fx.push_back(0.0);
  sxx.push_back(0.0);
  exx.push_back(0.0);
  eyy.push_back(0.0);

  // Solver loop
  mgis::size_type step = 0;
  mgis::size_type steps = 10;
  double dt = 0.1;
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
    ux.push_back(u_avg);
    const double force = assemble(*M_f);
    fx.push_back(force);
    sxx.push_back(m.s1.thermodynamic_forces[0]);
    exx.push_back(m.s1.gradients[0]);
    eyy.push_back(m.s1.gradients[1]);
  }

  auto status = EXIT_SUCCESS;
  auto nb_tests    = mgis::size_type{};
  auto nb_failures = mgis::size_type{};
  auto check = [&status,&nb_tests,&nb_failures](const bool c, const mgis::string_view msg){
    ++nb_tests;
    if(!c){
      std::cerr << msg << '\n';
      status = EXIT_FAILURE;
      ++nb_failures;
    }
  };
  for (mgis::size_type i = 0; i != ux.size(); ++i) {
    check(std::abs(ux[i]-exx[i])<1.e-14, "invalid axial strain value");
    check(std::abs(fx[i]-sxx[i])<1.e-14*yg, "invalid axial stress value");
    check(std::abs(sxx[i]-yg*exx[i])<1.e-14*yg, "invalid axial stress value");
    check(std::abs(eyy[i]+nu*exx[i])<1.e-14, "invalid orthoradial strain value");
    //    std::cout << u[i] << " " << f[i] << '\n';
  }
  // reporting
  std::cout << "Number of tests: " << nb_tests << '\n';
  std::cout << "Number of failures: " << nb_failures << '\n';
  return status;
}
