/*!
 * \file
 * bindings/fenics/tests/ElasticityUniaxialTensileTestImposedStrain3D-exy.cxx
 * \brief  This program tests the elastic response of an unit cube.
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
struct ImposedDisplacementValue : public dolfin::Expression {
  ImposedDisplacementValue(const double& t_) : dolfin::Expression(), t(t_) {}
  void eval(Eigen::Ref<Eigen::VectorXd> values,
            Eigen::Ref<const Eigen::VectorXd>) const override {
    values[0] = 1.e-3 * t;
  }
  ~ImposedDisplacementValue() override = default;

 private:
  const double& t;
};

int main() {
  // getting the path to the test library
  auto library = std::getenv("MGIS_TEST_BEHAVIOURS_LIBRARY");
  if (library == nullptr) {
    std::exit(EXIT_FAILURE);
  }
  // create mesh and boundaries
  auto mesh = std::make_shared<dolfin::UnitCubeMesh>(1, 1, 1);
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

  // boundary conditions
  auto zero = std::make_shared<dolfin::Constant>(0.0);
  std::vector<std::shared_ptr<const dolfin::DirichletBC>> bcs;
  bcs.push_back(std::make_shared<dolfin::DirichletBC>(V->sub(0), zero,
                                                      boundaries["sy1"]));
  bcs.push_back(std::make_shared<dolfin::DirichletBC>(V->sub(1), zero,
                                                      boundaries["sy1"]));
  bcs.push_back(std::make_shared<dolfin::DirichletBC>(V->sub(2), zero,
                                                      boundaries["sy1"]));
  bcs.push_back(std::make_shared<dolfin::DirichletBC>(
      V->sub(0), std::make_shared<ImposedDisplacementValue>(t),
      boundaries["sy2"]));
  bcs.push_back(std::make_shared<dolfin::DirichletBC>(V->sub(1), zero,
                                                      boundaries["sy2"]));
  bcs.push_back(std::make_shared<dolfin::DirichletBC>(V->sub(2), zero,
                                                      boundaries["sy2"]));

  // Solution function
  auto u = std::make_shared<dolfin::Function>(V);

  auto b = mgis::behaviour::load(library, "Elasticity",
                                 mgis::behaviour::Hypothesis::TRIDIMENSIONAL);
  auto m = mgis::fenics::NonLinearMaterial(u, element_t, element_s, b);
  constexpr const auto yg = 150e9;
  constexpr const auto nu = 0.3;
  setMaterialProperty(m.s0, "YoungModulus", yg);
  setMaterialProperty(m.s1, "YoungModulus", yg);
  setMaterialProperty(m.s0, "PoissonRatio", nu);
  setMaterialProperty(m.s1, "PoissonRatio", nu);
  setExternalStateVariable(m.s0, "Temperature", 293.15);
  setExternalStateVariable(m.s1, "Temperature", 293.15);

  // // Create forms and attach functions
  auto a = std::make_shared<MGISSmallStrainFormulation3D::BilinearForm>(V, V);
  a->t = m.getTangentOperatorFunction();
  auto L = std::make_shared<MGISSmallStrainFormulation3D::LinearForm>(V);
  L->f = std::make_shared<dolfin::Constant>(0.0, 0.0, 0.0);
  L->h = std::make_shared<dolfin::Constant>(0.0, 0.0, 0.0);
  L->s = m.getThermodynamicForcesFunction();

  // create non linear material problem
  mgis::fenics::NonLinearMechanicalProblem nonlinear_problem(a, L, u, bcs);

  // Create nonlinear solver and set parameters
  dolfin::NewtonSolver nonlinear_solver;
  nonlinear_solver.parameters["convergence_criterion"] = "incremental";
  nonlinear_solver.parameters["maximum_iterations"] = 50;
  nonlinear_solver.parameters["relative_tolerance"] = 1.0e-6;
  nonlinear_solver.parameters["absolute_tolerance"] = 1.0e-15;

  // post-processings data
  std::vector<std::array<double, 6>> s, e;
  e.push_back({0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
  s.push_back({0.0, 0.0, 0.0, 0.0, 0.0, 0.0});

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
    e.push_back(extract(m.s1.gradients));
    s.push_back(extract(m.s1.thermodynamic_forces));
  }
  /* tests */
  auto status = EXIT_SUCCESS;
  auto nb_tests = mgis::size_type{};
  auto nb_failures = mgis::size_type{};
  auto check = [&status, &nb_tests, &nb_failures](const bool c,
                                                  const mgis::string_view msg) {
    ++nb_tests;
    if (!c) {
      std::cerr << msg << '\n';
      status = EXIT_FAILURE;
      ++nb_failures;
    }
  };
  constexpr const auto eps = 1.e-14;
  constexpr const auto mu = yg / (2 * (1 + nu));
  for (mgis::size_type i = 0; i != e.size(); ++i) {
    // check strain values
    check(std::abs(e[i][0]) < eps, "invalid strain value");
    check(std::abs(e[i][1]) < eps, "invalid strain value");
    check(std::abs(e[i][2]) < eps, "invalid strain value");
    check(std::abs(e[i][4]) < eps, "invalid strain value");
    check(std::abs(e[i][5]) < eps, "invalid strain value");
    // check stress values
    check(std::abs(s[i][0]) < eps * yg, "invalid stress value");
    check(std::abs(s[i][1]) < eps * yg, "invalid stress value");
    check(std::abs(s[i][2]) < eps * yg, "invalid stress value");
    check(std::abs(s[i][4]) < eps * yg, "invalid stress value");
    check(std::abs(s[i][5]) < eps * yg, "invalid stress value");
    // checking behaviour
    check(std::abs(s[i][3] - 2 * mu * e[i][3]) < eps * yg,
          "invalid stress value");
  }
  // reporting
  std::cout << "Number of tests: " << nb_tests << '\n';
  std::cout << "Number of failures: " << nb_failures << '\n';
  return status;
}
