/*!
 * \file   bindings/fenics/tests/src/PlasticCylinderExpansion.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   19/12/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
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

#include "MGISSmallStrainFormulation2D.h"
// #include "MGIS/Behaviour/Behaviour.hxx"
// #include "MGIS/FEniCS/NonLinearMaterial.hxx"
// #include "MGIS/FEniCS/NonLinearMechanicalProblem.hxx"
// #include "MGIS/FEniCS/FEniCSTestingUtilities.hxx"

int main(){
  // getting the path to the test library
  auto getenv = [](const char* const n) -> std::string{
    const auto* const e = std::getenv(n);
    if(e==nullptr){
      std::exit(EXIT_FAILURE);
    }
    return e;
  };
  const auto data    = getenv("MGIS_FENICS_TEST_DATA");
  const auto library = getenv("MGIS_TEST_BEHAVIOURS_LIBRARY");
  // create mesh
  auto mesh = std::make_shared<dolfin::Mesh>(data+"/thick_cylinder.xml");
  auto facets = std::make_shared<dolfin::MeshFunction<size_t>>(
      mesh, data+"/thick_cylinder_facet_region.xml");
  // function space
  auto V = std::make_shared<MGISSmallStrainFormulation2D::FunctionSpace>(mesh);
<
  // Extract elements for stress and tangent
  MGISSmallStrainFormulation2D::BilinearForm::CoefficientSpace_t Vt(mesh);
  auto element_t = Vt.element();
  MGISSmallStrainFormulation2D::LinearForm::CoefficientSpace_s Vs(mesh);
  auto element_s = Vs.element();

  // Time parameter
  double t = 0.0;

  // boundary conditions
  // auto zero = std::make_shared<dolfin::Constant>(0.0);
  // std::vector<std::shared_ptr<const dolfin::DirichletBC>> bcs;
  // bcs.push_back(std::make_shared<dolfin::DirichletBC>(V->sub(1), zero, facets,1));
  // bcs.push_back(std::make_shared<dolfin::DirichletBC>(V->sub(0), zero, facets,3));


  return 0;
}
