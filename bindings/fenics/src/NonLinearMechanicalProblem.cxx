// Copyright (C) 2006-2011 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.
//
// First added:  2006-11-13
// Last changed: 2011-02-06

#include <dolfin/common/Timer.h>
#include <dolfin/fem/DirichletBC.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/Form.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/fem/SystemAssembler.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/la/GenericLinearSolver.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/nls/NewtonSolver.h>

#include "MGIS/FEniCS/NonLinearMechanicalProblem.hxx"

namespace mgis {

  namespace fenics {

    NonLinearMechanicalProblem::NonLinearMechanicalProblem(
        std::shared_ptr<const dolfin::Form> a,
        std::shared_ptr<const dolfin::Form> L,
        std::shared_ptr<dolfin::Function> u,
        const std::vector<std::shared_ptr<const dolfin::DirichletBC>> bcs)
        : assembler(a, L, bcs), unknowns(u) {}

    void NonLinearMechanicalProblem::F(dolfin::GenericVector&,
                                       const dolfin::GenericVector&) {
    }  // end of NonLinearMechanicalProblem::F

    void NonLinearMechanicalProblem::J(dolfin::GenericMatrix&,
                                       const dolfin::GenericVector&) {
    } // end of NonLinearMechanicalProblem::J

    void NonLinearMechanicalProblem::form(dolfin::GenericMatrix& A,
                                          dolfin::GenericMatrix&,
                                          dolfin::GenericVector& b,
                                          const dolfin::GenericVector& x) {
      dolfin::Timer timer("NonLinearMechanicalProblem form");
      // Update displacement ghost values
      this->unknowns->update();
      // Build A and b tensors
      this->form_tensors(A, b, x);
    }  // end of NonLinearMechanicalProblem::form

    void NonLinearMechanicalProblem::form_tensors(
        dolfin::GenericMatrix& A,
        dolfin::GenericVector& b,
        const dolfin::GenericVector& x) {
      // Assemble
      dolfin::Timer timer("Assemble non linear mechanical problem LHS/RHS");
      dolfin::set_log_active(false);
      this->assembler.assemble(A, b, x);
      dolfin::set_log_active(true);
    } // end of NonLinearMechanicalProblem::form_tensors

    NonLinearMechanicalProblem::~NonLinearMechanicalProblem() = default;

  }  // end of namespace fenics

}  // end of namespace mgis
