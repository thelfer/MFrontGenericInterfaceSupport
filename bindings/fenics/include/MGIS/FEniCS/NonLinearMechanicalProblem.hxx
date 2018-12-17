// Copyright (C) 2006-2011 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.
//
// First added:  2006-11-13
// Last changed: 2011-02-06

#ifndef LIB_MGIS_FENICS_NONLINEARMECHANICSPROBLEM_HXX
#define LIB_MGIS_FENICS_NONLINEARMECHANICSPROBLEM_HXX

#include <memory>
#include <vector>
#include <dolfin/nls/NonlinearProblem.h>
#include <dolfin/parameter/Parameters.h>
#include "MGIS/FEniCS/Config-FEniCS.hxx"

namespace mgis {

  namespace fenics {

    struct MGIS_FENICS_EXPORT NonLinearMechanicalProblem
      : public dolfin::NonlinearProblem {
      //! \brief Constructor
      NonLinearMechanicalProblem(
          std::shared_ptr<const dolfin::Form>,
          std::shared_ptr<const dolfin::Form>,
          std::shared_ptr<dolfin::Function>,
          const std::vector<std::shared_ptr<const dolfin::DirichletBC>>);
      //! \brief Loop quadrature points and compute local tangents and stresses
      void form(dolfin::GenericMatrix&,
                dolfin::GenericMatrix&,
                dolfin::GenericVector&,
                const dolfin::GenericVector&) override;
      //! \brief User defined assemble of residual vector
      void F(dolfin::GenericVector&, const dolfin::GenericVector&) override;
      //! \brief User defined assemble of Jacobian matrix
      void J(dolfin::GenericMatrix&, const dolfin::GenericVector&) override;
      //! \brief destructor
      ~NonLinearMechanicalProblem() override;

     private:

      // For system after constitutive update
      void form_tensors(dolfin::GenericMatrix&,
                        dolfin::GenericVector&,
                        const dolfin::GenericVector&);
      // assembler
      dolfin::SystemAssembler assembler;
      // unknowns
      std::shared_ptr<const dolfin::Function> unknowns;
      // disabling default constructors and assignement operators
      NonLinearMechanicalProblem() = delete;
      NonLinearMechanicalProblem(NonLinearMechanicalProblem&&) = delete;
      NonLinearMechanicalProblem(const NonLinearMechanicalProblem&) = delete;
      NonLinearMechanicalProblem& operator=(NonLinearMechanicalProblem&&) =
          delete;
      NonLinearMechanicalProblem& operator=(const NonLinearMechanicalProblem&) =
          delete;
};  // end of struct

}  // end of namespace fenics

}  // end of namespace mgis

#endif /* LIB_MGIS_FENICS_NONLINEARMECHANICSPROBLEM_HXX */
