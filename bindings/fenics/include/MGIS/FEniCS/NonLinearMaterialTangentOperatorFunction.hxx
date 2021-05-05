/*!
 * \file
 * bindings/fencis/include/MGIS/FEniCS/NonLinearMaterialTangentOperatorFunction.hxx
 * \brief
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

#ifndef LIB_MGIS_FENICS_NONLINEARMATERIALTANGENTOPERATORFUNCTION_HXX
#define LIB_MGIS_FENICS_NONLINEARMATERIALTANGENTOPERATORFUNCTION_HXX

#include <memory>
#include <vector>
#include "MGIS/FEniCS/Config-FEniCS.hxx"
#include "MGIS/FEniCS/NonLinearMaterialFunctionBase.hxx"

namespace mgis {

  namespace fenics {

    //! \brief function in charge of the export the flux to `FEniCS`.
    struct MGIS_FENICS_EXPORT NonLinearMaterialTangentOperatorFunction
        : NonLinearMaterialFunctionBase {
      // inheriting NonLinearMaterialFunctionBase constructor
      using NonLinearMaterialFunctionBase::NonLinearMaterialFunctionBase;

      void restrict(double*,
                    const dolfin::FiniteElement&,
                    const dolfin::Cell&,
                    const double*,
                    const ufc::cell&) const override;
      //! destructor
      ~NonLinearMaterialTangentOperatorFunction() override;

     private:
      NonLinearMaterialTangentOperatorFunction() = delete;
      NonLinearMaterialTangentOperatorFunction(
          NonLinearMaterialTangentOperatorFunction&&) = delete;
      NonLinearMaterialTangentOperatorFunction(
          const NonLinearMaterialTangentOperatorFunction&) = delete;
      NonLinearMaterialTangentOperatorFunction& operator=(
          NonLinearMaterialTangentOperatorFunction&&) = delete;
      NonLinearMaterialTangentOperatorFunction& operator=(
          const NonLinearMaterialTangentOperatorFunction&) = delete;
    };  // end of struct NonLinearMaterialFunctionFunction

  }  // end of namespace fenics

}  // end of namespace mgis

#endif /* LIB_MGIS_FENICS_NONLINEARMATERIALTANGENTOPERATORFUNCTION_HXX */
