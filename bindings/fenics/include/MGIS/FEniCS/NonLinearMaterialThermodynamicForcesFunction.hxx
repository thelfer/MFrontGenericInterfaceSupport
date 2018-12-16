/*!
 * \file   bindings/fencis/include/MGIS/FEniCS/NonLinearMaterialThermodynamicForcesFunction.hxx
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
 */

#ifndef LIB_MGIS_FENICS_NONLINEARMATERIALTHERMODYNAMICFORCESFUNCTION_HXX
#define LIB_MGIS_FENICS_NONLINEARMATERIALTHERMODYNAMICFORCESFUNCTION_HXX

#include <memory>
#include <vector>
#include "MGIS/FEniCS/Config-FEniCS.hxx"
#include "MGIS/FEniCS/NonLinearMaterialFunctionBase.hxx"

namespace mgis{
  
  namespace fenics {

    //! \brief function in charge of the export the flux to `FEniCS`.
    struct MGIS_FENICS_EXPORT NonLinearMaterialThermodynamicForcesFunction
      : NonLinearMaterialFunctionBase {
      // inheriting NonLinearMaterialFunctionBase constructor
      using NonLinearMaterialFunctionBase::NonLinearMaterialFunctionBase;

      void restrict(double*,
                    const dolfin::FiniteElement&,
                    const dolfin::Cell&,
                    const double*,
                    const ufc::cell&) const override;
      //! destructor
      ~NonLinearMaterialThermodynamicForcesFunction() override;
    private:
      NonLinearMaterialThermodynamicForcesFunction() = delete;
      NonLinearMaterialThermodynamicForcesFunction(NonLinearMaterialThermodynamicForcesFunction&&) = delete;
      NonLinearMaterialThermodynamicForcesFunction(const NonLinearMaterialThermodynamicForcesFunction&) = delete;
      NonLinearMaterialThermodynamicForcesFunction&
      operator=(NonLinearMaterialThermodynamicForcesFunction&&) = delete;
      NonLinearMaterialThermodynamicForcesFunction&
      operator=(const NonLinearMaterialThermodynamicForcesFunction&) = delete;
    }; // end of struct NonLinearMaterialFunctionFunction

  }  // end of namespace fenics

}  // end of namespace mgis

#endif /* LIB_MGIS_FENICS_NONLINEARMATERIALTHERMODYNAMICFORCESFUNCTION_HXX */
