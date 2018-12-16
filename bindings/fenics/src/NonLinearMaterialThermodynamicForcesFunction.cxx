/*!
 * \file   bindings/fenics/src/NonLinearMaterialThermodynamicForcesFunction.cxx
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

#include "MGIS/FEniCS/NonLinearMaterial.hxx"
#include "MGIS/FEniCS/NonLinearMaterialThermodynamicForcesFunction.hxx"

namespace mgis{
  
  namespace fenics {

    void NonLinearMaterialThermodynamicForcesFunction::restrict(double* const values,
								const dolfin::FiniteElement&,
								const dolfin::Cell&,
								const double*,
								const ufc::cell&) const{
      const auto& f = this->m.s1.thermodynamic_forces;
      std::copy(f.begin(),f.end(), values);
    } // end of NonLinearMaterialThermodynamicForcesFunction::restrict

    NonLinearMaterialThermodynamicForcesFunction::~NonLinearMaterialThermodynamicForcesFunction() = default;
    
  }  // end of namespace fenics

}  // end of namespace mgis
