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
 *
 * \note This file contains material strongly inspired by the
 * fenics-solid-materials by Kristian Oelgaard and Garth N. Wells.
 * See <https://bitbucket.org/fenics-apps/fenics-solid-mechanics> for
 * details.
 */

#include <dolfin/mesh/Cell.h>
#include <dolfin/fem/FiniteElement.h>
#include "MGIS/FEniCS/NonLinearMaterial.hxx"
#include "MGIS/FEniCS/NonLinearMaterialThermodynamicForcesFunction.hxx"

namespace mgis{
  
  namespace fenics {

    void NonLinearMaterialThermodynamicForcesFunction::restrict(
        double* const values,
        const dolfin::FiniteElement&,
        const dolfin::Cell& c,
        const double*,
        const ufc::cell&) const {
      const std::size_t cell_index = c.index();
      const std::size_t num_ip_dofs = this->elements->value_dimension(0);
      const std::size_t num_ip_per_cell =
          this->elements->space_dimension() / num_ip_dofs;
      const auto ths = this->m.s1.thermodynamic_forces_stride;
      dolfin_assert(num_ip_dofs == ths);
      for (std::size_t ip = 0; ip != num_ip_per_cell; ip++) {
        const auto& f = this->m.s1.thermodynamic_forces.data() +
                        ths * (cell_index * num_ip_per_cell + ip);
        for (mgis::size_type i = 0; i != ths; ++i) {
          values[num_ip_per_cell * i + ip] = f[i];
        }
      }
    } // end of NonLinearMaterialThermodynamicForcesFunction::restrict

    NonLinearMaterialThermodynamicForcesFunction::~NonLinearMaterialThermodynamicForcesFunction() = default;
    
  }  // end of namespace fenics

}  // end of namespace mgis
