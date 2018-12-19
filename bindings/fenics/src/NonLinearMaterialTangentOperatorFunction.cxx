/*!
 * \file   bindings/fenics/src/NonLinearMaterialTangentOperatorFunction.cxx
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
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MGIS/Behaviour/Integrate.hxx"
#include "MGIS/FEniCS/NonLinearMaterial.hxx"
#include "MGIS/FEniCS/NonLinearMaterialTangentOperatorFunction.hxx"

namespace mgis{
  
  namespace fenics {

    void NonLinearMaterialTangentOperatorFunction::restrict(
        double* const values,
        const dolfin::FiniteElement&,
        const dolfin::Cell& c,
        const double* nc,
        const ufc::cell&) const {
      // behaviour integration
      this->m.update(c,nc);
      const auto gs  = this->m.s1.gradients_stride;
      const auto ths = this->m.s1.thermodynamic_forces_stride;
      // updating the tangent operator
      const std::size_t cell_index = c.index();
      const std::size_t num_ip_dofs = this->elements->value_dimension(0);
      dolfin_assert(num_ip_dofs == gs * ths);
      const std::size_t num_ip_per_cell =
          this->elements->space_dimension() / num_ip_dofs;
      for (std::size_t ip = 0; ip != num_ip_per_cell; ++ip) {
        const auto gip = cell_index * num_ip_per_cell + ip;
        const auto Kb = this->m.K.data() + gip * num_ip_dofs;
	for (mgis::size_type i = 0; i != ths; ++i) {
          for (mgis::size_type j = 0; j != gs; ++j) {
            const std::size_t pos = i * gs + j;
            values[num_ip_per_cell * pos + ip] = Kb[pos];
          }
        }
      }
    } // end of NonLinearMaterialTangentOperatorFunction::restrict

    NonLinearMaterialTangentOperatorFunction::~NonLinearMaterialTangentOperatorFunction() = default;
    
  }  // end of namespace fenics

}  // end of namespace mgis
