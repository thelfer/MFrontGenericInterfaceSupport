/*!
 * \file   bindings/fencis/include/MGIS/FEniCS/NonLinearMaterial.hxx
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

#ifndef LIB_MGIS_FENICS_NONLINEARMATERIAL_HXX
#define LIB_MGIS_FENICS_NONLINEARMATERIAL_HXX

#include <memory>
#include <vector>
#include <boost/multi_array.hpp>
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MGIS/FEniCS/Config-FEniCS.hxx"
#include "MGIS/FEniCS/NonLinearMaterialTangentOperatorFunction.hxx"
#include "MGIS/FEniCS/NonLinearMaterialThermodynamicForcesFunction.hxx"

namespace mgis{

  namespace fenics {

    /*!
     * \brief main interface between `MGIS` and `FEniCS`.
     * This design has been inspired by the one of the
     * `fenics-solid-mechanics` project.
     */
    struct MGIS_FENICS_EXPORT NonLinearMaterial
      : public mgis::behaviour::MaterialDataManager {
      /*!
       * \param[in] u: unknowns space
       * \param[in] t: tangent operator finite elements
       * \param[in] e: stress finite elements
       * \param[in] bv: behaviour
       */
      NonLinearMaterial(std::shared_ptr<const dolfin::Function>,
			std::shared_ptr<const dolfin::FiniteElement>,
			std::shared_ptr<const dolfin::FiniteElement>,
			const mgis::behaviour::Behaviour&);
      //! set time increment
      void setTimeIncrement(const double);
      //! 
      void update(const dolfin::Cell&, const double*);
      //! \return a function able to evaluate the thermodynamic forces
      std::shared_ptr<NonLinearMaterialThermodynamicForcesFunction>
      getThermodynamicForcesFunction();
      //! \return a function able to evaluate the tangent operator
      std::shared_ptr<NonLinearMaterialTangentOperatorFunction>
      getTangentOperatorFunction();
      //! \brief displacements unknowns
      std::shared_ptr<const dolfin::Function> unknowns;
      //! \brief destructor
      ~NonLinearMaterial();
    private:
      //! \brief underlying elements for the thermodynamic forces
      std::shared_ptr<const dolfin::FiniteElement> tangent_operator_elements;
      //! \brief underlying elements for the tangent operator
      std::shared_ptr<const dolfin::FiniteElement> thf_elements;
      //! 
      void update_gradients(const dolfin::Cell&,
			    const double*);
      /*!
       * \brief basis function derivatives at integration points on
       * reference element
       */
      boost::multi_array<double, 3> dsf;
      //! \brief integration points in reference coordinates
      boost::multi_array<double, 2> ip_coordinates;
      //! \brief scratch data
      std::vector<double> expansion_coefficients;
      //! \brief current time increment
      double dt;
    }; // end of struct NonLinearMaterial

  }  // end of namespace fenics

}  // end of namespace mgis

#endif /* LIB_MGIS_FENICS_NONLINEARMATERIAL_HXX */
