/*!
 * \file   bindings/fenics/src/NonLinearMaterial.cxx
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

#include <dolfin/la/GenericVector.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/GenericDofMap.h>
#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/Integrate.hxx"
#include "MGIS/FEniCS/Utils.hxx"
#include "MGIS/FEniCS/NonLinearMaterial.hxx"

namespace mgis{
  
  namespace fenics {

    static void compute_strain2D(double * const e,
				 const boost::multi_array<double, 2>& dsf,
				 const std::vector<double>& ecs)
    {
      const std::size_t space_dim = dsf.shape()[0];
      e[0] = e[1] = e[2] = e[3] = 0;
      for (unsigned int dim = 0; dim < space_dim; dim++){
	// Ux,x (eps_xx)
	e[0] += dsf[dim][0]*ecs[dim];
	// Uy,y (eps_yy)
	e[1] += dsf[dim][1]*ecs[space_dim + dim];
	// Ux,y + Uy,x (gamma_xy)
	e[3] += (dsf[dim][1]*ecs[dim] +
		 dsf[dim][0]*ecs[space_dim + dim]);
      }
    }  // end of compute_strain2D
    
    static void compute_strain3D(double * const e,
				 const boost::multi_array<double, 2>& dsf,
				 const std::vector<double>& ecs){
      // Zero strain vector
      e[0] = e[1] = e[2] = e[3] = e[4] = e[5] = 0;
      const std::size_t space_dim = dsf.shape()[0];
      for (unsigned int dim = 0; dim < space_dim; dim++){
	// Ux,x (eps_xx)
	e[0] += dsf[dim][0]*ecs[dim];
	// Uy,y (eps_yy)
	e[1] += dsf[dim][1]*ecs[space_dim + dim];
	// Uz,z (eps_zz)
	e[2] += dsf[dim][2]*ecs[2*space_dim + dim];
	// Ux,y + Uy,x (gamma_xy)
	e[3] += (dsf[dim][1]*ecs[dim] +
		 dsf[dim][0]*ecs[space_dim + dim]);
	// Ux,z + Uz,x (gamma_xz)
	e[4] += (dsf[dim][2]*ecs[dim] +
		 dsf[dim][0]*ecs[2*space_dim + dim]);
	// Uy,z + Uz,y (gamma_yz)
	e[5] += (dsf[dim][2]*ecs[space_dim + dim] +
		 dsf[dim][1]*ecs[2*space_dim + dim]);
      }
    } // end of compute_strain3D
    
    NonLinearMaterial::NonLinearMaterial(std::shared_ptr<const dolfin::Function> u,
					 std::shared_ptr<const dolfin::FiniteElement> e,
					 const mgis::behaviour::Behaviour& bv)
      : mgis::behaviour::MaterialDataManager(bv,getNumberOfIntegrationPoints(*e)),
      unknowns(u),
      elements(e),
      dt(0)
    {
      const auto s =
	this->unknowns->function_space()->dofmap()->max_cell_dimension();
      this->expansion_coefficients.resize(s);
      // get stress UFC element
      auto ufc_element_sigma = this->elements->ufc_element();
      dolfin_assert(ufc_element_sigma);
      // Get stress dof dimension data
      const std::size_t dim = ufc_element_sigma->space_dimension();
      //      const std::size_t gdim = ufc_element_sigma->geometric_dimension();
      const std::size_t tdim = ufc_element_sigma->topological_dimension();
      // Get quadrature point coordinates on reference element
      this->ip_coordinates.resize(boost::extents[dim][tdim]);
      ufc_element_sigma->tabulate_reference_dof_coordinates(this->ip_coordinates.data());
      
      // Get displacement UFC element (single component)
      const dolfin::FiniteElement& u_element_new = *(*this->unknowns)[0].function_space()->element();
      auto ufc_element_u = u_element_new.ufc_element();
      dolfin_assert(ufc_element_u);
      // Compute basis function derivatives on reference element and store
      const std::size_t dim_u = ufc_element_u->space_dimension();
      //      const std::size_t gdim_u = ufc_element_u->geometric_dimension();
      const std::size_t num_points = this->ip_coordinates.shape()[0];
      //      const std::size_t value_size = ufc_element_u->reference_value_size();
      
      //boost::multi_array<double, 3> derivatives(boost::extents[num_points][dim_u][tdim]);
      this->dsf.resize(boost::extents[num_points][dim_u][tdim]);
      ufc_element_u->evaluate_reference_basis_derivatives(this->dsf.data(),
							  1, this->ip_coordinates.shape()[0],
							  this->ip_coordinates.data());
    } // end of NonLinearMaterial::NonLinearMaterial

    void NonLinearMaterial::setTimeIncrement(const double ndt){
      this->dt = ndt;
    } // end of NonLinearMaterial::setTimeIncrement

    std::shared_ptr<NonLinearMaterialThermodynamicForcesFunction>
    NonLinearMaterial::getThermodynamicForcesFunction(){
      return std::make_shared<NonLinearMaterialThermodynamicForcesFunction>(*this);
    } // end of NonLinearMaterial::getThermodynamicForcesFunction

    std::shared_ptr<NonLinearMaterialTangentOperatorFunction>
    NonLinearMaterial::getTangentOperatorFunction(){
      return std::make_shared<NonLinearMaterialTangentOperatorFunction>(*this);
    } // end of NonLinearMaterial::getTangentOperatorFunction
    
    void NonLinearMaterial::update_gradients(const dolfin::Cell& c,
					     const double* nc){
      const std::size_t cell_index = c.index();
      // function space
      const auto& fs = *(this->unknowns->function_space());
      // mesh
      const auto& m  = *(fs.mesh());
      // get solution dofs on cell
      const auto& dm = *(fs.dofmap());
      const auto dofs = dm.cell_dofs(cell_index);
      // get expansion coefficients on cell
      dolfin_assert(this->expansion_coefficients.size() ==  dm.max_cell_dimension());
      this->unknowns->vector()->get_local(this->expansion_coefficients.data(),
					  dofs.size(), dofs.data());
      // compute geometry mapping
      const auto gdim = m.geometry().dim();
      double lJ[9], lK[9], lDetJ;
      if (gdim == 2){
	compute_jacobian_triangle_2d(lJ, nc);
	compute_jacobian_inverse_triangle_2d(lK, lDetJ, lJ);
      } else if (gdim == 3){
	compute_jacobian_tetrahedron_3d(lJ, nc);
	compute_jacobian_inverse_tetrahedron_3d(lK, lDetJ, lJ);
      } else {
	mgis::raise("update_gradients: "
		    "unsupported dimension");
      }
      // duplicate data at each point
      auto num_points = this->dsf.shape()[0];
      boost::multi_array<double, 2> gJ(boost::extents[num_points][9]);
      std::vector<double> gDetJ(num_points, lDetJ);
      boost::multi_array<double, 2> gK(boost::extents[num_points][9]);
      for (std::size_t i = 0; i < num_points; ++i){
	for (std::size_t j = 0; j < 9; ++j){
	  gJ[i][j] = lJ[j];
	  gK[i][j] = lK[j];
	}
      }
      // Get displacement UFC element (single component)
      const dolfin::FiniteElement& u_element_new = *(*this->unknowns)[0].function_space()->element();
      auto ufc_element_u = u_element_new.ufc_element();
      dolfin_assert(ufc_element_u);
      // Push derivatives forward to current physical cell
      auto dim_u = this->dsf.shape()[1];
      boost::multi_array<double, 3> derivs_physical(boost::extents[num_points][dim_u][gdim]);
      ufc_element_u->transform_reference_basis_derivatives(derivs_physical.data(), 1, num_points,
							   this->dsf.data(),
							   this->ip_coordinates.data(),
							   gJ.data(), gDetJ.data(), gK.data(), 0);
      // loop over quadrature points
      const std::size_t num_ip_dofs = this->elements->value_dimension(0);
      const std::size_t num_ip_per_cell = this->elements->space_dimension()/num_ip_dofs;
      if (gdim == 2){
	const auto o = cell_index * num_ip_per_cell * 4;
	auto *const e = this->s1.gradients.data() + o;
	for (std::size_t ip = 0; ip < num_ip_per_cell; ++ip){
	  compute_strain2D(e + 4 * ip, derivs_physical[ip],
			   this->expansion_coefficients);
	}
      } else if (gdim == 3){
	const auto o = cell_index * num_ip_per_cell * 6;
	auto *const e = this->s1.gradients.data() + o;
	for (std::size_t ip = 0; ip < num_ip_per_cell; ++ip){
	  compute_strain3D(e + 6 * ip, derivs_physical[ip],
			   this->expansion_coefficients);
	}
      } else {
	mgis::raise("update_gradients: "
		    "unsupported dimension");
      }
    } // end of update_gradients
    
    void NonLinearMaterial::update(const dolfin::Cell& c,
				   const double* nc){
      const std::size_t cell_index = c.index();
      const std::size_t num_ip_dofs = this->elements->value_dimension(0);
      const std::size_t num_ip_per_cell = this->elements->space_dimension()/num_ip_dofs;
      const auto bc = cell_index       * num_ip_per_cell;
      const auto ec = (cell_index + 1) * num_ip_per_cell;
      this->update_gradients(c,nc);
      integrate(*this,this->dt,bc,ec);
    } // end of NonLinearMaterial::update
    
    NonLinearMaterial::~NonLinearMaterial() = default;

  }  // end of namespace fenics

}  // end of namespace mgis
