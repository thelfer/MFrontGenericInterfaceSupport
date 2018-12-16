/*!
 * \file   bindings/fenics/src/NonLinearMaterialFunctionBase.cxx
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

#include <dolfin/fem/FiniteElement.h>
#include "MGIS/FEniCS/NonLinearMaterial.hxx"
#include "MGIS/FEniCS/NonLinearMaterialFunctionBase.hxx"

namespace mgis{
  
  namespace fenics {

    NonLinearMaterialFunctionBase::NonLinearMaterialFunctionBase(NonLinearMaterial& nlm)
      : m(nlm)
    {} // end of NonLinearMaterialFunctionBase::NonLinearMaterialFunctionBase
    
    std::shared_ptr<const dolfin::FunctionSpace>
    NonLinearMaterialFunctionBase::function_space() const {
      return {};
    } // end of NonLinearMaterialFunctionBase::function_space
    
    std::size_t NonLinearMaterialFunctionBase::value_rank() const {
      return this->m.elements->value_rank();
    } // end of NonLinearMaterialFunctionBase::value_rank
    
    std::size_t NonLinearMaterialFunctionBase::value_dimension(const std::size_t i) const {
      return this->m.elements->value_dimension(i);
    } // end of NonLinearMaterialFunctionBase::value_dimension
    
    std::vector<std::size_t> NonLinearMaterialFunctionBase::value_shape() const {
      using size_type = std::vector<std::size_t>::size_type;
      std::vector<std::size_t> shape(this->value_rank(), 1);
      for (size_type i = 0; i != shape.size(); ++i){
        shape[i] = this->value_dimension(i);
      }
      return shape;
    } // end of NonLinearMaterialFunctionBase::value_shape
    
    void NonLinearMaterialFunctionBase::compute_vertex_values(std::vector<double>&,
							       const dolfin::Mesh&) const{
      dolfin::error("NonLinearMaterialFunctionBase::compute_vertex_values"
		    " not implemented");
    } // end of NonLinearMaterialFunctionBase::compute_vertex_values

    NonLinearMaterialFunctionBase::~NonLinearMaterialFunctionBase() = default;
    
  }  // end of namespace fenics

}  // end of namespace mgis
