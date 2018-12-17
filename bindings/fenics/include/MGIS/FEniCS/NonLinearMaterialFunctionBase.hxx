/*!
 * \file   bindings/fencis/include/MGIS/FEniCS/NonLinearMaterialFunctionBase.hxx
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

#ifndef LIB_MGIS_FENICS_NONLINEARMATERIALFUNCTIONBASE_HXX
#define LIB_MGIS_FENICS_NONLINEARMATERIALFUNCTIONBASE_HXX

#include <memory>
#include <vector>
#include <dolfin/function/GenericFunction.h>
#include "MGIS/FEniCS/Config-FEniCS.hxx"

namespace mgis{

  namespace fenics {

    // forward declaration
    struct NonLinearMaterial;
    
    /*!
     * \brief main interface between `MGIS` and `FEniCS`.
     * This design has been inspired by the one of the
     * `fenics-solid-mechanics` project.
     */
    struct MGIS_FENICS_EXPORT NonLinearMaterialFunctionBase
      : public dolfin::GenericFunction {
      /*! 
       * \brief constructor
       * \param[in] nlm: material 
       */
      NonLinearMaterialFunctionBase(NonLinearMaterial&,
				    std::shared_ptr<const dolfin::FiniteElement>);

      std::shared_ptr<const dolfin::FunctionSpace> function_space() const override;

      std::size_t value_rank() const override;

      std::size_t value_dimension(std::size_t) const override;

      std::vector<std::size_t> value_shape() const override;

      void compute_vertex_values(std::vector<double>&,
                                 const dolfin::Mesh&) const override;

      //! \brief destructor
      ~NonLinearMaterialFunctionBase() override;

    protected:
      NonLinearMaterial& m;
      std::shared_ptr<const dolfin::FiniteElement> elements;
    private:
      // disallow copying
      NonLinearMaterialFunctionBase(const NonLinearMaterialFunctionBase&) = delete;
      /// delete copy constructor and assignement
       NonLinearMaterialFunctionBase(NonLinearMaterialFunctionBase&&) = delete;
      // deleting assignment operator
      NonLinearMaterialFunctionBase&
      operator=(const NonLinearMaterialFunctionBase&) = delete;
      // deleting assignment operator
      NonLinearMaterialFunctionBase&
      operator=(NonLinearMaterialFunctionBase&&) = delete;
    }; // end of struct NonLinearMaterialFunctionBase

  }  // end of namespace fenics

}  // end of namespace mgis

#endif /* LIB_MGIS_FENICS_NONLINEARMATERIALFUNCTIONBASE_HXX */
