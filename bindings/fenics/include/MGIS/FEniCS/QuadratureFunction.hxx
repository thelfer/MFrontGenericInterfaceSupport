/*!
 * \file   bindings/fencis/include/MGIS/FEniCS/QuadratureFunction.hxx
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

#ifndef LIB_MGIS_FENICS_QUADRATUREFUNCTION_HXX
#define LIB_MGIS_FENICS_QUADRATUREFUNCTION_HXX

#include <memory>
#include <vector>

namespace dolfin {
  class Cell;
  class FiniteElement;
  class Function;
  class Mesh;
} // end of namespace dolphin

namespace mgis{

  namespace fenics {

    struct QuadratureFunction : public dolfin::GenericFunction {
      /// Constructor
      QuadratureFunction(std::shared_ptr<const dolfin::Mesh> mesh,
                         std::shared_ptr<const dolfin::FiniteElement> element,
                         const std::vector<double>& w);

      /// Constructor
      QuadratureFunction(std::shared_ptr<const dolfin::Mesh> mesh,
                         std::shared_ptr<const dolfin::FiniteElement> element,
                         std::shared_ptr<StateUpdate> state_updater,
                         const std::vector<double>& w);

      std::shared_ptr<const dolfin::FunctionSpace> function_space() const override;

      std::size_t value_rank() const override;

      std::size_t value_dimension(std::size_t) const override;

      std::vector<std::size_t> value_shape() const override;

      void restrict(double*,
                    const dolfin::FiniteElement&,
                    const dolfin::Cell&,
                    const double*,
                    const ufc::cell&) const override;

      void compute_vertex_values(std::vector<double>&,
                                 const dolfin::Mesh&) const override;

      std::shared_ptr<const dolfin::FiniteElement> element() const override;

      //! \brief destructor
      ~QuadratureFunction() override;

     private:

      // Finite element
      std::shared_ptr<const dolfin::FiniteElement> _element;
      // disallow copying
      QuadratureFunction(const QuadratureFunction&) = delete;
      /// delete copy constructor and assignement
       QuadratureFunction(QuadratureFunction&&) = delete;
      // deleting assignment operator
      QuadratureFunction& operator=(const QuadratureFunction&) =
          delete;
      // deleting assignment operator
      QuadratureFunction& operator=(QuadratureFunction&&) = delete;

    }; // end of struct QuadratureFunction

  }  // end of namespace fenics

}  // end of namespace mgis

#endif /* LIB_MGIS_FENICS_QUADRATUREFUNCTION_HXX */
