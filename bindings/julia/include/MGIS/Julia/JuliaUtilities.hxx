/*!
 * \file   JuliaUtilities.hxx
 * \brief
 * \author Thomas Helfer
 * \date   16/05/2019
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_JULIA_JULIAUTILITIES_HXX
#define LIB_MGIS_JULIA_JULIAUTILITIES_HXX

#include <vector>
#include <jlcxx/jlcxx.hpp>
#include <jlcxx/array.hpp>
#include <jlcxx/const_array.hpp>
#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/Variable.hxx"
#include "MGIS/Behaviour/Hypothesis.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/Integrate.hxx"

namespace jlcxx {

  //!
  template <>
  struct IsBits<mgis::behaviour::Variable::Type> : std::true_type {};
  //!
  template <>
  struct IsBits<mgis::behaviour::Hypothesis> : std::true_type {};
  //!
  template <>
  struct IsBits<mgis::behaviour::Behaviour::Symmetry> : std::true_type {};
  //!
  template <>
  struct IsBits<mgis::behaviour::Behaviour::BehaviourType> : std::true_type {};
  //!
  template <>
  struct IsBits<mgis::behaviour::Behaviour::Kinematic> : std::true_type {};
  //!
  template <>
  struct IsBits<mgis::behaviour::FiniteStrainBehaviourOptions::StressMeasure>
      : std::true_type {};
  //!
  template <>
  struct IsBits<mgis::behaviour::FiniteStrainBehaviourOptions::TangentOperator>
      : std::true_type {};
  //! \brief
  template <>
  struct IsBits<mgis::behaviour::IntegrationType> : std::true_type {};
}  // end of namespace jlcxx

namespace mgis {

  namespace julia {

    /*!
     * \brief expose std::vector<T>
     * \tparam T: underlying type
     * \param[in,out] m: julia module
     * \param[in] n: name of the exposed vector in Julia
     * \note: the type T must be exposed in Julia
     */
    template <typename T>
    void expose_std_vector(jlcxx::Module&, const char* const);

    template <typename T>
    jlcxx::ConstArray<T, 1u> make_view(const std::vector<T>&);

    template <typename T>
    void assign(std::vector<T>&, const jlcxx::ArrayRef<T>&);

  }  // end of namespace julia

}  // end of namespace mgis

#include "MGIS/Julia/JuliaUtilities.ixx"

#endif /* LIB_MGIS_JULIA_JULIAUTILITIES_HXX */
