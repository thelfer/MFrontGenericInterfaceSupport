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

#include <memory>
#include <vector>
#include <utility>
#include <jlcxx/jlcxx.hpp>
#include <jlcxx/array.hpp>
#include <jlcxx/const_array.hpp>
#include <jlcxx/stl.hpp>
#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/Variable.hxx"
#include "MGIS/Behaviour/Hypothesis.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/Integrate.hxx"

namespace jlcxx {

  //!
  template <>
  struct IsMirroredType<mgis::behaviour::Variable::Type> : std::true_type {};
  //!
  template <>
  struct IsMirroredType<mgis::behaviour::Hypothesis> : std::true_type {};
  //!
  template <>
  struct IsMirroredType<mgis::behaviour::Behaviour::Symmetry> : std::true_type {
  };
  //!
  template <>
  struct IsMirroredType<mgis::behaviour::Behaviour::BehaviourType>
      : std::true_type {};
  //!
  template <>
  struct IsMirroredType<mgis::behaviour::Behaviour::Kinematic>
      : std::true_type {};
  //!
  template <>
  struct IsMirroredType<
      mgis::behaviour::FiniteStrainBehaviourOptions::StressMeasure>
      : std::true_type {};
  //!
  template <>
  struct IsMirroredType<
      mgis::behaviour::FiniteStrainBehaviourOptions::TangentOperator>
      : std::true_type {};
  template <>
  struct IsMirroredType<mgis::behaviour::IntegrationType> : std::true_type {};
  template <>
  struct IsMirroredType<mgis::behaviour::BehaviourDataView> : std::false_type {
  };

}  // end of namespace jlcxx

namespace mgis {

  namespace julia {

    template <typename T>
    jlcxx::ConstArray<T, 1u> make_view(const std::vector<T>&);

    template <typename T>
    void assign(std::vector<T>&, const jlcxx::ArrayRef<T>&);

    /*!
     * \brief behaviour returned by `load`. It is shared by the objects built
     * from it (states and behaviour data), which thus keep it alive when the
     * Julia object of the behaviour is collected.
     */
    struct SharedBehaviour {
      //! \brief constructor from a loaded behaviour, which is moved
      explicit SharedBehaviour(mgis::behaviour::Behaviour&& b)
          : behaviour(
                std::make_shared<mgis::behaviour::Behaviour>(std::move(b))) {}
      //! \brief shared behaviour
      std::shared_ptr<mgis::behaviour::Behaviour> behaviour;
    };

    /*!
     * \brief object of type `T` (`State` or `BehaviourData`) built from a
     * shared behaviour, which it keeps alive.
     */
    template <typename T>
    struct Wrapper {
      //! \brief constructor: the other arguments are passed to the one of `T`
      template <typename... Arguments>
      Wrapper(const SharedBehaviour& b, Arguments&&... args)
          : behaviour(b.behaviour),
            value(*behaviour, std::forward<Arguments>(args)...) {}
      //! \brief behaviour, declared first to be destroyed last
      std::shared_ptr<const mgis::behaviour::Behaviour> behaviour;
      //! \brief wrapped object
      T value;
    };

  }  // end of namespace julia

}  // end of namespace mgis

namespace jlcxx {

  /*
   * `SharedBehaviour` and `Wrapper<T>` are seen from Julia as subtypes of
   * `Behaviour` and `T`, so that all the functions taking a `Behaviour` or a
   * `T` accept them. They do not derive from these types: the conversion
   * returns the shared behaviour or the wrapped object, and there is no
   * conversion in the other direction.
   */

  template <>
  struct SuperType<mgis::julia::SharedBehaviour> {
    using type = mgis::behaviour::Behaviour;
  };
  template <>
  struct UpCast<mgis::julia::SharedBehaviour> {
    static mgis::behaviour::Behaviour& apply(mgis::julia::SharedBehaviour& b) {
      return *(b.behaviour);
    }
  };
  template <>
  struct DownCast<mgis::behaviour::Behaviour, mgis::julia::SharedBehaviour> {
    static void apply(Module&) {}
  };
  template <typename T>
  struct SuperType<mgis::julia::Wrapper<T>> {
    using type = T;
  };
  template <typename T>
  struct UpCast<mgis::julia::Wrapper<T>> {
    static T& apply(mgis::julia::Wrapper<T>& w) { return w.value; }
  };
  template <typename T>
  struct DownCast<T, mgis::julia::Wrapper<T>> {
    static void apply(Module&) {}
  };

}  // end of namespace jlcxx

#include "MGIS/Julia/JuliaUtilities.ixx"

#endif /* LIB_MGIS_JULIA_JULIAUTILITIES_HXX */
