/*!
 * \file   MGIS/Function/EvaluatorModifierConcept.hxx
 * \brief
 * \author Thomas Helfer
 * \date   05/06/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_EVALUATORMODIFIERCONCEPT_HXX
#define LIB_MGIS_FUNCTION_EVALUATORMODIFIERCONCEPT_HXX

#include <tuple>
#include <concepts>
#include "MGIS/Function/EvaluatorConcept.hxx"

namespace mgis::function {

  /*!
   * \brief a small class to specify that a type is a Modifier
   *
   * This is required because such the following definition
   * of the EvaluatorModifierConcept concept is ill-formed:
   *
   * \code{.cpp}
   * template<typename EvaluatorModifierType>
   * concept EvaluatorModifierConcept = requires(const EvaluatorConcept auto& e,
   *                                             const EvaluatorModifierType& m)
   * { { m(e) } -> EvaluatorConcept;
   *  });
   * \endcode
   *
   * On the other side, defining the EvaluatorModifierConcept concept is
   * required to combine modifiers, so relying on EvaluatorModifierTag is
   * regarded as as a necessary evil.
   */
  struct EvaluatorModifierTag {};

  /*!
   * As explained in the description of the `EvaluatorModifierTag`,
   * the following definition of the EvaluatorModifierConcept concept
   * is ill-formed:
   *
   * \code{.cpp}
   * template<typename EvaluatorModifierType>
   * concept EvaluatorModifierConcept = requires(const EvaluatorConcept auto& e,
   *                                             const EvaluatorModifierType& m)
   * { { m(e) } -> EvaluatorConcept;
   *  });
   * \endcode
   *
   * When exposing an alias named Tag to be equal to the EvaluatorTag, a
   * class explictly states that it has a call operator able to transform
   * any evaluator to a new evaluator.
   */
  template <typename EvaluatorModifierType>
  concept EvaluatorModifierConcept =
      std::same_as<typename EvaluatorModifierType::Tag, EvaluatorModifierTag>;

  /*!
   * \return the evaluator resulting from appling the modifier to the evaluator
   * \param[in] e: evaluator
   * \param[in] m: modifier
   */
  template <EvaluatorConcept EvaluatorType,
            EvaluatorModifierConcept EvaluatorModifierType>
  constexpr auto operator|(const EvaluatorType&, EvaluatorModifierType);

  // forward declaration
  template <EvaluatorModifierConcept... Ms>
  requires(sizeof...(Ms) > 0) struct EvaluatorModifiersGroup;

  namespace internals {

    template <EvaluatorModifierConcept M0, EvaluatorModifierConcept... Ms>
    constexpr auto pop_front(const EvaluatorModifiersGroup<M0, Ms...>& m)  //
        requires(sizeof...(Ms) > 0) {
      return EvaluatorModifiersGroup<Ms...>(std::apply(
          [](auto, auto... rest) { return std::make_tuple(rest...); },
          m.group));
    }  // end of pop_front

    template <EvaluatorModifierConcept... Ms>
    constexpr auto apply(const EvaluatorModifiersGroup<Ms...>& m,
                         const EvaluatorConcept auto& e)  //
        requires(sizeof...(Ms) > 0) {
      if constexpr (sizeof...(Ms) == 1) {
        return std::get<0>(m.group)(e);
      } else {
        const auto ne = std::get<0>(m.group)(e);
        const auto nm = pop_front(m);
        return apply(nm, ne);
      }
    }  // end of apply

  }  // end of namespace internals

  /*!
   * \brief a class meant to group modifiers allowing modifiers
   * composability
   */
  template <EvaluatorModifierConcept... Ms>
  requires(sizeof...(Ms) > 0) struct EvaluatorModifiersGroup {
    //! \brief this alias allows to match the Evaluator Modifier concept
    using Tag = EvaluatorModifierTag;
    /*!
     * \brief constructor from an explicit list of modifiers
     */
    constexpr EvaluatorModifiersGroup(Ms... ms) : group(ms...) {}
    /*!
     * \brief constructor from a tuple of modifiers
     */
    constexpr EvaluatorModifiersGroup(std::tuple<Ms...> g) : group(g) {}

    /*!
     * \brief call operator applying each modifiers of the group
     * \param [in] e: evaluator
     */
    constexpr auto operator()(const EvaluatorConcept auto& e) const {
      return internals::apply(*this, e);
    }
    // list of grouped modifiers
    const std::tuple<Ms...> group;
  };

  template <EvaluatorModifierConcept... Ms>
  requires(sizeof...(Ms) > 0) EvaluatorModifiersGroup(Ms...)
  ->EvaluatorModifiersGroup<Ms...>;

  template <EvaluatorModifierConcept... Ms>
  requires(sizeof...(Ms) > 0) EvaluatorModifiersGroup(std::tuple<Ms...>)
  ->EvaluatorModifiersGroup<Ms...>;

  //! \brief group two modifiers
  template <EvaluatorModifierConcept M0, EvaluatorModifierConcept M1>
  constexpr auto operator|(M0 const& m0, M1 const& m1) {
    return EvaluatorModifiersGroup<M0, M1>{m0, m1};
  }

  template <EvaluatorModifierConcept M0, EvaluatorModifierConcept... M1>
  constexpr auto operator|(M0 const& m0,
                           EvaluatorModifiersGroup<M1...> const& m1) {
    return std::apply(
        [&](auto... m) {
          return EvaluatorModifiersGroup{m0, m...};
        },
        m1.group);
  }

  template <EvaluatorModifierConcept M0, EvaluatorModifierConcept... M1>
  constexpr auto operator|(EvaluatorModifiersGroup<M1...> const& m1,
                           M0 const& m0) {
    return std::apply(
        [&](auto... m) {
          return EvaluatorModifiersGroup{m..., m0};
        },
        m1.group);
  }

  template <EvaluatorModifierConcept... M0, EvaluatorModifierConcept... M1>
  constexpr auto operator|(EvaluatorModifiersGroup<M0...> const& m0,
                           EvaluatorModifiersGroup<M1...> const& m1) {
    return EvaluatorModifiersGroup{std::tuple_cat(m0.group, m1.group)};
  }

}  // namespace mgis::function

#include "MGIS/Function/EvaluatorModifierConcept.ixx"

#endif /* LIB_MGIS_FUNCTION_MODIFIERCONCEPT_HXX */
