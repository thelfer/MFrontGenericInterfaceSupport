/*!
 * \file   MGIS/Function/EvaluatorModifierConcept.hxx
 * \brief File which defines the concept an evaluator modifier must satisfy, the
 * `EvaluatorModifiersGroup` class used to combine modifiers and the `operator|`
 * used to apply them to an evaluator.
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
   * \brief a small class to specify that a type is a modifier
   *
   * This is required because the following definition of the
   * `EvaluatorModifierConcept` concept is ill-formed:
   *
   * \code{.cpp}
   * template<typename EvaluatorModifierType>
   * concept EvaluatorModifierConcept = requires(const EvaluatorConcept auto& e,
   *                                             const EvaluatorModifierType& m)
   * { { m(e) } -> EvaluatorConcept;
   *  });
   * \endcode
   *
   * On the other side, defining the `EvaluatorModifierConcept` concept is
   * required to combine modifiers, so relying on `EvaluatorModifierTag` is
   * regarded as a necessary evil.
   */
  struct EvaluatorModifierTag {};

  /*!
   * \brief concept defining the modifiers which can be applied to evaluators
   *
   * As explained in the description of the `EvaluatorModifierTag`, the
   * following definition of the `EvaluatorModifierConcept` concept is
   * ill-formed:
   *
   * \code{.cpp}
   * template<typename EvaluatorModifierType>
   * concept EvaluatorModifierConcept = requires(const EvaluatorConcept auto& e,
   *                                             const EvaluatorModifierType& m)
   * { { m(e) } -> EvaluatorConcept;
   *  });
   * \endcode
   *
   * When exposing an alias named `Tag` equal to the `EvaluatorModifierTag`, a
   * class explicitly states that its call operator can transform any evaluator
   * into a new evaluator.
   *
   * \tparam EvaluatorModifierType: type of the modifier
   */
  template <typename EvaluatorModifierType>
  concept EvaluatorModifierConcept =
      std::same_as<typename EvaluatorModifierType::Tag, EvaluatorModifierTag>;

  /*!
   * \brief applies a modifier to an evaluator
   *
   * \tparam EvaluatorType: evaluator type
   * \tparam EvaluatorModifierType: type of the modifier
   * \param[in] e: evaluator
   * \param[in] m: modifier
   * \return the evaluator resulting from applying the modifier to the evaluator
   */
  template <EvaluatorConcept EvaluatorType,
            EvaluatorModifierConcept EvaluatorModifierType>
  constexpr auto operator|(const EvaluatorType&, EvaluatorModifierType);

  // forward declaration
  template <EvaluatorModifierConcept... Ms>
  requires(sizeof...(Ms) > 0) struct EvaluatorModifiersGroup;

  namespace internals {

    /*!
     * \brief removes the first modifier of a group of modifiers
     *
     * \tparam M0: type of the first modifier, which is removed
     * \tparam Ms: types of the other modifiers
     * \param[in] m: group of modifiers
     * \return a group containing all the modifiers of the given group except
     * the first one
     */
    template <EvaluatorModifierConcept M0, EvaluatorModifierConcept... Ms>
    constexpr auto pop_front(const EvaluatorModifiersGroup<M0, Ms...>& m)  //
        requires(sizeof...(Ms) > 0) {
      return EvaluatorModifiersGroup<Ms...>(std::apply(
          [](auto, auto... rest) { return std::make_tuple(rest...); },
          m.group));
    }  // end of pop_front

    /*!
     * \brief applies a group of modifiers to an evaluator
     *
     * The modifiers are applied in order: the first modifier of the group is
     * applied to the evaluator, the second one to the evaluator resulting from
     * this first modification, and so on.
     *
     * \tparam Ms: types of the grouped modifiers
     * \param[in] m: group of modifiers
     * \param[in] e: evaluator
     * \return the evaluator resulting from applying all the modifiers to the
     * given evaluator
     */
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
   * \brief a class meant to group modifiers, allowing modifiers composability
   *
   * The modifiers are applied from left to right: the first modifier of the
   * group is applied to the evaluator, the second one to the resulting
   * evaluator, and so on.
   *
   * \tparam Ms: types of the grouped modifiers, at least one
   */
  template <EvaluatorModifierConcept... Ms>
  requires(sizeof...(Ms) > 0) struct EvaluatorModifiersGroup {
    //! \brief this alias allows to match the Evaluator Modifier concept
    using Tag = EvaluatorModifierTag;
    /*!
     * \brief constructor from an explicit list of modifiers
     * \param[in] ms: modifiers, which are applied from left to right
     */
    constexpr EvaluatorModifiersGroup(Ms... ms) : group(ms...) {}
    /*!
     * \brief constructor from a tuple of modifiers
     * \param[in] g: tuple of modifiers, which is applied from left to right
     */
    constexpr EvaluatorModifiersGroup(std::tuple<Ms...> g) : group(g) {}

    /*!
     * \brief call operator applying each modifier of the group
     * \param[in] e: evaluator the modifiers of the group are applied to
     * \return the evaluator resulting from applying all the modifiers of the
     * group to the given evaluator
     */
    constexpr auto operator()(const EvaluatorConcept auto& e) const {
      return internals::apply(*this, e);
    }
    //! \brief list of the grouped modifiers
    const std::tuple<Ms...> group;
  };

  /*!
   * \brief deduction guide allowing a group of modifiers to be built from an
   * explicit list of modifiers
   *
   * \tparam Ms: types of the grouped modifiers
   */
  template <EvaluatorModifierConcept... Ms>
  requires(sizeof...(Ms) > 0) EvaluatorModifiersGroup(Ms...)
  ->EvaluatorModifiersGroup<Ms...>;

  /*!
   * \brief deduction guide allowing a group of modifiers to be built from a
   * tuple of modifiers
   *
   * \tparam Ms: types of the grouped modifiers
   */
  template <EvaluatorModifierConcept... Ms>
  requires(sizeof...(Ms) > 0) EvaluatorModifiersGroup(std::tuple<Ms...>)
  ->EvaluatorModifiersGroup<Ms...>;

  /*!
   * \brief groups two modifiers into a new group
   *
   * \tparam M0: type of the first modifier, which is applied first
   * \tparam M1: type of the second modifier, which is applied to the result of
   * the first modification
   * \param[in] m0: first modifier
   * \param[in] m1: second modifier
   * \return a group containing the two given modifiers
   */
  template <EvaluatorModifierConcept M0, EvaluatorModifierConcept M1>
  constexpr auto operator|(M0 const& m0, M1 const& m1) {
    return EvaluatorModifiersGroup<M0, M1>{m0, m1};
  }

  /*!
   * \brief adds a modifier at the beginning of a group of modifiers
   *
   * \tparam M0: type of the modifier, which is applied first
   * \tparam M1: types of the modifiers of the group, which are applied after
   * the given modifier
   * \param[in] m0: modifier
   * \param[in] m1: group of modifiers
   * \return a new group of modifiers
   */
  template <EvaluatorModifierConcept M0, EvaluatorModifierConcept... M1>
  constexpr auto operator|(M0 const& m0,
                           EvaluatorModifiersGroup<M1...> const& m1) {
    return std::apply(
        [&](auto... m) {
          return EvaluatorModifiersGroup{m0, m...};
        },
        m1.group);
  }

  /*!
   * \brief adds a modifier at the end of a group of modifiers
   *
   * \tparam M0: type of the modifier, which is applied last
   * \tparam M1: types of the modifiers of the group, which are applied before
   * the given modifier
   * \param[in] m1: group of modifiers
   * \param[in] m0: modifier
   * \return a new group of modifiers
   */
  template <EvaluatorModifierConcept M0, EvaluatorModifierConcept... M1>
  constexpr auto operator|(EvaluatorModifiersGroup<M1...> const& m1,
                           M0 const& m0) {
    return std::apply(
        [&](auto... m) {
          return EvaluatorModifiersGroup{m..., m0};
        },
        m1.group);
  }

  /*!
   * \brief merges two groups of modifiers
   *
   * \tparam M0: types of the modifiers of the first group, which are applied
   * first
   * \tparam M1: types of the modifiers of the second group, which are applied
   * after the modifiers of the first group
   * \param[in] m0: first group of modifiers
   * \param[in] m1: second group of modifiers
   * \return a new group containing all the modifiers of the two given groups
   */
  template <EvaluatorModifierConcept... M0, EvaluatorModifierConcept... M1>
  constexpr auto operator|(EvaluatorModifiersGroup<M0...> const& m0,
                           EvaluatorModifiersGroup<M1...> const& m1) {
    return EvaluatorModifiersGroup{std::tuple_cat(m0.group, m1.group)};
  }

}  // namespace mgis::function

#include "MGIS/Function/EvaluatorModifierConcept.ixx"

#endif /* LIB_MGIS_FUNCTION_MODIFIERCONCEPT_HXX */
