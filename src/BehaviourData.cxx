/*!
 * \file   BehaviourData.cxx
 * \brief
 * \author Thomas Helfer
 * \date   02/08/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include "MGIS/Behaviour/BehaviourData.hxx"

namespace mgis {

  namespace behaviour {

    BehaviourData::BehaviourData() = default;
    BehaviourData::BehaviourData(BehaviourData&&) = default;
    BehaviourData::BehaviourData(const BehaviourData&) = default;
    BehaviourData& BehaviourData::operator=(BehaviourData&&) = default;
    BehaviourData& BehaviourData::operator=(const BehaviourData&) = default;

    // Note: initialiazing s1 with s0 is licit according to the C++ standard:
    //
    // 12.6.2.5
    // Initialization shall proceed in the following order:
    // ...
    // Then, nonstatic data members shall be initialized in the order they were
    // declared in the class definition (again regardless of the order of the
    // mem-initializers).
    BehaviourData::BehaviourData(const Behaviour& b)
        : dt(0), rdt(0), s0(b), s1(s0) {
      constexpr const auto zero = real{0};
      const auto ng = this->s0.gradients.size();
      const auto nth = this->s0.thermodynamic_forces.size();
      this->K.resize(ng * nth, zero);
    }  // end of Behaviour::Behaviour

  }  // end of namespace behaviour

}  // end of namespace mgis
