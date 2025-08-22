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

#include <ostream>
#include <algorithm>
#include "MGIS/Utilities/Markdown.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/BehaviourDataView.hxx"
#include "MGIS/Behaviour/BehaviourData.hxx"

namespace mgis::behaviour {

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
      : dt(0), rdt(1), speed_of_sound(0), s0(b), s1(s0) {
    static char error_message_buffer[512];
    this->error_message = error_message_buffer;
    this->K.resize(
        std::max(getTangentOperatorArraySize(b), mgis::size_type{1}));
  }  // end of Behaviour::Behaviour

  void update(BehaviourData& d) {
    std::fill(d.K.begin(), d.K.end(), real{0});
    d.rdt = 1;
    d.s0 = d.s1;
  }  // end of update

  void revert(BehaviourData& d) {
    std::fill(d.K.begin(), d.K.end(), real{0});
    d.rdt = 1;
    d.s1 = d.s0;
  }  // end of update

  void print_markdown(std::ostream& os,
                      const Behaviour& b,
                      const BehaviourData& d,
                      const mgis::size_type l) {
    os << mgis::utilities::get_heading_signs(l + 1)
       << " Behaviour description\n\n";
    print_markdown(os, b, l + 1);
    os << mgis::utilities::get_heading_signs(l + 1)
       << " State at the beginning of the time step\n\n";
    print_markdown(os, b, d.s0, l + 1);
    os << mgis::utilities::get_heading_signs(l + 1)
       << " State at the end of the time step\n\n";
    print_markdown(os, b, d.s1, l + 1);
  }  // end of print_markdown

}  // end of namespace mgis::behaviour
