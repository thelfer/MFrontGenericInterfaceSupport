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
#include <utility>
#include <algorithm>
#include "MGIS/Utilities/Markdown.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/BehaviourDataView.hxx"
#include "MGIS/Behaviour/BehaviourData.hxx"

namespace mgis::behaviour {

  /*!
   * \brief make the error message of a copy of the given object point to the
   * buffer of the copy if the given object uses its own buffer
   * \param[out] d: copy
   * \param[in] src: copied object
   */
  static void setErrorMessageBuffer(BehaviourData& d,
                                    const BehaviourData& src) noexcept {
    if (src.error_message == src.error_message_buffer.data()) {
      d.error_message = d.error_message_buffer.data();
    }
  }  // end of setErrorMessageBuffer

  BehaviourData::BehaviourData(BehaviourData&& src)
      : error_message(src.error_message),
        dt(src.dt),
        K(std::move(src.K)),
        rdt(src.rdt),
        speed_of_sound(src.speed_of_sound),
        s0(std::move(src.s0)),
        s1(std::move(src.s1)),
        error_message_buffer(src.error_message_buffer) {
    setErrorMessageBuffer(*this, src);
  }  // end of BehaviourData

  BehaviourData::BehaviourData(const BehaviourData& src)
      : error_message(src.error_message),
        dt(src.dt),
        K(src.K),
        rdt(src.rdt),
        speed_of_sound(src.speed_of_sound),
        s0(src.s0),
        s1(src.s1),
        error_message_buffer(src.error_message_buffer) {
    setErrorMessageBuffer(*this, src);
  }  // end of BehaviourData

  BehaviourData& BehaviourData::operator=(BehaviourData&& src) {
    if (&src != this) {
      this->dt = src.dt;
      this->K = std::move(src.K);
      this->rdt = src.rdt;
      this->speed_of_sound = src.speed_of_sound;
      this->s0 = std::move(src.s0);
      this->s1 = std::move(src.s1);
      this->error_message_buffer = src.error_message_buffer;
      // the error message is changed once the copies which may throw have
      // succeeded
      this->error_message = src.error_message;
      setErrorMessageBuffer(*this, src);
    }
    return *this;
  }  // end of operator=

  BehaviourData& BehaviourData::operator=(const BehaviourData& src) {
    if (&src != this) {
      this->dt = src.dt;
      this->K = src.K;
      this->rdt = src.rdt;
      this->speed_of_sound = src.speed_of_sound;
      this->s0 = src.s0;
      this->s1 = src.s1;
      this->error_message_buffer = src.error_message_buffer;
      // the error message is changed once the copies which may throw have
      // succeeded
      this->error_message = src.error_message;
      setErrorMessageBuffer(*this, src);
    }
    return *this;
  }  // end of operator=

  // Note: initialiazing s1 with s0 is licit according to the C++ standard:
  //
  // 12.6.2.5
  // Initialization shall proceed in the following order:
  // ...
  // Then, nonstatic data members shall be initialized in the order they were
  // declared in the class definition (again regardless of the order of the
  // mem-initializers).
  BehaviourData::BehaviourData(const Behaviour& b) : s0(b), s1(s0) {
    this->error_message = this->error_message_buffer.data();
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
