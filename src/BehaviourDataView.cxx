/*!
 * \file   BehaviourDataView.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   21/08/2025
 */

#include <ostream>
#include "MGIS/Utilities/Markdown.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/BehaviourDataView.hxx"

namespace mgis::behaviour {

  void print_markdown(std::ostream& os,
                      const Behaviour& b,
                      const BehaviourDataView& d,
                      const mgis::size_type l) {
    os << mgis::utilities::get_heading_signs(l + 1)
       << " Behaviour description\n\n";
    print_markdown(os, b, l + 1);
    os << mgis::utilities::get_heading_signs(l + 1)
       << " State at the beginning of the time step\n";
    print_markdown(os, b, d.s0, l + 1);
    os << mgis::utilities::get_heading_signs(l + 1)
       << " State at the end of the time step\n";
    print_markdown(os, b, d.s1, l + 1);
  }  // end of print_markdown

} // end of namespace mgis::behaviour