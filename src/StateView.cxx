/*!
 * \file   StateView.cxx
 * \brief
 * \author th202608
 * \date   21/08/2025
 */

#include <vector>
#include <ostream>
#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/Variable.hxx"
#include "MGIS/Utilities/Markdown.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/StateView.hxx"

namespace mgis::behaviour {

  static void print_variables(std::ostream& os,
                              const Behaviour& b,
                              const std::vector<Variable>& variables,
                              const mgis::real* const values) {
    auto o = mgis::size_type{};
    for (const auto& v : variables) {
      os << "- " << v.name << " (" << getVariableTypeAsString(v) << "): ";
      if (v.type == Variable::SCALAR) {
        os << values[o] << '\n';
        ++o;
      } else {
        const auto s = getVariableSize(v, b.hypothesis);
        os << '{';
        for (auto i = o; i != o + s;) {
          os << values[i];
          if (++i != o + s) {
            os << ", ";
          }
        }
        os << "}\n";
        o += s;
      }
    }
  }  // end of print_variables

  static void print_markdown_impl(std::ostream& os,
                                  const Behaviour& b,
                                  const auto& s,
                                  const mgis::size_type l) {
    if (!b.gradients.empty()) {
      os << mgis::utilities::get_heading_signs(l + 1) << " Gradients\n\n";
      print_variables(os, b, b.gradients, s.gradients);
      os << '\n';
    }
    if (!b.thermodynamic_forces.empty()) {
      os << mgis::utilities::get_heading_signs(l + 1)
         << " Thermodynamic forces\n\n";
      print_variables(os, b, b.thermodynamic_forces, s.thermodynamic_forces);
      os << '\n';
    }
    if (!b.mps.empty()) {
      os << mgis::utilities::get_heading_signs(l + 1)
         << " Material properties\n\n";
      print_variables(os, b, b.mps, s.material_properties);
      os << '\n';
    }
    if (!b.isvs.empty()) {
      os << mgis::utilities::get_heading_signs(l + 1)
         << " Internal state variables\n\n";
      print_variables(os, b, b.isvs, s.internal_state_variables);
      os << '\n';
    }
    if (!b.esvs.empty()) {
      os << mgis::utilities::get_heading_signs(l + 1)
         << " External state variables\n\n";
      os << "(" << b.esvs.size() << ")\n";
      print_variables(os, b, b.esvs, s.external_state_variables);
      os << '\n';
    }
  }  // end of print_markdown

  void print_markdown(std::ostream& os,
                      const Behaviour& b,
                      const StateView& s,
                      const mgis::size_type l) {
    print_markdown_impl(os, b, s, l);
  }

  void print_markdown(std::ostream& os,
                      const Behaviour& b,
                      const InitialStateView& s,
                      const mgis::size_type l) {
    print_markdown_impl(os, b, s, l);
  }
}  // end of namespace mgis::behaviour
