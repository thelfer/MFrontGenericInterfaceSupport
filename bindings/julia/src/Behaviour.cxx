/*!
 * \file   bindings/julia/src/Behaviour.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   14/05/2019
 */

#include <jlcxx/jlcxx.hpp>
#include <jlcxx/const_array.hpp>
#include "MGIS/Behaviour/Behaviour.hxx"

namespace jlcxx {
  template <>
  struct IsBits<mgis::behaviour::Hypothesis> : std::true_type {};
} // end of namespace jlcxx

void declareBehaviour(jlcxx::Module& m) {
  using mgis::behaviour::Behaviour;
  using mgis::behaviour::Hypothesis;

  Behaviour (*load)(const std::string&, const std::string&,
                    const Hypothesis) = &mgis::behaviour::load;
  m.add_type<Behaviour>("Behaviour")
      .method("get_material_properties",
              [](const Behaviour& b) { return b.mps; });
  m.method("load", load);
} // end of declareBehaviour
