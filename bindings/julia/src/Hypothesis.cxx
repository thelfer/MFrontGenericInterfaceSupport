/*!
 * \file   bindings/julia/src/Hypothesis.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   14/05/2019
 */

#include <jlcxx/jlcxx.hpp>
#include "MGIS/Behaviour/Hypothesis.hxx"

namespace jlcxx {
  template <>
  struct IsBits<mgis::behaviour::Hypothesis> : std::true_type {};
} // end of namespace jlcxx

void declareHypothesis(jlcxx::Module& m) {
  using mgis::behaviour::Hypothesis;

  mgis::behaviour::Hypothesis (*fromString)(const std::string&) =
      &mgis::behaviour::fromString;

  m.add_bits<Hypothesis>("Hypothesis");
  m.set_const("AXISYMMETRICALGENERALISEDPLANESTRAIN",
              Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN);
  m.set_const("AXISYMMETRICALGENERALISEDPLANESTRESS",
              Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS);
  m.set_const("AXISYMMETRICAL",
              Hypothesis::AXISYMMETRICAL);
  m.set_const("PLANESTRESS", Hypothesis::PLANESTRESS);
  m.set_const("PLANESTRAIN", Hypothesis::PLANESTRAIN);
  m.set_const("GENERALISEDPLANESTRAIN", Hypothesis::GENERALISEDPLANESTRAIN);
  m.set_const("TRIDIMENSIONAL", Hypothesis::TRIDIMENSIONAL);

  m.method("getSpaceDimension", &mgis::behaviour::getSpaceDimension);
  m.method("getStensorSize", &mgis::behaviour::getStensorSize);
  m.method("getTensorSize", &mgis::behaviour::getTensorSize);
  m.method("toString", &mgis::behaviour::toString);
  m.method("fromString", fromString);
}  // end of declareHypothesis

