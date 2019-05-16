/*!
 * \file   bindings/julia/src/behaviour-module.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   14/05/2019
 */

#include <jlcxx/jlcxx.hpp>

// forward declarations
void declareHypothesis(jlcxx::Module&);
void declareVariable(jlcxx::Module&);
void declareBehaviour(jlcxx::Module&);

JLCXX_MODULE define_mgis_behaviour_module(jlcxx::Module& m) {
  declareHypothesis(m);
  declareVariable(m);
  declareBehaviour(m);
}
