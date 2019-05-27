/*!
 * \file   bindings/julia/src/behaviour-module.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   14/05/2019
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <string>
#include <jlcxx/jlcxx.hpp>
#include "MGIS/Julia/JuliaUtilities.hxx"

// forward declarations
void declareHypothesis(jlcxx::Module&);
void declareVariable(jlcxx::Module&);
void declareBehaviour(jlcxx::Module&);
void declareState(jlcxx::Module&);
void declareBehaviourData(jlcxx::Module&);
void declareBehaviourDataView(jlcxx::Module&);
void declareMaterialStateManager(jlcxx::Module&);
void declareMaterialDataManager(jlcxx::Module&);
void declareIntegrate(jlcxx::Module&);


JLCXX_MODULE define_mgis_behaviour_module(jlcxx::Module& m) {
  mgis::julia::expose_std_vector<std::string>(m, "StringsVector");
  mgis::julia::expose_std_vector<mgis::real>(m, "RealsVector");
  declareHypothesis(m);
  declareVariable(m);
  declareBehaviour(m);
  declareState(m);
  declareBehaviourData(m);
  declareBehaviourDataView(m);
  declareMaterialStateManager(m);
  declareMaterialDataManager(m);
  declareIntegrate(m);
}
