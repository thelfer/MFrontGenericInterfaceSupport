/*!
 * \file   bindings/julia/src/mgis-module.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   14/05/2019
 */

#include <jlcxx/jlcxx.hpp>

// forward declarations
void declareThreadPool(jlcxx::Module& m);

JLCXX_MODULE define_mgis_module(jlcxx::Module& m){
  declareThreadPool(m);
}
