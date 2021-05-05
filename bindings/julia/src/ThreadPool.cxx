/*!
 * \file   bindings/julia/src/ThreadPool.cxx
 * \brief
 * \author Thomas Helfer
 * \date   14/05/2019
 */

#include <jlcxx/jlcxx.hpp>
#include "MGIS/ThreadPool.hxx"

void declareThreadPool(jlcxx::Module& m) {
  m.add_type<mgis::ThreadPool>("ThreadPool").constructor<mgis::size_type>();
}  // end of declareThreadPool
