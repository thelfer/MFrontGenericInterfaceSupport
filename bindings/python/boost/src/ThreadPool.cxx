/*!
 * \file   ThreadPool.cxx
 * \brief
 * \author Thomas Helfer
 * \date   12/11/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <boost/python/class.hpp>
#include "MGIS/ThreadPool.hxx"

void declareThreadPool() {
  boost::python::class_<mgis::ThreadPool, boost::noncopyable>(
      "ThreadPool", boost::python::init<mgis::size_type>());
}  // end of declareThreadPool
