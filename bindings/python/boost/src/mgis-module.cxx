/*!
 * \file   bindings/python/src/mgis-module.cxx
 * \brief
 * \author Thomas Helfer
 * \date   31/10/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <boost/python/module.hpp>

// forward declarations
void declareThreadPool();

BOOST_PYTHON_MODULE(_mgis) { declareThreadPool(); }  // end of module behaviour
