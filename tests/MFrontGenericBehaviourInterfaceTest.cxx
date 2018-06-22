/*!
 * \file   MFrontGenericBehaviourInterfaceTest.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   20/06/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <cstdlib>
#include "MFront/Behaviour/Description.hxx"

int main(const int argc, const char* const* argv) {
  constexpr const auto h = mfront::behaviour::Hypothesis::TRIDIMENSIONAL;
  if (argc != 3) {
    return EXIT_FAILURE;
  }
  const auto d = mfront::behaviour::load(argv[1], argv[2], h);
  return EXIT_SUCCESS;
}
