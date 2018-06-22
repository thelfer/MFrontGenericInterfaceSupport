/*!
 * \file   Behaviour.hxx
 * \brief    
 * \author Thomas Helfer
 * \date   19/06/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MFRONT_BEHAVIOUR_BEHAVIOUR_HXX
#define LIB_MFRONT_BEHAVIOUR_BEHAVIOUR_HXX

#include "MFront/Config.hxx"

namespace mfront {

  namespace behaviour {

    using Behaviour = void (*)(const real);

  }  // end of namespace behaviour

} // end of namespace mfront

#endif /* LIB_MFRONT_BEHAVIOUR_BEHAVIOUR_HXX */
