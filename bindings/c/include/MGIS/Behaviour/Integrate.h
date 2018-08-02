/*!
 * \file   Integrate.h
 * \brief    
 * \author Thomas Helfer
 * \date   02/08/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_BEHAVIOUR_INTEGRATE_H
#define LIB_MGIS_BEHAVIOUR_INTEGRATE_H

#include "MGIS/Status.h"
#include "MGIS/Behaviour/BehaviourDataView.h"
#include "MGIS/Behaviour/Behaviour.h"

#ifdef __cplusplus
extern "C" {
#endif /*  __cplusplus */

mgis_status integrate(mgis_bv_BehaviourDataView* const,
                      const mgis_bv_Behaviour* const);

#ifdef __cplusplus
}
#endif /*  __cplusplus */

#endif /*LIB_MGIS_BEHAVIOUR_INTEGRATE_H */
