/*!
 * \file   include/MGIS/Behaviour/BehaviourDataView.h
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

#ifndef LIB_MGIS_BEHAVIOUR_BEHAVIOURDATAVIEW_H
#define LIB_MGIS_BEHAVIOUR_BEHAVIOURDATAVIEW_H

#include "MGIS/Config.h"
#include "MGIS/Status.h"
#include "MGIS/Behaviour/BehaviourData.h"
#include "MGIS/Behaviour/BehaviourDataView.hxx"

#ifdef __cplusplus
extern "C" {
#endif /*  __cplusplus */

/*!
 * \brief create a view of a behaviour data
 * \param[in]  v: pointer to the view
 * \param[out] d: pointer to a pointer to the allocated data
 */
MGIS_C_EXPORT mgis_status mgis_bv_makeBehaviourDataView(
    mgis_bv_BehaviourDataView**, mgis_bv_BehaviourData* const);

#ifdef __cplusplus
}
#endif /*  __cplusplus */

#endif /* LIB_MGIS_BEHAVIOUR_BEHAVIOURDATAVIEW_H */
