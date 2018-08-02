/*!
 * \file   BehaviourFctPtr.h
 * \brief    
 * \author Thomas Helfer
 * \date   01/08/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_BEHAVIOUR_BEHAVIOURFCTPTR_H
#define LIB_MGIS_BEHAVIOUR_BEHAVIOURFCTPTR_H

#include "MGIS/Behaviour/BehaviourDataView.h"

typedef int(*MGIS_BV_BehaviourFctPtr)(const MGIS_BV_BehaviourDataView*);

#endif /* LIB_MGIS_BEHAVIOUR_BEHAVIOURFCTPTR_H */
