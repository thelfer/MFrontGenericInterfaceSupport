/*!
 * \file   BehaviourFctPtr.hxx
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

#ifndef LIB_MGIS_BEHAVIOUR_BEHAVIOURFCTPTR_HXX
#define LIB_MGIS_BEHAVIOUR_BEHAVIOURFCTPTR_HXX

/* forward declaration */
typedef struct mgis_bv_BehaviourDataView mgis_bv_BehaviourDataView;
/*! a simple alias */
typedef int (*mgis_bv_BehaviourFctPtr)(mgis_bv_BehaviourDataView* const);

#ifdef __cplusplus

namespace mgis {

  namespace behaviour {

    /*! \brief a simple alias */
    using BehaviourFctPtr = mgis_bv_BehaviourFctPtr;

  }  // end of namespace behaviour

} // end of namespace mgis

#endif /* __cplusplus */

#endif /* LIB_MGIS_BEHAVIOUR_BEHAVIOURFCTPTR_HXX */
