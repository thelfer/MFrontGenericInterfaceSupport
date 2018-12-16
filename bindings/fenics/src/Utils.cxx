/*!
 * \file   bindings/fencis/src/Utils.hxx
 * \brief    
 * \author Thomas Helfer
 * \date   14/12/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <dolfin/fem/FiniteElement.h>
#include "MGIS/FEniCS/Utils.hxx"


namespace mgis{

  namespace fenics {

    mgis::size_type getNumberOfIntegrationPoints(const dolfin::FiniteElement& e){
      return e.value_dimension(0);
    } // end of getNumberOfIntegrationPoints

    
    
  }  // end of namespace fenics

}  // end of namespace mgis
