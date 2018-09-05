/*!
 * \file   MatrixView.ixx
 * \brief    
 * \author Thomas Helfer
 * \date   03/09/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_MATRIXVIEW_IXX
#define LIB_MGIS_MATRIXVIEW_IXX

namespace mgis {

  real& MatrixView::operator()(const size_type i, const size_type j){
    return this->values[this->ioffset+i*this->row_stride+j];
  } // end of MatrixView::operator()

  const real& MatrixView::operator()(const size_type i, const size_type j) const {
    return this->values[this->ioffset+i*this->row_stride+j];
  } // end of MatrixView::operator()

}  // end of mgis

#endif /* LIB_MGIS_MATRIXVIEW_IXX */
