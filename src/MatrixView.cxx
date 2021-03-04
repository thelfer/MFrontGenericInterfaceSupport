/*!
 * \file   MatrixView.cxx
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

#include "MGIS/MatrixView.hxx"

namespace mgis {

  MatrixView::MatrixView(real* const v, const size_type o, const size_type s)
      : values(v),
        ioffset(o),
        row_stride(s) {}  // end of MatrixView::MatrixView

  MatrixView::MatrixView(MatrixView&&) = default;

  MatrixView::MatrixView(const MatrixView&) = default;

} // end of mgis
