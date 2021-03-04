/*!
 * \file   include/MGIS/MatrixView.hxx
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

#ifndef LIB_MGIS_MATRIXVIEW_HXX
#define LIB_MGIS_MATRIXVIEW_HXX

#include "MGIS/Config.hxx"

namespace mgis {

  /*!
   * \brief a structure in charge of creating a sub matrix view from a parent
   * matrix stored in a continous memory of data using a row major indexing
   * scheme. The number of of colums of the parent matrix is denoted the
   * row_stride here.
   */
  struct MGIS_EXPORT MatrixView {
    /*!
     * \brief constructor
     * \param[in] v: values
     * \param[in] o: offset of the first value
     * \param[in] s: row stride
     */
    MatrixView(real* const, const size_type, const size_type);
    //! \brief move constructor
    MatrixView(MatrixView&&);
    //! \brief copy constructor
    MatrixView(const MatrixView&);
    /*!
     * \brief access operator
     * \param[in] i: row index
     * \param[in] j: column index
     */
    inline real& operator()(const size_type, const size_type);
    /*!
     * \brief access operator
     * \param[in] i: row index
     * \param[in] j: column index
     */
    inline const real& operator()(const size_type, const size_type) const;

   private:
    //! underlying values
    real* const values;
    /*!
     * \brief offset of the first element of the sub matrix in the array holding
     * the parent matrix.
     */
    const size_type ioffset;
    //! number of columns of the matrix
    const size_type row_stride;
  };  // end of struct MatrixView

}  // end of mgis

#include "MGIS/MatrixView.ixx"

#endif /* LIB_MGIS_MATRIXVIEW_HXX */
