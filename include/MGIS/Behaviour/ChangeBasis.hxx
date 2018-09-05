/*!
 * \file   include/MGIS/Behaviour/ChangeBasis.hxx
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

#ifndef LIB_MGIS_BEHAVIOUR_CHANGEBASIS_HXX
#define LIB_MGIS_BEHAVIOUR_CHANGEBASIS_HXX

#include "MGIS/Raise.hxx"
#include "MGIS/Config.hxx"
#include "MGIS/Behaviour/Variable.hxx"
#include "MGIS/Behaviour/Hypothesis.hxx"

namespace mgis {

  // forward declaration
  struct MatrixView;

  namespace behaviour {

    /*!
     * \brief a structure in charge of rotating objects in 2D
     */
    struct Rotation2D {
      /*!
       * \brief m: rotation matrix in row-major ordering
       */
      Rotation2D(const real* const);
      /*!
       * \brief 
       */
      Rotation2D(const real, const real, const real, const real);
      //! \brief invert the rotation
      Rotation2D transpose() const;
      /*!
       * \brief rotate a symmetric tensor
       * \param[out] o: rotated values
       * \param[out] i: initial values
       */
      void rotateVector(real* const, const real* const) const;
      /*!
       * \brief rotate a symmetric tensor
       * \param[out] o: rotated values
       * \param[out] i: initial values
       */
      void rotateStensor(real* const, const real* const) const;
      /*!
       * \brief rotate a tensor
       * \param[out] o: rotated values
       * \param[out] i: initial values
       */
      void rotateTensor(real* const, const real* const) const;
      //! \brief
      void buildVectorRotationOperator(mgis::MatrixView&) const;
      //! \brief
      void buildStensorRotationOperator(mgis::MatrixView&) const;
      //! \brief
      void buildTensorRotationOperator(mgis::MatrixView&) const;

     private:
      // coefficients of the rotation matrix
      const real m00, m01, m10, m11;
    };  // end of struct Rotation3D

    /*!
     * \brief a structure in charge of rotating objects in 3D
     */
    struct Rotation3D {
      /*!
       * \brief m: rotation matrix in row-major ordering
       */
      Rotation3D(const real* const);
      /*!
       * \brief 
       */
      Rotation3D(const real,
                 const real,
                 const real,
                 const real,
                 const real,
                 const real,
                 const real,
                 const real,
                 const real);
      //! \brief invert the rotation
      Rotation3D transpose() const;
      /*!
       * \brief rotate a symmetric tensor
       * \param[out] o: rotated values
       * \param[out] i: initial values
       */
      void rotateVector(real* const, const real* const) const;
      /*!
       * \brief rotate a symmetric tensor
       * \param[out] o: rotated values
       * \param[out] i: initial values
       */
      void rotateStensor(real* const, const real* const) const;
      /*!
       * \brief rotate a tensor
       * \param[out] o: rotated values
       * \param[out] i: initial values
       */
      void rotateTensor(real* const, const real* const) const;
      //! \brief
      void buildVectorRotationOperator(mgis::MatrixView&) const;
      //! \brief
      void buildStensorRotationOperator(mgis::MatrixView&) const;
      //! \brief
      void buildTensorRotationOperator(mgis::MatrixView&) const;

     private:
      // coefficients of the rotation matrix
      const real m00, m01, m02, m10, m11, m12, m20, m21, m22;
    };  // end of struct Rotation3D

    /*!
     * \brief change the basis of gradients or thermodynamic fluxes
     * \param[in,out] v: values
     * \param[in] vs: list of variables
     * \param[in] h: hypothesis
     * \param[in] r: rotation matrix, stored in row-major ordering (C ordering)
     */
    MGIS_EXPORT void changeBasis(real* const,
                                 const std::vector<Variable>&,
                                 const Hypothesis,
                                 const real* const);

    /*!
     * \brief change the basis of gradients or thermodynamic fluxes
     * \param[out] o: rotated values
     * \param[in] i: initial values
     * \param[in] vs: list of variables
     * \param[in] h: hypothesis
     * \param[in] r: rotation matrix, stored in row-major ordering (C ordering)
     */
    MGIS_EXPORT void changeBasis(real* const,
                                 const real* const,
                                 const std::vector<Variable>&,
                                 const Hypothesis,
                                 const real* const);
    /*!
     * \brief change the basis of gradients or thermodynamic fluxes
     * \param[out] o: rotated values
     * \param[in] i: initial values
     * \param[in] vs: list of variables
     * \param[in] h: hypothesis
     * \param[in] r: rotation matrix
     */
    template <typename Rotation>
    void changeBasis(real* const,
                     const real* const,
                     const std::vector<Variable>&,
                     const Hypothesis,
                     const Rotation&);

  }  // end of namespace behaviour

}  // end of namespace mgis

#include "MGIS/Behaviour/ChangeBasis.ixx"

#endif /* LIB_MGIS_BEHAVIOUR_CHANGEBASIS_HXX */
