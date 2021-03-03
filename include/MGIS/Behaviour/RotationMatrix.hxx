/*!
 * \file   include/MGIS/Behaviour/RotationMatrix.hxx
 * \brief
 * \author Thomas Helfer
 * \date   15/02/2021
 */

#ifndef LIB_MGIS_BEHAVIOUR_ROTATIONMATRIX_HXX
#define LIB_MGIS_BEHAVIOUR_ROTATIONMATRIX_HXX

#include <array>
#include <vector>
#include "MGIS/Span.hxx"
#include "MGIS/StorageMode.hxx"

namespace mgis::behaviour {

  /*!
   * \brief an helper structure to store a material axis.
   *
   * In 2D, a material axis is represented either by:
   *
   * - an array of 2 values
   * - an array of 2 x n values where n is the number of integration points.
   *
   * In 3D, a material axis is represented either by:
   *
   * - an array of 3 values
   * - an array of 3 x n values where n is the number of integration points.
   */
  struct MaterialAxisStorage {
    //! \brief constructor from an external array
    MaterialAxisStorage(const mgis::span<const mgis::real> &,
                        const mgis::StorageMode &);
    /*!
     * \brief constructor from a temporary
     * \param[in] values: values to be stored
     */
    MaterialAxisStorage(std::vector<mgis::real> &&);
    //! \brief destructor
    ~MaterialAxisStorage();

   private:
    //! \brief internal storage, if required
    std::vector<mgis::real> a_values;

   public:
    //! \brief description of the material axis
    const mgis::span<const mgis::real> a;
  };  // end of struct MaterialAxisStorage

  /*!
   * \brief a structure to represent a 2D rotation matrix
   * \note a rotation matrix can be represented by:
   * - an array of 2 values
   * - an array of 2 x n values where n is the number of integration points
   */
  struct MGIS_EXPORT RotationMatrix2D : public MaterialAxisStorage {
    /*!
     * \brief constructor from an external array
     * \param[in] v: values defining the first material axis
     * \param[in] s: storage mode used for the first material axis
     */
    RotationMatrix2D(const mgis::span<const mgis::real> &,
                     const mgis::StorageMode &);
    //! \brief destructor
    ~RotationMatrix2D();
  };  // end of struct RotationMatrix2D

  //! \brief a structure to represent a 3D rotation matrix
  struct MGIS_EXPORT RotationMatrix3D {
    /*!
     * \brief default constructor
     * \param[in] v1: first material axis
     * \param[in] s1: storage mode used for the first material axis
     * \param[in] v2: second material axis
     * \param[in] s1: storage mode used for the second material axis
     */
    RotationMatrix3D(const mgis::span<const mgis::real> &,
                     const mgis::StorageMode &,
                     const mgis::span<const mgis::real> &,
                     const mgis::StorageMode &);
    /*!
     * \brief default constructor
     * \param[in] v1: first material axis
     * \param[in] v2: second material axis
     * \param[in] s: storage mode used for the second material axes
     */
    RotationMatrix3D(const mgis::span<const mgis::real> &,
                     const mgis::span<const mgis::real> &,
                     const mgis::StorageMode &);
    //! \brief first material axis
    const MaterialAxisStorage a1;
    //! \brief second material axis
    const MaterialAxisStorage a2;
    //! \brief destructor
    ~RotationMatrix3D();
  };  // end of struct RotationMatrix2D

  /*!
   * \brief an helper function to build a 2D rotation matrix for an in-plane
   * unit vector.
   * \param[in] a: material axis
   */
  std::array<mgis::real, 9u> buildRotationMatrix(
      const mgis::span<const mgis::real, 2u> &);
  /*!
   * \brief an helper function to build a 3D rotation matrix for an in-plane
   * unit vector.
   * \param[in] a1: first material axis
   * \param[in] a2: second material axis
   */
  std::array<mgis::real, 9u> buildRotationMatrix(
      const mgis::span<const mgis::real, 3u> &,
      const mgis::span<const mgis::real, 3u> &);

}  // end of namespace mgis::behaviour

#include "MGIS/Behaviour/RotationMatrix.ixx"

#endif /* LIB_MGIS_BEHAVIOUR_ROTATIONMATRIX_HXX */
