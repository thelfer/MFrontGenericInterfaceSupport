/*!
 * \file   include/MGIS/StorageMode.hxx
 * \brief    
 * \author Thomas Helfer
 * \date   17/02/2021
 */

#ifndef LIB_MGIS_STORAGEMODE_HXX
#define LIB_MGIS_STORAGEMODE_HXX

namespace mgis{

  /*!
   * \brief storage option for a non uniform material property or non
   * uniform external state variable.
   */
  enum struct StorageMode {
    LOCAL_STORAGE,    //!< \brief use `std::vector` to store the data
    EXTERNAL_STORAGE  //!< \brief use `mgis::span`  to store the data
  };                  // end of StorageMode

}  // namespace mgis

#endif /* LIB_MGIS_STORAGEMODE_HXX */
