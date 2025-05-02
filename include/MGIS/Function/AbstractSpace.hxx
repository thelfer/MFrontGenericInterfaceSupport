/*!
 * \file   MGIS/Function/AbstractSpace.hxx
 * \brief
 * \author Thomas Helfer
 * \date   01/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_ABSTRACTSPACE_HXX
#define LIB_MGIS_FUNCTION_ABSTRACTSPACE_HXX

#include "MGIS/Config.hxx"

namespace mgis::function {

  /*!
   * \brief a abstract class describing the discretization space on
   * which is built a function.
   */
  struct MGIS_EXPORT AbstractSpace {
    /*!
     * \brief return the size of the space
     *
     * This can be the number of nodes (nodal space),
     * the number of integrations points (quadrature space), etc.
     */
    virtual size_type getSpaceSize() const = 0;
    //! \brief destructor
    virtual ~AbstractSpace() noexcept;
  };

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_ABSTRACTSPACE_HXX */
