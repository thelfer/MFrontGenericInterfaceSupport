/*!
 * \file   MGIS/QuadratureFunction/AbstractQuadratureSpace.hxx
 * \brief    
 * \author Thomas Helfer
 * \date   01/05/2025
 */

#ifndef LIB_MGIS_QUADRATUREFUNCTION_ABSTRACTQUADRATURESPACE_HXX
#define LIB_MGIS_QUADRATUREFUNCTION_ABSTRACTQUADRATURESPACE_HXX

#include "MGIS/Config.hxx"

namespace mgis::quadrature_function {


  struct MGIS_EXPORT AbstractQuadratureSpace {
    virtual size_type getNumberOfIntegrationPoints() const = 0;
    //! \brief destructor
    virtual ~AbstractQuadratureSpace();
  };

} // end of namespace mgis::quadrature_function

#endif /* LIB_MGIS_QUADRATUREFUNCTION_ABSTRACTQUADRATURESPACE_HXX */
