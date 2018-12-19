/*!
 * \file   bindings/fencis/tests/src/FEniCSTestingUtilities.cxx
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
 *
 * \note This file contains material strongly inspired by the
 * fenics-solid-materials by Kristian Oelgaard and Garth N. Wells.
 * See <https://bitbucket.org/fenics-apps/fenics-solid-mechanics> for
 * details.
 */

#include <dolfin/common/Array.h>
#include "MGIS/FEniCS/FEniCSTestingUtilities.hxx"

namespace mgis{

  namespace fenics{

    struct sx1 : public dolfin::SubDomain
    {
      bool inside(const dolfin::Array<double>& x, bool) const {
	return std::abs(x[0]) < DOLFIN_EPS;
      }
    }; // end of sx1

    struct sx2 : public dolfin::SubDomain
    {
      bool inside(const dolfin::Array<double>& x, bool) const {
	return std::abs(x[0] - 1.0) < DOLFIN_EPS;
      }
    }; // end of sx2

    struct sy1 : public dolfin::SubDomain
    {
      bool inside(const dolfin::Array<double>& x, bool) const {
	return std::abs(x[1]) < DOLFIN_EPS;
      }
    }; // end of sy1

    struct sy2 : public dolfin::SubDomain
    {
      bool inside(const dolfin::Array<double>& x, bool) const {
	return std::abs(x[1] - 1.0) < DOLFIN_EPS;
      }
    }; // end of sy2

    struct sz1 : public dolfin::SubDomain
    {
      bool inside(const dolfin::Array<double>& x, bool) const {
	return std::abs(x[2]) < DOLFIN_EPS;
      }
    }; // end of sz1

    struct sz2 : public dolfin::SubDomain
    {
      bool inside(const dolfin::Array<double>& x, bool) const {
	return std::abs(x[2] - 1.0) < DOLFIN_EPS;
      }
    }; // end of sz2
    
    std::map<std::string,std::shared_ptr<dolfin::SubDomain>>
    getUnitCubeBoundaries(){
      auto boundaries = std::map<std::string,std::shared_ptr<dolfin::SubDomain>>{};
      boundaries.insert({"sx1",std::make_shared<sx1>()});
      boundaries.insert({"sx2",std::make_shared<sx2>()});
      boundaries.insert({"sy1",std::make_shared<sy1>()});
      boundaries.insert({"sy2",std::make_shared<sy2>()});
      boundaries.insert({"sz1",std::make_shared<sz1>()});
      boundaries.insert({"sz2",std::make_shared<sz2>()});
      return boundaries;
    } // end of getUnitCubeBoundaries

  } // end of namespace fenics

} // end of namespace mgis
