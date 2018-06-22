/*!
 * \file   Variable.cxx
 * \brief
 * \author Thomas Helfer
 * \date   19/06/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include "MFront/Raise.hxx"
#include "MFront/Behaviour/Variable.hxx"

namespace mfront {

  namespace behaviour {

    static size_type getStensorSize(const Hypothesis h) {
      if ((h == Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN) ||
          (h == Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS)) {
        return 3u;
      } else if ((h == Hypothesis::AXISYMMETRICAL) ||
                 (h == Hypothesis::PLANESTRESS) ||
                 (h == Hypothesis::PLANESTRAIN) ||
                 (h == Hypothesis::GENERALISEDPLANESTRAIN)) {
        return 4u;
      } else if (h == Hypothesis::TRIDIMENSIONAL) {
        return 6u;
      }
      mfront::raise("getStensorSize: unsupported modelling hypothesis");
    }  // end of getStensorSize

    size_type getTensorSize(const Hypothesis h) {
      if ((h == Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN) ||
          (h == Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS)) {
        return 3u;
      } else if ((h == Hypothesis::AXISYMMETRICAL) ||
                 (h == Hypothesis::PLANESTRESS) ||
                 (h == Hypothesis::PLANESTRAIN) ||
                 (h == Hypothesis::GENERALISEDPLANESTRAIN)) {
        return 5u;
      } else if (h == Hypothesis::TRIDIMENSIONAL) {
        return 9u;
      }
      mfront::raise("getTensorSize: unsupported modelling hypothesis");
    }  // end of getTensorSize

    size_type getArraySize(const std::vector<Variable>& vs,
                           const Hypothesis h) {
      auto getSize = [h](const Variable &v) -> size_type {
        if (v.type == Variable::SCALAR) {
          return 1;
        } else if (v.type == Variable::STENSOR) {
          return getStensorSize(h);
        }
        if (v.type != Variable::TENSOR) {
          mfront::raise("getArraySize: unsupported variable type");
        }
        return getTensorSize(h);
      };
      auto s = size_type{};
      for (const auto& v : vs) {
        s += getSize(v);
      }
      return s;
    }  // end of getArraySize

  }  // end of namespace behaviour

}  // end of namespace mfront
