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

#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/Variable.hxx"

namespace mgis {

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

    static size_type getTensorSize(const Hypothesis h) {
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

    size_type getVariableSize(const Variable &v, const Hypothesis h) {
      if (v.type == Variable::SCALAR) {
        return 1;
      } else if (v.type == Variable::STENSOR) {
        return getStensorSize(h);
      }
      if (v.type != Variable::TENSOR) {
        mfront::raise("getArraySize: unsupported variable type");
      }
      return getTensorSize(h);
    }  // end of getVariableSize

    size_type getArraySize(const std::vector<Variable> &vs,
                           const Hypothesis h) {
      auto s = size_type{};
      for (const auto &v : vs) {
        s += getVariableSize(v, h);
      }
      return s;
    }  // end of getArraySize

    size_type getVariableOffset(const std::vector<Variable> &vs,
                                const std::string &n,
                                const Hypothesis h) {
      auto o = size_type{};
      for (const auto &v : vs) {
        if (v.name == n) {
          return o;
        }
        o += getVariableSize(v, h);
      }
      raise("getVariableOffset: no variable named '" + n + "'");
    }  // end of getVariableOffset

  }  // end of namespace behaviour

}  // end of namespace mgis
