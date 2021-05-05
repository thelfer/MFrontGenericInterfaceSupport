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

    std::string getVariableTypeAsString(const Variable &v) {
      if (v.type == Variable::SCALAR) {
        return "Scalar";
      } else if (v.type == Variable::VECTOR) {
        return "Vector";
      } else if (v.type == Variable::STENSOR) {
        return "Stensor";
      }
      if (v.type != Variable::TENSOR) {
        mgis::raise("getVariableTypeAsString: unsupported variable type");
      }
      return "Tensor";
    }

    size_type getVariableSize(const Variable &v, const Hypothesis h) {
      if (v.type == Variable::SCALAR) {
        return 1;
      } else if (v.type == Variable::VECTOR) {
        return getSpaceDimension(h);
      } else if (v.type == Variable::STENSOR) {
        return getStensorSize(h);
      }
      if (v.type != Variable::TENSOR) {
        mgis::raise("getArraySize: unsupported variable type");
      }
      return getTensorSize(h);
    }  // end of getVariableSize

    bool contains(const std::vector<Variable> &vs, const string_view n) {
      return std::find_if(vs.begin(), vs.end(), [&n](const Variable &v) {
               return v.name == n;
             }) != vs.end();
    }  // end of contains

    const Variable &getVariable(const std::vector<Variable> &vs,
                                const string_view n) {
      const auto p =
          std::find_if(vs.begin(), vs.end(),
                       [&n](const Variable &v) { return v.name == n; });
      if (p == vs.end()) {
        mgis::raise("getVariable: no variable named '" + std::string(n) + "'");
      }
      return *p;
    }  // end of getVariable

    size_type getArraySize(const std::vector<Variable> &vs,
                           const Hypothesis h) {
      auto s = size_type{};
      for (const auto &v : vs) {
        s += getVariableSize(v, h);
      }
      return s;
    }  // end of getArraySize

    size_type getVariableOffset(const std::vector<Variable> &vs,
                                const string_view n,
                                const Hypothesis h) {
      auto o = size_type{};
      for (const auto &v : vs) {
        if (v.name == n) {
          return o;
        }
        o += getVariableSize(v, h);
      }
      raise("getVariableOffset: no variable named '" + std::string(n) + "'");
    }  // end of getVariableOffset

  }  // end of namespace behaviour

}  // end of namespace mgis
