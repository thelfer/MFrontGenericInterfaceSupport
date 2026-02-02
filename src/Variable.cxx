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

#include <bitset>
#include <climits>
#include <algorithm>
#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/Variable.hxx"

namespace mgis::behaviour::internals {

  static int extractAndShift(int &v, const std::size_t s) {
    using bits = std::bitset<sizeof(int) * CHAR_BIT>;
    const auto m = (bits{}.set() << s).flip();
    bits value(v);
    const auto r = value & m;
    value >>= s;
    v = static_cast<int>(value.to_ulong());
    return static_cast<int>(r.to_ulong());
  }  // end of extractAndShift

  static mgis::behaviour::Variable::Type getTinyVectorVariableType(int &v) {
    // tiny vector
    const auto N = extractAndShift(v, 2);
    if (N == 0) {
      return mgis::behaviour::Variable::VECTOR;
    } else if (N == 1) {
      return mgis::behaviour::Variable::VECTOR_1D;
    } else if (N == 2) {
      return mgis::behaviour::Variable::VECTOR_2D;
    } else if (N == 3) {
      return mgis::behaviour::Variable::VECTOR_3D;
    } else {
      mgis::raise("getTinyVectorVariableType: invalid space dimension");
    }
  }  // end of getTinyVectorVariableType

  static mgis::behaviour::Variable::Type getSymmetricTensorVariableType(
      int &v) {
    const auto N = extractAndShift(v, 2);
    if (N == 0) {
      return mgis::behaviour::Variable::STENSOR;
    } else if (N == 1) {
      return mgis::behaviour::Variable::STENSOR_1D;
    } else if (N == 2) {
      return mgis::behaviour::Variable::STENSOR_2D;
    } else if (N == 3) {
      return mgis::behaviour::Variable::STENSOR_3D;
    } else {
      mgis::raise("getSymmetricTensorVariableType: invalid space dimension");
    }
  }  // end of getSymmetricTensorVariableType

  static mgis::behaviour::Variable::Type getUnSymmetricTensorVariableType(
      int &v) {
    const auto N = extractAndShift(v, 2);
    if (N == 0) {
      return mgis::behaviour::Variable::TENSOR;
    } else if (N == 1) {
      return mgis::behaviour::Variable::TENSOR_1D;
    } else if (N == 2) {
      return mgis::behaviour::Variable::TENSOR_2D;
    } else if (N == 3) {
      return mgis::behaviour::Variable::TENSOR_3D;
    } else {
      mgis::raise("getUnSymmetricTensorVariableType: invalid space dimension");
    }
  }  // end of getUnSymmetricTensorVariableType

  static mgis::behaviour::Variable::Type getVariableType(const int t) {
    int v = t;
    const auto type = extractAndShift(v, 3);
    if (type == 0) {
      return Variable::SCALAR;
    } else if (type == 1) {
      return getSymmetricTensorVariableType(v);
    } else if (type == 2) {
      return getTinyVectorVariableType(v);
    } else if (type == 3) {
      return getUnSymmetricTensorVariableType(v);
    } else if (type == 4) {
      return Variable::HIGHER_ORDER_TENSOR;
    } else if (type != 5) {
      mgis::raise("getVariableType: unsupported variable type");
    }
    return Variable::ARRAY;
  }  // end of getVariableType

  static size_t getVariableSize(int &, const mgis::behaviour::Hypothesis);

  static size_t getTinyVectorVariableSize(int &v,
                                          const mgis::behaviour::Hypothesis h) {
    // tiny vector
    const auto N = extractAndShift(v, 2);
    if (N < 0) {
      mgis::raise("invalid tensorial object dimension");
    }
    if (N == 0) {
      return mgis::behaviour::getSpaceDimension(h);
    } else if (N > 3) {
      mgis::raise("invalid space dimension");
    }
    return static_cast<size_t>(N);
  }  // end of getTinyVectorVariableSize

  static size_t getSymmetricTensorVariableSize(
      int &v, const mgis::behaviour::Hypothesis h) {
    const auto N = extractAndShift(v, 2);
    if (N == 0) {
      return mgis::behaviour::getStensorSize(h);
    } else if (N == 1) {
      return 3u;
    } else if (N == 2) {
      return 4u;
    } else if (N == 3) {
      return 6u;
    } else {
      mgis::raise("invalid space dimension");
    }
  }  // end of getSymmetricTensorVariableSize

  static size_t getUnSymmetricTensorVariableSize(
      int &v, const mgis::behaviour::Hypothesis h) {
    const auto N = extractAndShift(v, 2);
    if (N == 0) {
      return mgis::behaviour::getTensorSize(h);
    } else if (N == 1) {
      return 3u;
    } else if (N == 2) {
      return 5u;
    } else if (N == 3) {
      return 9u;
    } else {
      mgis::raise("invalid space dimension");
    }
  }  // end of getUnSymmetricTensorVariableSize

  static size_t getArrayVariableSize(int &v,
                                     const mgis::behaviour::Hypothesis h) {
    // number of dimension of the array
    const auto a = extractAndShift(v, 3);
    if (a == 0) {
      mgis::raise("invalid array arity");
    }
    auto n = size_t{1};
    for (int i = 0; i != a; ++i) {
      const auto d = extractAndShift(v, 7);
      if (d < 1) {
        mgis::raise("invalid array dimension");
      }
      n *= d;
    }
    return n * getVariableSize(v, h);
  }  // end of getArrayVariableSize

  static size_t getVariableSize(int &v, const mgis::behaviour::Hypothesis h) {
    const auto type = extractAndShift(v, 3);
    if (type == 0) {
      return 1u;
    } else if (type == 1) {
      return getSymmetricTensorVariableSize(v, h);
    } else if (type == 2) {
      return getTinyVectorVariableSize(v, h);
    } else if (type == 3) {
      return getUnSymmetricTensorVariableSize(v, h);
    } else if (type == 4) {
      // derivative type
      const auto s1 = getVariableSize(v, h);
      const auto s2 = getVariableSize(v, h);
      return s1 * s2;
    } else if (type != 5) {
      mgis::raise("unsupported variable type");
    }
    return getArrayVariableSize(v, h);
  }

  std::string getVariableTypeSymbolicRepresentation(int &);

  static std::string getTinyVectorVariableTypeSymbolicRepresentation(int &v) {
    // tiny vector
    const auto N = extractAndShift(v, 2);
    if (N == 0) {
      return "tvector<N, real>";
    }
    return "tvector<" + std::to_string(N) + ", real>";
  }  // end of getVariableTypeSymbolicRepresentation

  static std::string getSymmetricTensorVariableTypeSymbolicRepresentation(
      int &v) {
    const auto N = extractAndShift(v, 2);
    if (N == 0) {
      return "stensor<N, real>";
    }
    return "stensor<" + std::to_string(N) + ", real>";
  }  // end of getSymmetricTensorVariableTypeSymbolicRepresentation

  static std::string getUnSymmetricTensorVariableTypeSymbolicRepresentation(
      int &v) {
    const auto N = extractAndShift(v, 2);
    if (N == 0) {
      return "tensor<N, real>";
    }
    return "tensor<" + std::to_string(N) + ", real>";
  }  // end of getUnSymmetricTensorVariableTypeSymbolicRepresentation

  static std::string getArrayVariableTypeSymbolicRepresentation(int &v) {
    // number of dimension of the array
    const auto a = extractAndShift(v, 3);
    if (a == 0) {
      mgis::raise("invalid array arity");
    }
    auto r = std::string{"array<"};
    for (int i = 0; i != a; ++i) {
      const auto d = extractAndShift(v, 7);
      r += std::to_string(d) + ", ";
    }
    return r + getVariableTypeSymbolicRepresentation(v) + ">";
  }  // end of getArrayVariableTypeSymbolicRepresentation

  std::string getVariableTypeSymbolicRepresentation(int &v) {
    const auto type = extractAndShift(v, 3);
    if (type == 0) {
      return "real";
    } else if (type == 1) {
      return getSymmetricTensorVariableTypeSymbolicRepresentation(v);
    } else if (type == 2) {
      return getTinyVectorVariableTypeSymbolicRepresentation(v);
    } else if (type == 3) {
      return getUnSymmetricTensorVariableTypeSymbolicRepresentation(v);
    } else if (type == 4) {
      // derivative type
      const auto s1 = getVariableTypeSymbolicRepresentation(v);
      const auto s2 = getVariableTypeSymbolicRepresentation(v);
      return "derivative_type<" + s1 + ", " + s2 + ">";
    } else if (type != 5) {
      mgis::raise("unsupported variable type");
    }
    return getArrayVariableTypeSymbolicRepresentation(v);
  }  // end of getVariableTypeSymbolicRepresentation

}  // end of namespace mgis::behaviour::internals

namespace mgis::behaviour {

  std::string getVariableTypeAsString(const Variable &v) {
    if (v.type == Variable::SCALAR) {
      return "Scalar";
    } else if (v.type == Variable::VECTOR) {
      return "Vector";
    } else if (v.type == Variable::VECTOR_1D) {
      return "Vector_1D";
    } else if (v.type == Variable::VECTOR_2D) {
      return "Vector_2D";
    } else if (v.type == Variable::VECTOR_3D) {
      return "Vector_3D";
    } else if (v.type == Variable::STENSOR) {
      return "Stensor";
    } else if (v.type == Variable::STENSOR_1D) {
      return "Stensor_1D";
    } else if (v.type == Variable::STENSOR_2D) {
      return "Stensor_2D";
    } else if (v.type == Variable::STENSOR_3D) {
      return "Stensor_3D";
    } else if (v.type == Variable::TENSOR) {
      return "Tensor";
    } else if (v.type == Variable::TENSOR_1D) {
      return "Tensor_1D";
    } else if (v.type == Variable::TENSOR_2D) {
      return "Tensor_2D";
    } else if (v.type == Variable::TENSOR_3D) {
      return "Tensor_3D";
    } else if (v.type == Variable::HIGHER_ORDER_TENSOR) {
      return "HigherOrderTensor";
    } else if (v.type == Variable::ARRAY) {
      return "Array";
    }
    mgis::raise("getVariableTypeAsString: unsupported variable type");
  }

  Variable::Type getVariableType(const int id) {
    return internals::getVariableType(id);
  }

  size_type getVariableSize(const Variable &v, const Hypothesis h) {
    auto id = v.type_identifier;
    const auto s = internals::getVariableSize(id, h);
    if (id != 0) {
      mgis::raise("getVariableSize: invalid type identifier '" +
                  std::to_string(id) + "'");
    }
    return s;
  }  // end of getVariableSize

  std::optional<size_type> getVariableSize(Context &ctx,
                                           const Variable &v,
                                           const Hypothesis h) noexcept {
    try {
      return getVariableSize(v, h);
    } catch (...) {
      std::ignore = registerExceptionInErrorBacktrace(ctx);
    }
    return {};
  }  // end of getVariableSize

  bool contains(const std::vector<Variable> &vs,
                const std::string_view n) noexcept {
    return std::find_if(vs.begin(), vs.end(), [&n](const Variable &v) {
             return v.name == n;
           }) != vs.end();
  }  // end of contains

  const Variable &getVariable(const std::vector<Variable> &vs,
                              const std::string_view n) {
    const auto p = std::find_if(
        vs.begin(), vs.end(), [&n](const Variable &v) { return v.name == n; });
    if (p == vs.end()) {
      mgis::raise("getVariable: no variable named '" + std::string(n) + "'");
    }
    return *p;
  }  // end of getVariable

  std::optional<const Variable *> getVariable(
      Context &ctx,
      const std::vector<Variable> &vs,
      const std::string_view n) noexcept {
    const auto p = std::find_if(
        vs.begin(), vs.end(), [&n](const Variable &v) { return v.name == n; });
    if (p == vs.end()) {
      return ctx.registerErrorMessage("getVariable: no variable named '" +
                                      std::string(n) + "'");
    }
    return &(*p);
  }  // end of getVariable

  size_type getArraySize(const std::vector<Variable> &vs, const Hypothesis h) {
    auto s = size_type{};
    for (const auto &v : vs) {
      s += getVariableSize(v, h);
    }
    return s;
  }  // end of getArraySize

  size_type getVariableOffset(const std::vector<Variable> &vs,
                              const std::string_view n,
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

  std::string getVariableTypeSymbolicRepresentation(const int id) {
    auto t = id;
    const auto s = internals::getVariableTypeSymbolicRepresentation(t);
    if (t != 0) {
      mgis::raise(
          "getVariableTypeSymbolicRepresentation: "
          "invalid type identifier");
    }
    return s;
  }  // end of getVariableTypeSymbolicRepresentation

}  // end of namespace mgis::behaviour
