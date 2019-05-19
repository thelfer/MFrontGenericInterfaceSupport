/*!
 * \file   JuliaUtilities.ixx
 * \brief    
 * \author Thomas Helfer
 * \date   17/05/2019
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */


#ifndef LIB_MGIS_JULIA_ARRAYVIEW_HXX
#define LIB_MGIS_JULIA_ARRAYVIEW_HXX

#include <cstdint>
#include <jlcxx/jlcxx.hpp>
#include <jlcxx/tuple.hpp>

namespace mgis {

  namespace julia {

    typedef std::int64_t index_t;

    namespace detail {
      // Helper to make a C++ tuple of longs based on the number of elements
      template <index_t N, typename... TypesT>
      struct LongNTuple {
        typedef typename LongNTuple<N - 1, index_t, TypesT...>::type type;
      };

      template <typename... TypesT>
      struct LongNTuple<0, TypesT...> {
        typedef std::tuple<TypesT...> type;
      };
    }

    /// Wrap a const pointer
    template <typename T>
    struct Ptr {
      const T* ptr;
    };

    template <typename T>
    struct IsBits<Ptr<T>> : std::true_type {};

    template <typename T>
    struct InstantiateParametricType<Ptr<T>> {
      int operator()(Module&) const {
        // Register the Julia type if not already instantiated
        if (!static_type_mapping<Ptr<T>>::has_julia_type()) {
          jl_datatype_t* dt = (jl_datatype_t*)apply_type(
              (jl_value_t*)julia_type("Ptr"),
              jl_svec1(static_type_mapping<T>::julia_type()));
          set_julia_type<Ptr<T>>(dt);
          protect_from_gc(dt);
        }
        return 0;
      }
    };

    /// Wrap a pointer, providing the Julia array interface for it
    /// The parameter N represents the number of dimensions
    template <typename T, index_t N>
    class ArrayView {
     public:
      typedef typename detail::LongNTuple<N>::type size_t;

      template <typename... SizesT>
      ArrayView(const T* ptr, const SizesT... sizes)
          : m_arr(ptr), m_sizes(sizes...) {}

      T getindex(const std::int64_t i) const { return m_arr[i - 1]; }

      size_t size() const { return m_sizes; }

      T* ptr() const { return m_arr; }

     private:
      T* m_arr;
      const size_t m_sizes;
    };

    template <typename T, index_t N>
    struct InstantiateParametricType<ArrayView<T, N>>
        : InstantiateParametricType<Ptr<T>> {};

    template <typename T, typename... SizesT>
    ArrayView<T, sizeof...(SizesT)> make_array_view(const T* p,
                                                    const SizesT... sizes) {
      return ArrayView<T, sizeof...(SizesT)>(p, sizes...);
    }

    template <typename T, index_t N>
    struct IsImmutable<ArrayView<T, N>> : std::false_type {};

    template <typename T, index_t N>
    struct ConvertToJulia<ArrayView<T, N>> {
      jl_value_t* operator()(const ArrayView<T, N>& arr) {
        jl_value_t* result = nullptr;
        jl_value_t* ptr = nullptr;
        jl_value_t* size = nullptr;
        JL_GC_PUSH3(&result, &ptr, &size);
        ptr = box(Ptr<T>({arr.ptr()}));
        size = convert_to_julia(arr.size());
        result = jl_new_struct(julia_type<ArrayView<T, N>>(), ptr, size);
        JL_GC_POP();
        return result;
      }
    };

    template <typename T, index_t N>
    struct static_type_mapping<ArrayView<T, N>> {
      typedef jl_value_t* type;
      static jl_datatype_t* julia_type() {
        static jl_datatype_t* app_dt = nullptr;
        if (app_dt == nullptr) {
          jl_datatype_t* pdt =
              (jl_datatype_t*)::jlcxx::julia_type("ArrayView");
          jl_value_t* boxed_n = box(N);
          JL_GC_PUSH1(&boxed_n);
          app_dt = (jl_datatype_t*)apply_type(
              (jl_value_t*)pdt, jl_svec2(::jlcxx::julia_type<T>(), boxed_n));
          protect_from_gc(app_dt);
          JL_GC_POP();
        }
        return app_dt;
      }
    };

  }  // namespace julia

}  // namespace mgis

#endif /* LIB_MGIS_JULIA_ARRAYVIEW_HXX */
