/*!
 * \file   MaterialDataManager.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   21/08/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"

namespace mgis {

  namespace behaviour {

    MaterialDataManager::MaterialDataManager(MaterialDataManager&&) = default;
    MaterialDataManager::MaterialDataManager(const MaterialDataManager&) =
        default;

    MaterialDataManager& MaterialDataManager::operator=(
        MaterialDataManager&& src) {
      mgis::raise_if(&src.b != &this->b,
                     "MaterialDataManager::operator=: unmatched behaviour");
      mgis::raise_if(src.K_stride != this->K_stride,
                     "MaterialDataManager::operator=: "
                     "unmatched stiffness matrix size");
      if (&src != this) {
        this->s0 = std::move(src.s0);
        this->s1 = std::move(src.s1);
        this->K = std::move(src.K);
      }
      return *this;
    }  // end of MaterialDataManager::operator=

    MaterialDataManager& MaterialDataManager::operator=(
        const MaterialDataManager& src) {
      mgis::raise_if(&src.b != &this->b,
                     "MaterialDataManager::operator=: unmatched behaviour");
      mgis::raise_if(src.K_stride != this->K_stride,
                     "MaterialDataManager::operator=: "
                     "unmatched stiffness matrix size");
      if (&src != this) {
        this->s0 = src.s0;
        this->s1 = src.s1;
        this->K  = src.K;
      }
      return *this;
    }  // end of MaterialDataManager::operator=

    // Note: initialiazing s1 with s0 is licit according to the C++ standard:
    //
    // 12.6.2.5
    // Initialization shall proceed in the following order:
    // ...
    // Then, nonstatic data members shall be initialized in the order they were
    // declared in the class definition (again regardless of the order of the
    // mem-initializers).
    MaterialDataManager::MaterialDataManager(const Behaviour& behaviour,
                                             const size_type s)
        : s0(behaviour, s),
          s1(s0),
          K_stride(s0.gradients.size() * s0.thermodynamic_forces.size()),
          b(behaviour) {
      constexpr const auto zero = real{0};
      this->K.resize(s * this->K_stride, zero);
    }  // end of Behaviour::Behaviour

    MaterialDataManager::~MaterialDataManager() = default;

    void update(MaterialDataManager& d) {
      std::fill(d.K.begin(), d.K.end(), real{0});
      d.s0 = d.s1;
    }  // end of update

    void revert(MaterialDataManager& d) {
      std::fill(d.K.begin(), d.K.end(), real{0});
      d.s1 = d.s0;
    }  // end of update

  }  // end of namespace behaviour

}  // end of namespace mgis
