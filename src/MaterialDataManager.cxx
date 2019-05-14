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
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"

namespace mgis {

  namespace behaviour {

    MaterialDataManager::MaterialDataManager(const Behaviour& behaviour,
                                             const size_type s)
        : s0(behaviour, s),
          s1(behaviour, s),
          n(s),
          K_stride(getTangentOperatorArraySize(behaviour)),
          b(behaviour) {
      constexpr const auto zero = real{0};
      this->K_values.resize(this->n * this->K_stride, zero);
      this->K = mgis::span<real>(this->K_values);
    }  // end of Behaviour::Behaviour

    MaterialDataManager::MaterialDataManager(
        const Behaviour& behaviour,
        const size_type s,
        const MaterialDataManagerInitializer& i)
        : s0(behaviour, s, i.s0),
          s1(behaviour, s, i.s1),
          n(s),
          K_stride(getTangentOperatorArraySize(behaviour)),
          b(behaviour) {
      if (i.K.empty()) {
        constexpr const auto zero = real{0};
        this->K_values.resize(this->n * this->K_stride, zero);
        this->K = mgis::span<real>(this->K_values);
      } else {
        if (i.K.size() != this->n * this->K_stride) {
          mgis::raise(
              "MaterialDataManager::MaterialDataManager: "
              "the consistent tangent operator has not been allocated "
              "properly");
        }
        this->K = i.K;
      }
    }  // end of Behaviour::Behaviour

    MaterialDataManager::~MaterialDataManager() = default;

    void update(MaterialDataManager& m) {
      std::fill(m.K.begin(), m.K.end(), real{0});
      update_values(m.s0, m.s1);
    }  // end of update

    void revert(MaterialDataManager& m) {
      std::fill(m.K.begin(), m.K.end(), real{0});
      update_values(m.s1, m.s0);
    }  // end of update

  }  // end of namespace behaviour

}  // end of namespace mgis
