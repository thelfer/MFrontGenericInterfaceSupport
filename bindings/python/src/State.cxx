/*!
 * \file   State.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   06/11/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <boost/python/class.hpp>
#include "MGIS/Behaviour/State.hxx"

void declareState();

void declareState() {
  using mgis::behaviour::Behaviour;
  using mgis::behaviour::State;
  boost::python::class_<State>("State", boost::python::no_init)
      .add_property("stored_energy", &State::stored_energy)
      .add_property("dissipated_energy", &State::dissipated_energy);
}  // end of declareState
