/*!
 * \file   IntegrationFailureDebug.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   21/08/2025
 */

#include <cstdlib>
#include <stdexcept>
#include "MGIS/Behaviour/BehaviourIntegrationFailureAnalyser.hxx"

int main() {
  using namespace mgis::behaviour::debug;
  const auto f1 =
      getBehaviourIntegrationFailureAnalysisFileName("Norton", "md");
  const auto f2 =
      getBehaviourIntegrationFailureAnalysisFileName("Norton", "md");
  if (f1 != "Norton-0.md") {
    throw(std::runtime_error("invalid first file name (got '" + f1 +
                             "', expected 'Norton-0.md')"));
  }
  if (f2 != "Norton-1.md") {
    throw(std::runtime_error("invalid first file name (got '" + f2 +
                             "', expected 'Norton-1.md')"));
  }
  return EXIT_SUCCESS;
}