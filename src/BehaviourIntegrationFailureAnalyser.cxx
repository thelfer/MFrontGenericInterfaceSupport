/*!
 * \file   BehaviourIntegrationFailureAnalyser.cxx
 * \brief
 * \author Thomas Helfer
 * \date   21/08/2025
 */

#include <cstddef>
#include <fstream>
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/BehaviourData.hxx"
#include "MGIS/Behaviour/BehaviourDataView.hxx"
#include "MGIS/Behaviour/BehaviourIntegrationFailureAnalyser.hxx"

namespace mgis::behaviour::debug {

  static std::size_t getIdentifier() {
    static auto id = std::size_t{};
    return id++;
  }

  BehaviourIntegrationFailureAnalyser::
      ~BehaviourIntegrationFailureAnalyser() noexcept = default;

  static std::function<
      std::string(std::string_view, std::string_view, std::string_view)>&
  getBehaviourIntegrationFailureAnalysisFileNameGenerator() {
    static std::function<std::string(std::string_view, std::string_view,
                                     std::string_view)>
        g = [](std::string_view b, std::string_view id, std::string_view ext) {
          return std::string{b} + '-' + std::string(id) + "." +
                 std::string{ext};
        };
    return g;
  }

  std::string getBehaviourIntegrationFailureAnalysisFileName(
      std::string_view b, std::string_view ext) {
    const auto& g = getBehaviourIntegrationFailureAnalysisFileNameGenerator();
    return g(b, std::to_string(getIdentifier()), ext);
  }  // end of getBehaviourIntegrationFailureAnalysisFileName

  void setBehaviourIntegrationFailureAnalysisFileNameGenerator(
      std::function<std::string(
          std::string_view, std::string_view, std::string_view)> g) {
    getBehaviourIntegrationFailureAnalysisFileNameGenerator() = g;
  }  // end of setBehaviourIntegrationFailureAnalysisFileNameGenerator

  struct DefaultBehaviourIntegrationFailureAnalyser final
      : BehaviourIntegrationFailureAnalyser {
    //! \brief default constructor
    DefaultBehaviourIntegrationFailureAnalyser() = default;
    //
    void analyse(const Behaviour& b,
                 const BehaviourData& d) const noexcept override {
      std::ofstream f(
          getBehaviourIntegrationFailureAnalysisFileName(b.function, "md"));
      f.precision(14);
      if (f) {
        print_markdown(f, b, d, 0);
      }
    }
    void analyse(const Behaviour& b,
                 const BehaviourDataView& d) const noexcept override {
      // write inputs to output file
      std::ofstream f(
          getBehaviourIntegrationFailureAnalysisFileName(b.function, "md"));
      f.precision(14);
      if (f) {
        print_markdown(f, b, d, 0);
      }
    }
    bool shallCopyBehaviourDataBeforeIntegration() const noexcept override {
      return true;
    }
    ~DefaultBehaviourIntegrationFailureAnalyser() override = default;
  };

  const BehaviourIntegrationFailureAnalyser&
  getDefaultBehaviourIntegrationFailureAnalyser() {
    static DefaultBehaviourIntegrationFailureAnalyser analyser;
    return analyser;
  }  // end of getDefaultBehaviourFailureAnalysisgingOptions

}  // end of namespace mgis::behaviour::debug
