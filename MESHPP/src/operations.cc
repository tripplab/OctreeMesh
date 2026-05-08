#include "operations.h"

#include <sstream>

namespace meshpp {
namespace {

class ScaleOperation : public MeshOperation {
 public:
  const char* Name() const override { return "scale"; }

  ValidationReport Configure(const std::string& spec) override {
    ValidationReport report;
    std::istringstream iss(spec);
    if (!(iss >> factor_)) {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: scale expects numeric factor (e.g. scale:2.0)", 0});
    }
    return report;
  }

  ValidationReport Apply(MeshData* mesh) const override {
    for (auto& node : mesh->nodes) {
      node.xyz[0] *= factor_;
      node.xyz[1] *= factor_;
      node.xyz[2] *= factor_;
    }
    return {};
  }

 private:
  double factor_ = 1.0;
};

class TranslateOperation : public MeshOperation {
 public:
  const char* Name() const override { return "translate"; }

  ValidationReport Configure(const std::string& spec) override {
    ValidationReport report;
    std::istringstream iss(spec);
    char comma1 = 0;
    char comma2 = 0;
    if (!(iss >> dx_ >> comma1 >> dy_ >> comma2 >> dz_) || comma1 != ',' || comma2 != ',') {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: translate expects dx,dy,dz (e.g. translate:1,2,3)", 0});
    }
    return report;
  }

  ValidationReport Apply(MeshData* mesh) const override {
    for (auto& node : mesh->nodes) {
      node.xyz[0] += dx_;
      node.xyz[1] += dy_;
      node.xyz[2] += dz_;
    }
    return {};
  }

 private:
  double dx_ = 0;
  double dy_ = 0;
  double dz_ = 0;
};

}  // namespace

std::unique_ptr<MeshOperation> CreateOperation(const std::string& name) {
  if (name == "scale") {
    return std::unique_ptr<MeshOperation>(new ScaleOperation());
  }
  if (name == "translate") {
    return std::unique_ptr<MeshOperation>(new TranslateOperation());
  }
  return nullptr;
}

ValidationReport ApplyOperationPipeline(const std::vector<std::string>& specs, MeshData* mesh) {
  ValidationReport report;
  for (const auto& spec : specs) {
    const auto pos = spec.find(':');
    const std::string name = pos == std::string::npos ? spec : spec.substr(0, pos);
    const std::string args = pos == std::string::npos ? "" : spec.substr(pos + 1);

    auto op = CreateOperation(name);
    if (!op) {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: unknown operation: " + name, 0});
      return report;
    }
    auto config_report = op->Configure(args);
    if (!config_report.ok()) {
      return config_report;
    }
    auto apply_report = op->Apply(mesh);
    if (!apply_report.ok()) {
      return apply_report;
    }
  }
  return report;
}

}  // namespace meshpp
