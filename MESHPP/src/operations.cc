#include "operations.h"

#include <cmath>
#include <iomanip>
#include <iostream>
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

class StatsOperation : public MeshOperation {
 public:
  const char* Name() const override { return "mesh_stats"; }

  ValidationReport Configure(const std::string& spec) override {
    ValidationReport report;
    if (spec.empty()) {
      return report;
    }

    const std::string kPrefix = "format=";
    if (spec.compare(0, kPrefix.size(), kPrefix) != 0) {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: mesh_stats supports optional format=text (e.g. mesh_stats:format=text)", 0});
      return report;
    }

    const std::string value = spec.substr(kPrefix.size());
    if (value != "text") {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: mesh_stats format must be text", 0});
    }
    return report;
  }

  ValidationReport Apply(MeshData* mesh) const override {
    ValidationReport report;
    if (mesh->nodes.empty()) {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: mesh_stats requires at least one node", 0});
      return report;
    }

    double min_x = mesh->nodes[0].xyz[0];
    double min_y = mesh->nodes[0].xyz[1];
    double min_z = mesh->nodes[0].xyz[2];
    double max_x = min_x;
    double max_y = min_y;
    double max_z = min_z;
    double sum_x = 0.0;
    double sum_y = 0.0;
    double sum_z = 0.0;

    for (const auto& node : mesh->nodes) {
      const double x = node.xyz[0];
      const double y = node.xyz[1];
      const double z = node.xyz[2];

      if (x < min_x) min_x = x;
      if (y < min_y) min_y = y;
      if (z < min_z) min_z = z;
      if (x > max_x) max_x = x;
      if (y > max_y) max_y = y;
      if (z > max_z) max_z = z;

      sum_x += x;
      sum_y += y;
      sum_z += z;
    }

    const double count = static_cast<double>(mesh->nodes.size());
    const double cx = sum_x / count;
    const double cy = sum_y / count;
    const double cz = sum_z / count;
    const double dx = max_x - min_x;
    const double dy = max_y - min_y;
    const double dz = max_z - min_z;
    const double diag = std::sqrt(dx * dx + dy * dy + dz * dz);

    const std::streamsize old_precision = std::cout.precision();
    const auto old_flags = std::cout.flags();
    std::cout << std::fixed << std::setprecision(6);

    std::cout << "mesh.stats.nodes=" << mesh->nodes.size() << "\n";
    std::cout << "mesh.stats.elements=" << mesh->elements.size() << "\n";
    std::cout << "mesh.stats.min.x=" << min_x << "\n";
    std::cout << "mesh.stats.min.y=" << min_y << "\n";
    std::cout << "mesh.stats.min.z=" << min_z << "\n";
    std::cout << "mesh.stats.max.x=" << max_x << "\n";
    std::cout << "mesh.stats.max.y=" << max_y << "\n";
    std::cout << "mesh.stats.max.z=" << max_z << "\n";
    std::cout << "mesh.stats.center.x=" << cx << "\n";
    std::cout << "mesh.stats.center.y=" << cy << "\n";
    std::cout << "mesh.stats.center.z=" << cz << "\n";
    std::cout << "mesh.stats.bbox.dx=" << dx << "\n";
    std::cout << "mesh.stats.bbox.dy=" << dy << "\n";
    std::cout << "mesh.stats.bbox.dz=" << dz << "\n";
    std::cout << "mesh.stats.bbox.diag=" << diag << "\n";

    std::cout.flags(old_flags);
    std::cout.precision(old_precision);
    return report;
  }
};

}  // namespace

std::unique_ptr<MeshOperation> CreateOperation(const std::string& name) {
  if (name == "scale") {
    return std::unique_ptr<MeshOperation>(new ScaleOperation());
  }
  if (name == "translate") {
    return std::unique_ptr<MeshOperation>(new TranslateOperation());
  }
  if (name == "mesh_stats") {
    return std::unique_ptr<MeshOperation>(new StatsOperation());
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
