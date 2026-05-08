#include "operations.h"

#include <cmath>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <unordered_set>

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

class OctetOperation : public MeshOperation {
 public:
  const char* Name() const override { return "octet"; }

  ValidationReport Configure(const std::string& spec) override {
    ValidationReport report;
    if (spec.size() != 6 || std::tolower(spec[1]) != 'x' || std::tolower(spec[3]) != 'y' || std::tolower(spec[5]) != 'z' ||
        (spec[0] != '+' && spec[0] != '-') || (spec[2] != '+' && spec[2] != '-') || (spec[4] != '+' && spec[4] != '-')) {
      report.issues.push_back({ExitCode::kUsageError,
                               "E_USAGE: octet expects one of (+|-)x(+|-)y(+|-)z (e.g. octet:+x+y+z or octet:+x+y-z)", 0});
      return report;
    }

    keep_positive_x_ = spec[0] == '+';
    keep_positive_y_ = spec[2] == '+';
    keep_positive_z_ = spec[4] == '+';
    spec_ = spec;
    return report;
  }

  ValidationReport Apply(MeshData* mesh) const override {
    std::unordered_set<std::size_t> octet_node_ids;
    octet_node_ids.reserve(mesh->nodes.size());

    auto in_axis = [](double value, bool keep_positive) { return keep_positive ? value >= 0.0 : value < 0.0; };

    for (const auto& node : mesh->nodes) {
      if (in_axis(node.xyz[0], keep_positive_x_) && in_axis(node.xyz[1], keep_positive_y_) && in_axis(node.xyz[2], keep_positive_z_)) {
        octet_node_ids.insert(node.id);
      }
    }

    std::vector<HexElement> kept_elements;
    kept_elements.reserve(mesh->elements.size());
    std::unordered_set<std::size_t> referenced_node_ids;

    for (const auto& element : mesh->elements) {
      bool keep = true;
      for (std::size_t node_id : element.node_ids) {
        if (octet_node_ids.find(node_id) == octet_node_ids.end()) {
          keep = false;
          break;
        }
      }
      if (keep) {
        kept_elements.push_back(element);
        for (std::size_t node_id : element.node_ids) {
          referenced_node_ids.insert(node_id);
        }
      }
    }

    std::vector<Node> kept_nodes;
    kept_nodes.reserve(mesh->nodes.size());
    for (const auto& node : mesh->nodes) {
      if (referenced_node_ids.find(node.id) != referenced_node_ids.end()) {
        kept_nodes.push_back(node);
      }
    }

    std::unordered_map<std::size_t, std::size_t> node_id_to_index;
    node_id_to_index.reserve(kept_nodes.size());
    for (std::size_t i = 0; i < kept_nodes.size(); ++i) {
      node_id_to_index[kept_nodes[i].id] = i;
    }

    const std::size_t original_nodes = mesh->nodes.size();
    const std::size_t original_elements = mesh->elements.size();
    mesh->nodes = std::move(kept_nodes);
    mesh->elements = std::move(kept_elements);
    mesh->node_id_to_index = std::move(node_id_to_index);

    std::cout << "mesh.octet.spec=" << spec_ << "\n";
    std::cout << "mesh.octet.nodes.kept=" << mesh->nodes.size() << "\n";
    std::cout << "mesh.octet.nodes.dropped=" << (original_nodes - mesh->nodes.size()) << "\n";
    std::cout << "mesh.octet.elements.kept=" << mesh->elements.size() << "\n";
    std::cout << "mesh.octet.elements.dropped=" << (original_elements - mesh->elements.size()) << "\n";
    if (mesh->elements.empty()) {
      std::cout << "mesh.octet.warning=selection produced no elements\n";
    }

    return {};
  }

 private:
  bool keep_positive_x_ = true;
  bool keep_positive_y_ = true;
  bool keep_positive_z_ = true;
  std::string spec_;
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
  if (name == "octet") {
    return std::unique_ptr<MeshOperation>(new OctetOperation());
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
