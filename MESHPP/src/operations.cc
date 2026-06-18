#include "operations.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

namespace meshpp {
namespace {


struct EigenSystem3x3 {
  std::array<double, 3> values;
  std::array<std::array<double, 3>, 3> vectors;
};

double Dot(const std::array<double, 3>& a, const std::array<double, 3>& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

std::array<double, 3> Cross(const std::array<double, 3>& a, const std::array<double, 3>& b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

double Norm(const std::array<double, 3>& v) {
  return std::sqrt(Dot(v, v));
}

std::array<double, 3> Subtract(const std::array<double, 3>& a, const std::array<double, 3>& b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

std::array<double, 3> ScaleVector(const std::array<double, 3>& v, double scale) {
  return {v[0] * scale, v[1] * scale, v[2] * scale};
}

void PrintVectorReport(const std::string& name, const std::array<double, 3>& v) {
  std::cout << "mesh.align_axes." << name << ".x=" << v[0] << "\n";
  std::cout << "mesh.align_axes." << name << ".y=" << v[1] << "\n";
  std::cout << "mesh.align_axes." << name << ".z=" << v[2] << "\n";
}

void PrintPointReport(const std::string& name, const Node& node) {
  std::cout << "mesh.align_axes." << name << ".id=" << node.id << "\n";
  PrintVectorReport(name, node.xyz);
}

double DeterminantColumns(const std::array<double, 3>& c0, const std::array<double, 3>& c1, const std::array<double, 3>& c2) {
  return c0[0] * (c1[1] * c2[2] - c1[2] * c2[1]) - c1[0] * (c0[1] * c2[2] - c0[2] * c2[1]) +
         c2[0] * (c0[1] * c1[2] - c0[2] * c1[1]);
}

void NormalizeSign(std::array<double, 3>* v) {
  std::size_t max_index = 0;
  double max_abs = std::fabs((*v)[0]);
  for (std::size_t i = 1; i < v->size(); ++i) {
    const double candidate = std::fabs((*v)[i]);
    if (candidate > max_abs) {
      max_abs = candidate;
      max_index = i;
    }
  }
  if ((*v)[max_index] < 0.0) {
    for (double& component : *v) {
      component = -component;
    }
  }
}

EigenSystem3x3 JacobiEigenDecomposition(std::array<std::array<double, 3>, 3> a) {
  EigenSystem3x3 result{{a[0][0], a[1][1], a[2][2]}, {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}}};

  for (int sweep = 0; sweep < 64; ++sweep) {
    int p = 0;
    int q = 1;
    double max_offdiag = std::fabs(a[0][1]);
    if (std::fabs(a[0][2]) > max_offdiag) {
      max_offdiag = std::fabs(a[0][2]);
      p = 0;
      q = 2;
    }
    if (std::fabs(a[1][2]) > max_offdiag) {
      max_offdiag = std::fabs(a[1][2]);
      p = 1;
      q = 2;
    }
    if (max_offdiag < 1e-12) {
      break;
    }

    const double app = a[p][p];
    const double aqq = a[q][q];
    const double apq = a[p][q];
    const double tau = (aqq - app) / (2.0 * apq);
    const double t = (tau >= 0.0 ? 1.0 : -1.0) / (std::fabs(tau) + std::sqrt(1.0 + tau * tau));
    const double c = 1.0 / std::sqrt(1.0 + t * t);
    const double s = t * c;

    for (int k = 0; k < 3; ++k) {
      if (k == p || k == q) {
        continue;
      }
      const double akp = a[k][p];
      const double akq = a[k][q];
      a[k][p] = c * akp - s * akq;
      a[p][k] = a[k][p];
      a[k][q] = s * akp + c * akq;
      a[q][k] = a[k][q];
    }

    a[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
    a[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
    a[p][q] = 0.0;
    a[q][p] = 0.0;

    for (int k = 0; k < 3; ++k) {
      const double vkp = result.vectors[k][p];
      const double vkq = result.vectors[k][q];
      result.vectors[k][p] = c * vkp - s * vkq;
      result.vectors[k][q] = s * vkp + c * vkq;
    }
  }

  result.values = {a[0][0], a[1][1], a[2][2]};
  return result;
}

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

class RotateOperation : public MeshOperation {
 public:
  const char* Name() const override { return "rotate"; }

  ValidationReport Configure(const std::string& spec) override {
    ValidationReport report;
    std::istringstream iss(spec);
    char comma1 = 0;
    char comma2 = 0;
    char comma3 = 0;
    if (!(iss >> ux_ >> comma1 >> uy_ >> comma2 >> uz_ >> comma3 >> angle_degrees_) || comma1 != ',' || comma2 != ',' || comma3 != ',') {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: rotate expects ux,uy,uz,degrees (e.g. rotate:0,0,1,90)", 0});
      return report;
    }
    char extra = 0;
    if (iss >> extra) {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: rotate expects ux,uy,uz,degrees (e.g. rotate:0,0,1,90)", 0});
      return report;
    }
    if (!std::isfinite(ux_) || !std::isfinite(uy_) || !std::isfinite(uz_) || !std::isfinite(angle_degrees_)) {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: rotate values must be finite numbers", 0});
      return report;
    }

    const double axis_length = std::sqrt(ux_ * ux_ + uy_ * uy_ + uz_ * uz_);
    if (!std::isfinite(axis_length)) {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: rotate axis vector length must be finite", 0});
      return report;
    }
    if (axis_length == 0.0) {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: rotate axis vector must be non-zero", 0});
      return report;
    }

    ux_ /= axis_length;
    uy_ /= axis_length;
    uz_ /= axis_length;
    return report;
  }

  ValidationReport Apply(MeshData* mesh) const override {
    constexpr double kPi = 3.141592653589793238462643383279502884;
    const double radians = angle_degrees_ * kPi / 180.0;
    const double c = std::cos(radians);
    const double s = std::sin(radians);
    const double one_minus_c = 1.0 - c;

    for (auto& node : mesh->nodes) {
      const double x = node.xyz[0];
      const double y = node.xyz[1];
      const double z = node.xyz[2];
      const double dot = ux_ * x + uy_ * y + uz_ * z;
      const double cross_x = uy_ * z - uz_ * y;
      const double cross_y = uz_ * x - ux_ * z;
      const double cross_z = ux_ * y - uy_ * x;

      node.xyz[0] = x * c + cross_x * s + ux_ * dot * one_minus_c;
      node.xyz[1] = y * c + cross_y * s + uy_ * dot * one_minus_c;
      node.xyz[2] = z * c + cross_z * s + uz_ * dot * one_minus_c;
    }
    return {};
  }

 private:
  double ux_ = 0.0;
  double uy_ = 0.0;
  double uz_ = 1.0;
  double angle_degrees_ = 0.0;
};


class AlignAxesOperation : public MeshOperation {
 public:
  const char* Name() const override { return "align_axes"; }

  ValidationReport Configure(const std::string& spec) override {
    ValidationReport report;
    if (spec == "pca") {
      method_ = Method::kPca;
    } else if (spec == "simple") {
      method_ = Method::kSimple;
    } else {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: align_axes expects pca or simple (e.g. align_axes:pca or align_axes:simple)", 0});
    }
    return report;
  }

  ValidationReport Apply(MeshData* mesh) const override {
    if (method_ == Method::kSimple) {
      return ApplySimple(mesh);
    }
    return ApplyPca(mesh);
  }

 private:
  enum class Method { kPca, kSimple };

  ValidationReport ApplySimple(MeshData* mesh) const {
    constexpr std::size_t kMaxElements = 10;
    constexpr double kTolerance = 1e-9;
    ValidationReport report;
    const std::size_t original_nodes = mesh->nodes.size();
    const std::size_t original_elements = mesh->elements.size();

    if (mesh->elements.empty()) {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: align_axes:simple requires at least one element", 0});
      return report;
    }

    if (mesh->elements.size() > kMaxElements) {
      mesh->elements.resize(kMaxElements);
    }

    const HexElement first_element = mesh->elements.front();
    const auto p0_it = mesh->node_id_to_index.find(first_element.node_ids[0]);
    if (p0_it == mesh->node_id_to_index.end()) {
      report.issues.push_back({ExitCode::kTopologyError, "E_TOPOLOGY: align_axes:simple first element references missing first node", 0});
      return report;
    }
    const Node p0_node = mesh->nodes[p0_it->second];

    struct NeighborCandidate {
      std::size_t element_offset;
      Node node;
      double distance;
    };

    std::vector<NeighborCandidate> candidates;
    candidates.reserve(first_element.node_ids.size() - 1);
    for (std::size_t i = 1; i < first_element.node_ids.size(); ++i) {
      const auto node_it = mesh->node_id_to_index.find(first_element.node_ids[i]);
      if (node_it == mesh->node_id_to_index.end()) {
        report.issues.push_back({ExitCode::kTopologyError, "E_TOPOLOGY: align_axes:simple first element references missing neighbor node", 0});
        return report;
      }
      const Node candidate_node = mesh->nodes[node_it->second];
      candidates.push_back({i, candidate_node, Norm(Subtract(candidate_node.xyz, p0_node.xyz))});
    }

    std::stable_sort(candidates.begin(), candidates.end(), [](const NeighborCandidate& lhs, const NeighborCandidate& rhs) {
      if (lhs.distance != rhs.distance) {
        return lhs.distance < rhs.distance;
      }
      return lhs.element_offset < rhs.element_offset;
    });

    if (candidates.size() < 3 || candidates[0].distance <= kTolerance || candidates[1].distance <= kTolerance || candidates[2].distance <= kTolerance) {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: align_axes:simple requires three non-zero edge neighbors", 0});
      return report;
    }

    Node p1_node = candidates[0].node;
    Node p2_node = candidates[1].node;
    Node p3_node = candidates[2].node;
    std::array<double, 3> e1 = Subtract(p1_node.xyz, p0_node.xyz);
    std::array<double, 3> e2 = Subtract(p2_node.xyz, p0_node.xyz);
    std::array<double, 3> e3 = Subtract(p3_node.xyz, p0_node.xyz);
    double e1_length = Norm(e1);
    double e2_length = Norm(e2);
    double e3_length = Norm(e3);
    std::array<double, 3> u1 = ScaleVector(e1, 1.0 / e1_length);
    std::array<double, 3> u2 = ScaleVector(e2, 1.0 / e2_length);
    std::array<double, 3> u3 = ScaleVector(e3, 1.0 / e3_length);

    std::array<double, 3> u1_cross_u2 = Cross(u1, u2);
    if (Dot(u1_cross_u2, u3) < 0.0) {
      std::swap(p2_node, p3_node);
      std::swap(e2, e3);
      std::swap(e2_length, e3_length);
      u2 = ScaleVector(e2, 1.0 / e2_length);
      u3 = ScaleVector(e3, 1.0 / e3_length);
      u1_cross_u2 = Cross(u1, u2);
    }

    const double u1_u2 = Dot(u1, u2);
    const double u1_u3 = Dot(u1, u3);
    const double u2_u3 = Dot(u2, u3);
    const double cross_alignment = Norm(Subtract(u3, u1_cross_u2));
    if (std::fabs(u1_u2) > kTolerance || std::fabs(u1_u3) > kTolerance || std::fabs(u2_u3) > kTolerance || cross_alignment > kTolerance) {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: align_axes:simple edge neighbors do not form an orthonormal right-handed frame", 0});
      return report;
    }

    const std::array<std::array<double, 3>, 3> rotation{{u1, u2, u3}};
    const double determinant = DeterminantColumns(rotation[0], rotation[1], rotation[2]);
    if (std::fabs(determinant - 1.0) > kTolerance) {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: align_axes:simple rotation matrix determinant is not +1", 0});
      return report;
    }

    std::unordered_set<std::size_t> referenced_node_ids;
    referenced_node_ids.reserve(mesh->elements.size() * 8);
    for (const auto& element : mesh->elements) {
      for (std::size_t node_id : element.node_ids) {
        referenced_node_ids.insert(node_id);
      }
    }

    std::vector<Node> kept_nodes;
    kept_nodes.reserve(referenced_node_ids.size());
    for (auto node : mesh->nodes) {
      if (referenced_node_ids.find(node.id) != referenced_node_ids.end()) {
        const std::array<double, 3> p{node.xyz[0], node.xyz[1], node.xyz[2]};
        node.xyz[0] = Dot(rotation[0], p);
        node.xyz[1] = Dot(rotation[1], p);
        node.xyz[2] = Dot(rotation[2], p);
        kept_nodes.push_back(node);
      }
    }

    std::unordered_map<std::size_t, std::size_t> node_id_to_index;
    node_id_to_index.reserve(kept_nodes.size());
    for (std::size_t i = 0; i < kept_nodes.size(); ++i) {
      node_id_to_index[kept_nodes[i].id] = i;
    }

    mesh->nodes = std::move(kept_nodes);
    mesh->node_id_to_index = std::move(node_id_to_index);

    const std::array<double, 3> ru1{Dot(rotation[0], u1), Dot(rotation[1], u1), Dot(rotation[2], u1)};
    const std::array<double, 3> ru2{Dot(rotation[0], u2), Dot(rotation[1], u2), Dot(rotation[2], u2)};
    const std::array<double, 3> ru3{Dot(rotation[0], u3), Dot(rotation[1], u3), Dot(rotation[2], u3)};

    const std::streamsize old_precision = std::cout.precision();
    const auto old_flags = std::cout.flags();
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "mesh.align_axes.method=simple\n";
    PrintPointReport("p0", p0_node);
    PrintPointReport("p1", p1_node);
    PrintPointReport("p2", p2_node);
    PrintPointReport("p3", p3_node);
    PrintVectorReport("e1", e1);
    PrintVectorReport("e2", e2);
    PrintVectorReport("e3", e3);
    std::cout << "mesh.align_axes.e1.length=" << e1_length << "\n";
    std::cout << "mesh.align_axes.e2.length=" << e2_length << "\n";
    std::cout << "mesh.align_axes.e3.length=" << e3_length << "\n";
    PrintVectorReport("u1", u1);
    PrintVectorReport("u2", u2);
    PrintVectorReport("u3", u3);
    std::cout << "mesh.align_axes.u1.length=" << Norm(u1) << "\n";
    std::cout << "mesh.align_axes.u2.length=" << Norm(u2) << "\n";
    std::cout << "mesh.align_axes.u3.length=" << Norm(u3) << "\n";
    PrintVectorReport("u1_cross_u2", u1_cross_u2);
    for (std::size_t r = 0; r < 3; ++r) {
      for (std::size_t c = 0; c < 3; ++c) {
        std::cout << "mesh.align_axes.matrix.r" << r << c << "=" << rotation[r][c] << "\n";
      }
    }
    PrintVectorReport("R_u1", ru1);
    PrintVectorReport("R_u2", ru2);
    PrintVectorReport("R_u3", ru3);
    std::cout << "mesh.align_axes.det=" << determinant << "\n";
    std::cout << "mesh.align_axes.u1_dot_u2=" << u1_u2 << "\n";
    std::cout << "mesh.align_axes.u1_dot_u3=" << u1_u3 << "\n";
    std::cout << "mesh.align_axes.u2_dot_u3=" << u2_u3 << "\n";
    std::cout << "mesh.align_axes.elements.exported=" << mesh->elements.size() << "\n";
    std::cout << "mesh.align_axes.elements.dropped=" << (original_elements - mesh->elements.size()) << "\n";
    std::cout << "mesh.align_axes.nodes.exported=" << mesh->nodes.size() << "\n";
    std::cout << "mesh.align_axes.nodes.dropped=" << (original_nodes - mesh->nodes.size()) << "\n";
    std::cout.flags(old_flags);
    std::cout.precision(old_precision);
    return {};
  }

  ValidationReport ApplyPca(MeshData* mesh) const {
    ValidationReport report;
    if (mesh->nodes.size() < 3) {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: align_axes:pca requires at least three nodes", 0});
      return report;
    }

    std::array<double, 3> center{0.0, 0.0, 0.0};
    for (const auto& node : mesh->nodes) {
      center[0] += node.xyz[0];
      center[1] += node.xyz[1];
      center[2] += node.xyz[2];
    }
    const double node_count = static_cast<double>(mesh->nodes.size());
    center[0] /= node_count;
    center[1] /= node_count;
    center[2] /= node_count;

    std::array<std::array<double, 3>, 3> covariance{{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
    for (const auto& node : mesh->nodes) {
      const std::array<double, 3> d{node.xyz[0] - center[0], node.xyz[1] - center[1], node.xyz[2] - center[2]};
      for (std::size_t r = 0; r < 3; ++r) {
        for (std::size_t c = r; c < 3; ++c) {
          covariance[r][c] += d[r] * d[c] / node_count;
        }
      }
    }
    covariance[1][0] = covariance[0][1];
    covariance[2][0] = covariance[0][2];
    covariance[2][1] = covariance[1][2];

    EigenSystem3x3 eigen = JacobiEigenDecomposition(covariance);
    std::array<int, 3> order{0, 1, 2};
    std::sort(order.begin(), order.end(), [&](int lhs, int rhs) { return eigen.values[lhs] > eigen.values[rhs]; });

    const double largest_abs = std::max(std::fabs(eigen.values[order[0]]), 1.0);
    if ((std::fabs(eigen.values[order[0]] - eigen.values[order[1]]) / largest_abs) < 1e-9 ||
        (std::fabs(eigen.values[order[1]] - eigen.values[order[2]]) / largest_abs) < 1e-9) {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: align_axes:pca principal axes are ambiguous for this mesh", 0});
      return report;
    }

    std::array<double, 3> v0{eigen.vectors[0][order[0]], eigen.vectors[1][order[0]], eigen.vectors[2][order[0]]};
    std::array<double, 3> v1{eigen.vectors[0][order[1]], eigen.vectors[1][order[1]], eigen.vectors[2][order[1]]};
    std::array<double, 3> v2{eigen.vectors[0][order[2]], eigen.vectors[1][order[2]], eigen.vectors[2][order[2]]};
    NormalizeSign(&v0);
    NormalizeSign(&v1);
    NormalizeSign(&v2);

    // Column order maps the dominant PCA axis v0 to global Z, v1 to global X, and v2 to global Y.
    std::array<std::array<double, 3>, 3> basis{{v1, v2, v0}};
    if (DeterminantColumns(basis[0], basis[1], basis[2]) < 0.0) {
      for (double& component : basis[2]) {
        component = -component;
      }
    }

    std::array<std::array<double, 3>, 3> rotation{{{basis[0][0], basis[0][1], basis[0][2]},
                                                   {basis[1][0], basis[1][1], basis[1][2]},
                                                   {basis[2][0], basis[2][1], basis[2][2]}}};

    for (auto& node : mesh->nodes) {
      const std::array<double, 3> p{node.xyz[0], node.xyz[1], node.xyz[2]};
      node.xyz[0] = Dot(rotation[0], p);
      node.xyz[1] = Dot(rotation[1], p);
      node.xyz[2] = Dot(rotation[2], p);
    }

    const std::streamsize old_precision = std::cout.precision();
    const auto old_flags = std::cout.flags();
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "mesh.align_axes.method=pca\n";
    const std::array<double, 3> sorted_eigenvalues{eigen.values[order[0]], eigen.values[order[1]], eigen.values[order[2]]};
    const std::array<std::array<double, 3>, 3> sorted_eigenvectors{v0, v1, v2};
    for (std::size_t i = 0; i < sorted_eigenvalues.size(); ++i) {
      std::cout << "mesh.align_axes.eigenvalue." << i << "=" << sorted_eigenvalues[i] << "\n";
      std::cout << "mesh.align_axes.eigenvector." << i << ".x=" << sorted_eigenvectors[i][0] << "\n";
      std::cout << "mesh.align_axes.eigenvector." << i << ".y=" << sorted_eigenvectors[i][1] << "\n";
      std::cout << "mesh.align_axes.eigenvector." << i << ".z=" << sorted_eigenvectors[i][2] << "\n";
    }
    for (std::size_t r = 0; r < 3; ++r) {
      for (std::size_t c = 0; c < 3; ++c) {
        std::cout << "mesh.align_axes.matrix.r" << r << c << "=" << rotation[r][c] << "\n";
      }
    }
    std::cout.flags(old_flags);
    std::cout.precision(old_precision);
    return report;
  }

  Method method_ = Method::kPca;
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

class CylinderOperation : public MeshOperation {
 public:
  const char* Name() const override { return "cylinder"; }

  ValidationReport Configure(const std::string& spec) override {
    ValidationReport report;
    std::istringstream iss(spec);
    if (!(iss >> radius_) || !std::isfinite(radius_) || radius_ <= 0.0) {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: cylinder expects positive radius in angstroms, e.g. cylinder:35", 0});
      return report;
    }
    char extra = 0;
    if (iss >> extra) {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: cylinder expects positive radius in angstroms, e.g. cylinder:35", 0});
    }
    return report;
  }

  ValidationReport Apply(MeshData* mesh) const override {
    const double radius2 = radius_ * radius_;
    std::unordered_set<std::size_t> candidate_node_ids;
    candidate_node_ids.reserve(mesh->nodes.size());

    for (const auto& node : mesh->nodes) {
      const double x = node.xyz[0];
      const double y = node.xyz[1];
      const double z = node.xyz[2];
      const double r2 = x * x + y * y;
      if (r2 <= radius2 && z > 0.0) {
        candidate_node_ids.insert(node.id);
      }
    }

    std::vector<HexElement> kept_elements;
    kept_elements.reserve(mesh->elements.size());
    std::unordered_set<std::size_t> referenced_node_ids;
    referenced_node_ids.reserve(candidate_node_ids.size());
    for (const auto& element : mesh->elements) {
      bool keep = true;
      for (std::size_t node_id : element.node_ids) {
        if (candidate_node_ids.find(node_id) == candidate_node_ids.end()) {
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
    kept_nodes.reserve(referenced_node_ids.size());
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

    std::cout << "mesh.cylinder.radius=" << radius_ << "\n";
    std::cout << "mesh.cylinder.nodes_kept=" << mesh->nodes.size() << "\n";
    std::cout << "mesh.cylinder.nodes_dropped=" << (original_nodes - mesh->nodes.size()) << "\n";
    std::cout << "mesh.cylinder.elements_kept=" << mesh->elements.size() << "\n";
    std::cout << "mesh.cylinder.elements_dropped=" << (original_elements - mesh->elements.size()) << "\n";
    return {};
  }

 private:
  double radius_ = 0.0;
};

}  // namespace

std::unique_ptr<MeshOperation> CreateOperation(const std::string& name) {
  if (name == "scale") {
    return std::unique_ptr<MeshOperation>(new ScaleOperation());
  }
  if (name == "translate") {
    return std::unique_ptr<MeshOperation>(new TranslateOperation());
  }
  if (name == "rotate") {
    return std::unique_ptr<MeshOperation>(new RotateOperation());
  }
  if (name == "align_axes") {
    return std::unique_ptr<MeshOperation>(new AlignAxesOperation());
  }
  if (name == "mesh_stats") {
    return std::unique_ptr<MeshOperation>(new StatsOperation());
  }
  if (name == "octet") {
    return std::unique_ptr<MeshOperation>(new OctetOperation());
  }
  if (name == "cylinder") {
    return std::unique_ptr<MeshOperation>(new CylinderOperation());
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
