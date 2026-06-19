#include "operations.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <map>
#include <set>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

namespace meshpp {
namespace {

using Vec3 = std::array<double, 3>;
using Vec3i = std::array<int, 3>;
using Mat3 = std::array<std::array<double, 3>, 3>;

struct SymmetricSpectralSystem3x3 {
  Vec3 values;
  Mat3 vectors;
};

double Dot(const Vec3& a, const Vec3& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vec3 Cross(const Vec3& a, const Vec3& b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

double Norm(const Vec3& v) {
  return std::sqrt(Dot(v, v));
}

Vec3 Add(const Vec3& a, const Vec3& b) {
  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

Vec3 Subtract(const Vec3& a, const Vec3& b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

Vec3 ScaleVector(const Vec3& v, double scale) {
  return {v[0] * scale, v[1] * scale, v[2] * scale};
}

Vec3 Normalize(const Vec3& v) {
  return ScaleVector(v, 1.0 / Norm(v));
}

void PrintVectorReport(const std::string& name, const Vec3& v) {
  std::cout << "mesh.align_axes." << name << ".x=" << v[0] << "\n";
  std::cout << "mesh.align_axes." << name << ".y=" << v[1] << "\n";
  std::cout << "mesh.align_axes." << name << ".z=" << v[2] << "\n";
}

void PrintPointReport(const std::string& name, const Node& node) {
  std::cout << "mesh.align_axes." << name << ".id=" << node.id << "\n";
  PrintVectorReport(name, node.xyz);
}

void PrintMatrixReport(const std::string& name, const Mat3& matrix) {
  for (std::size_t r = 0; r < 3; ++r) {
    for (std::size_t c = 0; c < 3; ++c) {
      std::cout << "mesh.align_axes." << name << ".r" << r << c << "=" << matrix[r][c] << "\n";
    }
  }
}


std::vector<std::string> SplitColonTokens(const std::string& spec) {
  std::vector<std::string> tokens;
  std::size_t start = 0;
  while (start <= spec.size()) {
    const std::size_t pos = spec.find(':', start);
    tokens.push_back(spec.substr(start, pos == std::string::npos ? std::string::npos : pos - start));
    if (pos == std::string::npos) {
      break;
    }
    start = pos + 1;
  }
  return tokens;
}

double DeterminantColumns(const Vec3& c0, const Vec3& c1, const Vec3& c2) {
  return c0[0] * (c1[1] * c2[2] - c1[2] * c2[1]) - c1[0] * (c0[1] * c2[2] - c0[2] * c2[1]) +
         c2[0] * (c0[1] * c1[2] - c0[2] * c1[1]);
}

Mat3 IdentityMatrix() {
  return {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
}

Mat3 ZeroMatrix() {
  return {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
}

Vec3 MatrixColumn(const Mat3& m, std::size_t c) {
  return {m[0][c], m[1][c], m[2][c]};
}

void SetMatrixColumn(Mat3* m, std::size_t c, const Vec3& v) {
  for (std::size_t r = 0; r < 3; ++r) {
    (*m)[r][c] = v[r];
  }
}

Mat3 Transpose(const Mat3& m) {
  Mat3 result = ZeroMatrix();
  for (std::size_t r = 0; r < 3; ++r) {
    for (std::size_t c = 0; c < 3; ++c) {
      result[r][c] = m[c][r];
    }
  }
  return result;
}

Mat3 Multiply(const Mat3& a, const Mat3& b) {
  Mat3 result = ZeroMatrix();
  for (std::size_t r = 0; r < 3; ++r) {
    for (std::size_t c = 0; c < 3; ++c) {
      for (std::size_t k = 0; k < 3; ++k) {
        result[r][c] += a[r][k] * b[k][c];
      }
    }
  }
  return result;
}

Vec3 Multiply(const Mat3& m, const Vec3& v) {
  return {Dot(m[0], v), Dot(m[1], v), Dot(m[2], v)};
}

Mat3 OuterProduct(const Vec3& a, const Vec3& b) {
  Mat3 result = ZeroMatrix();
  for (std::size_t r = 0; r < 3; ++r) {
    for (std::size_t c = 0; c < 3; ++c) {
      result[r][c] = a[r] * b[c];
    }
  }
  return result;
}

void AddMatrixInPlace(Mat3* a, const Mat3& b) {
  for (std::size_t r = 0; r < 3; ++r) {
    for (std::size_t c = 0; c < 3; ++c) {
      (*a)[r][c] += b[r][c];
    }
  }
}

double Determinant(const Mat3& m) {
  return DeterminantColumns(MatrixColumn(m, 0), MatrixColumn(m, 1), MatrixColumn(m, 2));
}

double MaxAbsCoeff(const Mat3& m) {
  double result = 0.0;
  for (const auto& row : m) {
    for (double value : row) {
      result = std::max(result, std::fabs(value));
    }
  }
  return result;
}

Mat3 Subtract(const Mat3& a, const Mat3& b) {
  Mat3 result = ZeroMatrix();
  for (std::size_t r = 0; r < 3; ++r) {
    for (std::size_t c = 0; c < 3; ++c) {
      result[r][c] = a[r][c] - b[r][c];
    }
  }
  return result;
}

void NormalizeSign(Vec3* v) {
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

SymmetricSpectralSystem3x3 JacobiEigenDecomposition(Mat3 a) {
  SymmetricSpectralSystem3x3 result{{a[0][0], a[1][1], a[2][2]}, IdentityMatrix()};

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

Mat3 OrthonormalizeColumns(Mat3 m) {
  Vec3 c0 = Normalize(MatrixColumn(m, 0));
  Vec3 c1 = Subtract(MatrixColumn(m, 1), ScaleVector(c0, Dot(MatrixColumn(m, 1), c0)));
  if (Norm(c1) < 1e-12) {
    c1 = std::fabs(c0[0]) < 0.9 ? Vec3{1.0, 0.0, 0.0} : Vec3{0.0, 1.0, 0.0};
    c1 = Subtract(c1, ScaleVector(c0, Dot(c1, c0)));
  }
  c1 = Normalize(c1);
  Vec3 c2 = Cross(c0, c1);
  if (Dot(c2, MatrixColumn(m, 2)) < 0.0) {
    c2 = ScaleVector(c2, -1.0);
  }
  SetMatrixColumn(&m, 0, c0);
  SetMatrixColumn(&m, 1, c1);
  SetMatrixColumn(&m, 2, c2);
  return m;
}

Mat3 KabschRotationFromCovariance(const Mat3& h) {
  Mat3 hth = Multiply(Transpose(h), h);
  SymmetricSpectralSystem3x3 eigen = JacobiEigenDecomposition(hth);
  std::array<int, 3> order{0, 1, 2};
  std::sort(order.begin(), order.end(), [&](int lhs, int rhs) { return eigen.values[lhs] > eigen.values[rhs]; });

  Mat3 v = ZeroMatrix();
  Mat3 u = ZeroMatrix();
  for (std::size_t sorted = 0; sorted < order.size(); ++sorted) {
    const int eig_index = order[sorted];
    Vec3 v_col{eigen.vectors[0][eig_index], eigen.vectors[1][eig_index], eigen.vectors[2][eig_index]};
    v_col = Normalize(v_col);
    SetMatrixColumn(&v, sorted, v_col);

    const double singular_value = std::sqrt(std::max(0.0, eigen.values[eig_index]));
    if (singular_value > 1e-12) {
      SetMatrixColumn(&u, sorted, ScaleVector(Multiply(h, v_col), 1.0 / singular_value));
    }
  }

  u = OrthonormalizeColumns(u);
  v = OrthonormalizeColumns(v);

  Mat3 d = IdentityMatrix();
  if (Determinant(Multiply(v, Transpose(u))) < 0.0) {
    d[2][2] = -1.0;
  }
  return Multiply(Multiply(v, d), Transpose(u));
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
    const auto tokens = SplitColonTokens(spec);
    if (tokens.empty() || tokens[0].empty()) {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: align_axes expects pca, simple, or kabsch (e.g. align_axes:pca, align_axes:simple, or align_axes:kabsch)", 0});
      return report;
    }

    if (tokens[0] == "pca") {
      method_ = Method::kPca;
    } else if (tokens[0] == "simple") {
      method_ = Method::kSimple;
    } else if (tokens[0] == "kabsch") {
      method_ = Method::kKabsch;
    } else {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: align_axes expects pca, simple, or kabsch (e.g. align_axes:pca, align_axes:simple, or align_axes:kabsch)", 0});
      return report;
    }

    for (std::size_t i = 1; i < tokens.size(); ++i) {
      if (tokens[i] == "snap") {
        snap_enabled_ = true;
      } else if (tokens[i] == "global") {
        global_refit_ = true;
      } else {
        report.issues.push_back({ExitCode::kUsageError, "E_USAGE: align_axes unknown modifier: " + tokens[i], 0});
        return report;
      }
    }
    if (snap_enabled_ && method_ != Method::kKabsch) {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: align_axes:snap requires a successful rotation (use kabsch)", 0});
    }
    if (global_refit_ && method_ != Method::kKabsch) {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: align_axes:global requires kabsch", 0});
    }
    return report;
  }

  ValidationReport Apply(MeshData* mesh) const override {
    if (method_ == Method::kSimple) {
      return ApplySimple(mesh);
    }
    if (method_ == Method::kKabsch) {
      return ApplyKabsch(mesh);
    }
    return ApplyPca(mesh);
  }

 private:
  enum class Method { kPca, kSimple, kKabsch };

  struct SnapAxisResult {
    std::vector<long> k;
    double phase = 0.0;
    double max_residual = 0.0;
    double max_plane_spread = 0.0;
    int unique_planes = 0;
    std::vector<double> snapped;
  };

  struct SnapResult {
    SnapAxisResult axis[3];
    int cells_collapsed = 0;
    int cells_nonunit = 0;
    bool ok = false;
  };

  static SnapAxisResult SnapAxis(const std::vector<Node>& nodes, int axis, double L) {
    SnapAxisResult result;
    const std::size_t n = nodes.size();
    result.k.resize(n);
    result.snapped.resize(n);

    double cmin = std::numeric_limits<double>::infinity();
    for (const auto& node : nodes) {
      cmin = std::min(cmin, node.xyz[axis]);
    }

    // Assign every rotated coordinate to the nearest uniform plane index, then
    // average out residual noise into one global lattice phase for this axis.
    double phase_accumulator = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      result.k[i] = std::lround((nodes[i].xyz[axis] - cmin) / L);
      phase_accumulator += nodes[i].xyz[axis] - static_cast<double>(result.k[i]) * L;
    }
    result.phase = n == 0 ? 0.0 : phase_accumulator / static_cast<double>(n);

    std::map<long, std::pair<double, double>> bands;
    for (std::size_t i = 0; i < n; ++i) {
      const double coordinate = nodes[i].xyz[axis];
      const double snapped = result.phase + static_cast<double>(result.k[i]) * L;
      result.snapped[i] = snapped;
      result.max_residual = std::max(result.max_residual, std::fabs(coordinate - snapped));

      auto it = bands.find(result.k[i]);
      if (it == bands.end()) {
        bands[result.k[i]] = {coordinate, coordinate};
      } else {
        it->second.first = std::min(it->second.first, coordinate);
        it->second.second = std::max(it->second.second, coordinate);
      }
    }
    for (const auto& band : bands) {
      result.max_plane_spread = std::max(result.max_plane_spread, band.second.second - band.second.first);
    }
    result.unique_planes = static_cast<int>(bands.size());
    return result;
  }


  static SnapAxisResult SnapAxisVariable(const std::vector<Node>& nodes, int axis, double L) {
    SnapAxisResult result;
    const std::size_t n = nodes.size();
    result.k.resize(n);
    result.snapped.resize(n);
    if (n == 0) {
      return result;
    }

    std::vector<std::pair<double, std::size_t>> coords;
    coords.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
      coords.push_back({nodes[i].xyz[axis], i});
    }
    std::sort(coords.begin(), coords.end());

    const double split_gap = std::max(1e-12, 0.25 * L);
    long plane = 0;
    std::size_t cluster_begin = 0;
    while (cluster_begin < coords.size()) {
      std::size_t cluster_end = cluster_begin + 1;
      double sum = coords[cluster_begin].first;
      while (cluster_end < coords.size() &&
             coords[cluster_end].first - coords[cluster_end - 1].first <= split_gap) {
        sum += coords[cluster_end].first;
        ++cluster_end;
      }
      const double center = sum / static_cast<double>(cluster_end - cluster_begin);
      const double spread = coords[cluster_end - 1].first - coords[cluster_begin].first;
      result.max_plane_spread = std::max(result.max_plane_spread, spread);
      for (std::size_t j = cluster_begin; j < cluster_end; ++j) {
        const std::size_t node_index = coords[j].second;
        result.k[node_index] = plane;
        result.snapped[node_index] = center;
        result.max_residual = std::max(result.max_residual, std::fabs(nodes[node_index].xyz[axis] - center));
      }
      ++plane;
      cluster_begin = cluster_end;
    }
    result.unique_planes = static_cast<int>(plane);
    return result;
  }

  bool snap_enabled_ = false;
  bool global_refit_ = false;

  struct RefitResult {
    Mat3 R = IdentityMatrix();
    double L[3] = {0.0, 0.0, 0.0};
    int iterations = 0;
    int edges_used = 0;
    int last_bucket_changes = 0;
    bool nonuniform[3] = {false, false, false};
    double resid_mean_deg = 0.0;
    double resid_p95_deg = 0.0;
    double resid_max_deg = 0.0;
    double min_dominance = 1.0;
    bool ok = false;
  };

  struct KabschResult {
    Mat3 R = IdentityMatrix();
    Vec3 centroidP{0.0, 0.0, 0.0};
    Vec3 centroidQ{0.0, 0.0, 0.0};
    double L = 0.0;
    double rmsd = 0.0;
    double rmsd_norm = 0.0;
    double det = 0.0;
    double ortho_err_max = 0.0;
    std::array<Vec3i, 8> bits;
    std::array<int, 8> ids{};
    bool ok = false;
  };

  struct KabschSeed {
    KabschResult result;
    std::size_t element_index = 0;
    std::size_t element_id = 0;
    std::size_t candidates_checked = 0;
    std::size_t candidates_rejected = 0;
    bool ok = false;
  };

  static bool TryBuildKabschSeed(const MeshData& mesh, const HexElement& element, KabschResult* result) {
    constexpr double kTolerance = 1e-9;
    std::array<Vec3, 8> corners;
    std::array<int, 8> ids{};
    for (std::size_t i = 0; i < element.node_ids.size(); ++i) {
      const auto node_it = mesh.node_id_to_index.find(element.node_ids[i]);
      if (node_it == mesh.node_id_to_index.end()) {
        return false;
      }
      const Node& node = mesh.nodes[node_it->second];
      corners[i] = Vec3{node.xyz[0], node.xyz[1], node.xyz[2]};
      ids[i] = static_cast<int>(node.id);
    }

    const Vec3& P0 = corners[0];
    const Vec3 e1 = Subtract(corners[4], P0);
    const Vec3 e2 = Subtract(corners[1], P0);
    const Vec3 e3 = Subtract(corners[3], P0);
    const double e1_norm = Norm(e1);
    const double e2_norm = Norm(e2);
    const double e3_norm = Norm(e3);
    if (e1_norm < kTolerance || e2_norm < kTolerance || e3_norm < kTolerance) {
      return false;
    }

    const Vec3 w1 = Normalize(e1);
    if (std::fabs(Dot(w1, Normalize(e2))) > 0.999) {
      return false;
    }
    const Vec3 w2_seed = Subtract(e2, ScaleVector(w1, Dot(e2, w1)));
    if (Norm(w2_seed) < kTolerance) {
      return false;
    }
    const Vec3 w2 = Normalize(w2_seed);
    const Vec3 w3 = Cross(w1, w2);

    KabschResult candidate;
    candidate.L = (e1_norm + e2_norm + e3_norm) / 3.0;
    candidate.ids = ids;

    std::array<Vec3, 8> P;
    std::array<Vec3, 8> Q;
    std::set<std::array<int, 3>> seen;
    for (std::size_t i = 0; i < corners.size(); ++i) {
      const Vec3 r = Subtract(corners[i], P0);
      const int bx = (Dot(r, w1) > 0.5 * candidate.L) ? 1 : 0;
      const int by = (Dot(r, w2) > 0.5 * candidate.L) ? 1 : 0;
      const int bz = (Dot(r, w3) > 0.5 * candidate.L) ? 1 : 0;
      candidate.bits[i] = Vec3i{bx, by, bz};
      P[i] = corners[i];
      Q[i] = Vec3{bx * candidate.L, by * candidate.L, bz * candidate.L};
      seen.insert({bx, by, bz});
    }
    if (seen.size() != 8) {
      return false;
    }

    for (std::size_t i = 0; i < P.size(); ++i) {
      candidate.centroidP = Add(candidate.centroidP, P[i]);
      candidate.centroidQ = Add(candidate.centroidQ, Q[i]);
    }
    candidate.centroidP = ScaleVector(candidate.centroidP, 1.0 / 8.0);
    candidate.centroidQ = ScaleVector(candidate.centroidQ, 1.0 / 8.0);

    Mat3 H = ZeroMatrix();
    for (std::size_t i = 0; i < P.size(); ++i) {
      AddMatrixInPlace(&H, OuterProduct(Subtract(P[i], candidate.centroidP), Subtract(Q[i], candidate.centroidQ)));
    }

    candidate.R = KabschRotationFromCovariance(H);
    candidate.det = Determinant(candidate.R);

    double sse = 0.0;
    for (std::size_t i = 0; i < P.size(); ++i) {
      const Vec3 error = Subtract(Multiply(candidate.R, Subtract(P[i], candidate.centroidP)), Subtract(Q[i], candidate.centroidQ));
      sse += Dot(error, error);
    }
    candidate.rmsd = std::sqrt(sse / 8.0);
    candidate.rmsd_norm = candidate.rmsd / candidate.L;
    candidate.ortho_err_max = MaxAbsCoeff(Subtract(Multiply(candidate.R, Transpose(candidate.R)), IdentityMatrix()));
    candidate.ok = true;
    *result = candidate;
    return true;
  }

  static KabschSeed FindKabschSeed(const MeshData& mesh) {
    KabschSeed seed;
    for (std::size_t i = 0; i < mesh.elements.size(); ++i) {
      ++seed.candidates_checked;
      KabschResult result;
      if (TryBuildKabschSeed(mesh, mesh.elements[i], &result)) {
        seed.result = result;
        seed.element_index = i;
        seed.element_id = mesh.elements[i].id;
        seed.ok = true;
        return seed;
      }
      ++seed.candidates_rejected;
    }
    return seed;
  }


  static double RotationDeltaRadians(const Mat3& a, const Mat3& b) {
    const Mat3 rel = Multiply(Transpose(a), b);
    const double value = std::max(-1.0, std::min(1.0, (rel[0][0] + rel[1][1] + rel[2][2] - 1.0) * 0.5));
    return std::acos(value);
  }

  static int DominantAxis(const Vec3& v) {
    int axis = 0;
    double best = std::fabs(v[0]);
    for (int c = 1; c < 3; ++c) {
      if (std::fabs(v[c]) > best) {
        best = std::fabs(v[c]);
        axis = c;
      }
    }
    return axis;
  }

  struct Hist {
    double lo;
    double hi;
    int nb;
    std::vector<uint64_t> bin;
    uint64_t under = 0;
    uint64_t over = 0;
    uint64_t n = 0;
    double sum = 0.0;
    double exact_min = 1e300;
    double exact_max = -1e300;

    Hist(double lo_, double hi_, int nb_) : lo(lo_), hi(hi_), nb(nb_), bin(static_cast<std::size_t>(nb_), 0) {}

    void Add(double x) {
      ++n;
      sum += x;
      exact_min = std::min(exact_min, x);
      exact_max = std::max(exact_max, x);
      if (x < lo) {
        ++under;
        return;
      }
      if (x >= hi) {
        ++over;
        return;
      }
      int i = static_cast<int>((x - lo) / (hi - lo) * static_cast<double>(nb));
      if (i < 0) i = 0;
      if (i >= nb) i = nb - 1;
      ++bin[static_cast<std::size_t>(i)];
    }

    double Mean() const { return n ? sum / static_cast<double>(n) : std::numeric_limits<double>::quiet_NaN(); }

    double Quantile(double q) const {
      if (n == 0) return std::numeric_limits<double>::quiet_NaN();
      const uint64_t target = static_cast<uint64_t>(std::ceil(q * static_cast<double>(n)));
      uint64_t cum = under;
      if (cum >= target) return lo;
      for (int i = 0; i < nb; ++i) {
        cum += bin[static_cast<std::size_t>(i)];
        if (cum >= target) return lo + (hi - lo) * (static_cast<double>(i) + 0.5) / static_cast<double>(nb);
      }
      return exact_max;
    }
  };

  ValidationReport RefineRotationGlobal(const MeshData& mesh, const Mat3& seed_R, double seed_L, RefitResult* result) const {
    constexpr int kRefitMaxIters = 5;
    constexpr double kRefitAngleTolRad = 1e-7;
    constexpr double kEdgeAxisMinDominance = 0.80;
    constexpr double kMinEdgeLen = 1e-9;
    constexpr int kHexEdges[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                                       {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};

    ValidationReport report;
    if (mesh.elements.empty()) {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: align_axes:global requires hex connectivity", 0});
      return report;
    }

    std::cout << "meshpp apply: align_axes:kabsch:global refitting rotation from all cell edge directions\n";

    const std::size_t edge_count = mesh.elements.size() * 12;
    std::vector<int8_t> curr_axis(edge_count, -1);
    std::vector<int8_t> prev_axis(edge_count, -1);
    Mat3 rcur = seed_R;
    int iters = 0;
    int bucket_changes = 0;
    int edges_used = 0;
    double min_dom_final = 1.0;

    for (int iter = 0; iter < kRefitMaxIters; ++iter) {
      Mat3 H = ZeroMatrix();
      double min_dom = 1.0;
      int used = 0;
      int changes = 0;

      for (std::size_t c = 0; c < mesh.elements.size(); ++c) {
        const auto& h = mesh.elements[c];
        for (int e = 0; e < 12; ++e) {
          const std::size_t idx = c * 12 + static_cast<std::size_t>(e);
          const auto a_it = mesh.node_id_to_index.find(h.node_ids[kHexEdges[e][0]]);
          const auto b_it = mesh.node_id_to_index.find(h.node_ids[kHexEdges[e][1]]);
          if (a_it == mesh.node_id_to_index.end() || b_it == mesh.node_id_to_index.end()) {
            curr_axis[idx] = -1;
            continue;
          }
          const Vec3 ev = Subtract(mesh.nodes[b_it->second].xyz, mesh.nodes[a_it->second].xyz);
          const double len = Norm(ev);
          if (len < kMinEdgeLen) {
            curr_axis[idx] = -1;
            continue;
          }
          const Vec3 u = ScaleVector(ev, 1.0 / len);
          const Vec3 f = Multiply(rcur, u);
          const int ax = DominantAxis(f);
          min_dom = std::min(min_dom, std::fabs(f[ax]) / std::max(Norm(f), kMinEdgeLen));
          const double sign = f[ax] >= 0.0 ? 1.0 : -1.0;
          Vec3 target{0.0, 0.0, 0.0};
          target[ax] = sign;
          AddMatrixInPlace(&H, OuterProduct(u, target));
          const int8_t sid = static_cast<int8_t>(ax + (sign < 0.0 ? 3 : 0));
          if (iter > 0 && prev_axis[idx] >= 0 && prev_axis[idx] != sid) {
            ++changes;
          }
          curr_axis[idx] = sid;
          ++used;
        }
      }

      if (used == 0) {
        report.issues.push_back({ExitCode::kUsageError, "E_USAGE: align_axes:global no usable edges", 0});
        return report;
      }

      const Mat3 rnew = KabschRotationFromCovariance(H);
      const double dtheta = RotationDeltaRadians(rnew, rcur);
      prev_axis = curr_axis;
      rcur = rnew;
      min_dom_final = min_dom;
      bucket_changes = changes;
      edges_used = used;
      ++iters;
      std::cout << "meshpp apply: align_axes:kabsch:global iteration=" << (iter + 1) << " edges=" << used << " bucket_changes=" << changes
                << " dtheta_rad=" << dtheta << "\n";
      if (iter > 0 && dtheta < kRefitAngleTolRad) {
        break;
      }
    }

    result->R = rcur;
    result->iterations = iters;
    result->edges_used = edges_used;
    result->last_bucket_changes = bucket_changes;
    result->min_dominance = min_dom_final;

    if (min_dom_final < kEdgeAxisMinDominance) {
      report.issues.push_back({ExitCode::kUsageError,
                               "E_USAGE: align_axes:global edges exceed ~37deg off-axis after refit; mesh is not a single coherently-oriented grid (multi-orientation/body-fitted input)",
                               0});
      return report;
    }

    Hist hL[3] = {Hist(0.25 * seed_L, 4.0 * seed_L, 4000), Hist(0.25 * seed_L, 4.0 * seed_L, 4000),
                   Hist(0.25 * seed_L, 4.0 * seed_L, 4000)};
    Hist hR(0.0, 10.0, 4000);

    for (std::size_t c = 0; c < mesh.elements.size(); ++c) {
      const auto& h = mesh.elements[c];
      for (int e = 0; e < 12; ++e) {
        const int8_t sid = curr_axis[c * 12 + static_cast<std::size_t>(e)];
        if (sid < 0) continue;
        const int ax = sid % 3;
        const double sign = sid >= 3 ? -1.0 : 1.0;
        const auto a_it = mesh.node_id_to_index.find(h.node_ids[kHexEdges[e][0]]);
        const auto b_it = mesh.node_id_to_index.find(h.node_ids[kHexEdges[e][1]]);
        if (a_it == mesh.node_id_to_index.end() || b_it == mesh.node_id_to_index.end()) continue;
        const Vec3 ev = Subtract(mesh.nodes[b_it->second].xyz, mesh.nodes[a_it->second].xyz);
        const double len = Norm(ev);
        if (len < kMinEdgeLen) continue;
        hL[ax].Add(len);
        const Vec3 ru = Multiply(result->R, ScaleVector(ev, 1.0 / len));
        const double d = std::max(-1.0, std::min(1.0, ru[ax] * sign));
        hR.Add(std::acos(d) * 180.0 / 3.141592653589793238462643383279502884);
      }
    }

    for (int a = 0; a < 3; ++a) {
      result->L[a] = hL[a].n ? hL[a].Quantile(0.5) : seed_L;
    }
    result->resid_mean_deg = hR.Mean();
    result->resid_p95_deg = hR.Quantile(0.95);
    result->resid_max_deg = hR.n ? hR.exact_max : 0.0;

    for (int a = 0; a < 3; ++a) {
      const double med = hL[a].Quantile(0.5);
      const double q1 = hL[a].Quantile(0.25);
      const double q3 = hL[a].Quantile(0.75);
      result->nonuniform[a] = (med > 0.0) && ((q3 - q1) / med > 0.05);
      if (result->nonuniform[a]) {
        std::cerr << "W_QUALITY: align_axes:global non-uniform cell sizes on axis " << a << "; route to :snapvar\n";
      }
    }

    if (result->resid_max_deg > 1.0) {
      std::cerr << "W_QUALITY: align_axes:global residual max > 1deg\n";
    }
    result->ok = true;
    return report;
  }

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
    const std::array<std::array<double, 3>, 3> rotation{{u1, u2, u3}};
    const std::array<double, 3> ru1{Dot(rotation[0], u1), Dot(rotation[1], u1), Dot(rotation[2], u1)};
    const std::array<double, 3> ru2{Dot(rotation[0], u2), Dot(rotation[1], u2), Dot(rotation[2], u2)};
    const std::array<double, 3> ru3{Dot(rotation[0], u3), Dot(rotation[1], u3), Dot(rotation[2], u3)};
    const double determinant = DeterminantColumns(rotation[0], rotation[1], rotation[2]);

    const std::streamsize old_precision = std::cout.precision();
    const auto old_flags = std::cout.flags();
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "mesh.align_axes.method=simple\n";
    PrintPointReport("p0", p0_node);
    PrintPointReport("p1", p1_node);
    PrintPointReport("p2", p2_node);
    PrintPointReport("p3", p3_node);
    std::cout << "mesh.align_axes.e1.length=" << e1_length << "\n";
    std::cout << "mesh.align_axes.e2.length=" << e2_length << "\n";
    std::cout << "mesh.align_axes.e3.length=" << e3_length << "\n";
    PrintVectorReport("e1", e1);
    PrintVectorReport("e2", e2);
    PrintVectorReport("e3", e3);
    std::cout << "mesh.align_axes.u1.length=" << Norm(u1) << "\n";
    std::cout << "mesh.align_axes.u2.length=" << Norm(u2) << "\n";
    std::cout << "mesh.align_axes.u3.length=" << Norm(u3) << "\n";
    PrintVectorReport("u1", u1);
    PrintVectorReport("u2", u2);
    PrintVectorReport("u3", u3);
    PrintVectorReport("u1_cross_u2", u1_cross_u2);
    PrintMatrixReport("matrix", rotation);
    PrintVectorReport("R_u1", ru1);
    PrintVectorReport("R_u2", ru2);
    PrintVectorReport("R_u3", ru3);
    std::cout << "mesh.align_axes.det=" << determinant << "\n";
    std::cout << "mesh.align_axes.u1_dot_u2=" << u1_u2 << "\n";
    std::cout << "mesh.align_axes.u1_dot_u3=" << u1_u3 << "\n";
    std::cout << "mesh.align_axes.u2_dot_u3=" << u2_u3 << "\n";
    std::cout.flush();
    std::cout.flags(old_flags);
    std::cout.precision(old_precision);

    if (std::fabs(u1_u2) > kTolerance || std::fabs(u1_u3) > kTolerance || std::fabs(u2_u3) > kTolerance || cross_alignment > kTolerance) {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: align_axes:simple edge neighbors do not form an orthonormal right-handed frame", 0});
      return report;
    }

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

    const std::streamsize old_export_precision = std::cout.precision();
    const auto old_export_flags = std::cout.flags();
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "mesh.align_axes.snap=off\n";
    std::cout << "mesh.align_axes.elements.exported=" << mesh->elements.size() << "\n";
    std::cout << "mesh.align_axes.elements.dropped=" << (original_elements - mesh->elements.size()) << "\n";
    std::cout << "mesh.align_axes.nodes.exported=" << mesh->nodes.size() << "\n";
    std::cout << "mesh.align_axes.nodes.dropped=" << (original_nodes - mesh->nodes.size()) << "\n";
    std::cout.flags(old_export_flags);
    std::cout.precision(old_export_precision);
    return {};
  }


  ValidationReport ApplyKabsch(MeshData* mesh) const {
    constexpr std::size_t kMaxElements = 10;
    ValidationReport report;
    const std::size_t original_nodes = mesh->nodes.size();
    const std::size_t original_elements = mesh->elements.size();
    if (mesh->elements.empty()) {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: align_axes:kabsch requires at least one element", 0});
      return report;
    }

    if (global_refit_) {
      std::cout << "meshpp apply: align_axes:kabsch computing seed rotation, then global refit over " << mesh->elements.size() << " elements\n";
    }

    KabschSeed seed = FindKabschSeed(*mesh);
    if (!seed.ok) {
      report.issues.push_back({ExitCode::kUsageError, "E_USAGE: align_axes:kabsch no usable hexahedral seed found in current mesh", 0});
      return report;
    }
    KabschResult result = seed.result;
    const auto ids = result.ids;

    RefitResult refit;
    if (global_refit_) {
      auto refit_report = RefineRotationGlobal(*mesh, result.R, result.L, &refit);
      if (!refit_report.ok()) {
        return refit_report;
      }
      result.R = refit.R;
      result.det = Determinant(result.R);
      result.ortho_err_max = MaxAbsCoeff(Subtract(Multiply(result.R, Transpose(result.R)), IdentityMatrix()));
    } else {
      refit.L[0] = result.L;
      refit.L[1] = result.L;
      refit.L[2] = result.L;
    }

    if (!global_refit_ && mesh->elements.size() > kMaxElements) {
      std::cout << "meshpp apply: align_axes:kabsch seed-only export filter keeping first " << kMaxElements << " elements\n";
      mesh->elements.resize(kMaxElements);
    }

    const std::streamsize old_precision = std::cout.precision();
    const auto old_flags = std::cout.flags();
    std::cout << std::fixed << std::setprecision(12);
    std::cout << "mesh.align_axes.method=kabsch\n";
    std::cout << "mesh.align_axes.seed.element.id=" << seed.element_id << "\n";
    std::cout << "mesh.align_axes.seed.element.index=" << seed.element_index << "\n";
    std::cout << "mesh.align_axes.seed.candidates_checked=" << seed.candidates_checked << "\n";
    std::cout << "mesh.align_axes.seed.candidates_rejected=" << seed.candidates_rejected << "\n";
    std::cout << "mesh.align_axes.L=" << result.L << "\n";
    std::cout << "mesh.align_axes.centroidP.x=" << result.centroidP[0] << "\n";
    std::cout << "mesh.align_axes.centroidP.y=" << result.centroidP[1] << "\n";
    std::cout << "mesh.align_axes.centroidP.z=" << result.centroidP[2] << "\n";
    std::cout << "mesh.align_axes.centroidQ.x=" << result.centroidQ[0] << "\n";
    std::cout << "mesh.align_axes.centroidQ.y=" << result.centroidQ[1] << "\n";
    std::cout << "mesh.align_axes.centroidQ.z=" << result.centroidQ[2] << "\n";
    for (int r = 0; r < 3; ++r) {
      for (int c = 0; c < 3; ++c) {
        std::cout << "mesh.align_axes.matrix.r" << r << c << "=" << result.R[r][c] << "\n";
      }
    }
    std::cout << "mesh.align_axes.det=" << result.det << "\n";
    std::cout << "mesh.align_axes.ortho_err_max=" << result.ortho_err_max << "\n";
    std::cout << "mesh.align_axes.rmsd=" << result.rmsd << "\n";
    std::cout << "mesh.align_axes.rmsd_norm=" << result.rmsd_norm << "\n";
    for (std::size_t i = 0; i < ids.size(); ++i) {
      std::cout << "mesh.align_axes.corner." << ids[i] << ".bits=" << result.bits[i][0] << result.bits[i][1] << result.bits[i][2] << "\n";
    }
    std::cout << "mesh.align_axes.global=" << (global_refit_ ? "on" : "off") << "\n";
    if (global_refit_) {
      std::cout << "mesh.align_axes.global.model=single_orientation\n";
      std::cout << "mesh.align_axes.global.iterations=" << refit.iterations << "\n";
      std::cout << "mesh.align_axes.global.edges_used=" << refit.edges_used << "\n";
      std::cout << "mesh.align_axes.global.bucket_changes_final=" << refit.last_bucket_changes << "\n";
      std::cout << "mesh.align_axes.global.min_dominance=" << refit.min_dominance << "\n";
      std::cout << "mesh.align_axes.global.resid_mean_deg=" << refit.resid_mean_deg << "\n";
      std::cout << "mesh.align_axes.global.resid_p95_deg=" << refit.resid_p95_deg << "\n";
      std::cout << "mesh.align_axes.global.resid_max_deg=" << refit.resid_max_deg << "\n";
      std::cout << "mesh.align_axes.global.L.x=" << refit.L[0] << "\n";
      std::cout << "mesh.align_axes.global.L.y=" << refit.L[1] << "\n";
      std::cout << "mesh.align_axes.global.L.z=" << refit.L[2] << "\n";
      std::cout << "mesh.align_axes.global.nonuniform.x=" << (refit.nonuniform[0] ? 1 : 0) << "\n";
      std::cout << "mesh.align_axes.global.nonuniform.y=" << (refit.nonuniform[1] ? 1 : 0) << "\n";
      std::cout << "mesh.align_axes.global.nonuniform.z=" << (refit.nonuniform[2] ? 1 : 0) << "\n";
    }
    std::cout.flags(old_flags);
    std::cout.precision(old_precision);

    if (result.rmsd_norm > 0.02) {
      std::cerr << "W_QUALITY: align_axes:kabsch input deviates from ideal cube\n";
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
        const Vec3 p{node.xyz[0], node.xyz[1], node.xyz[2]};
        const Vec3 rotated = Multiply(result.R, p);
        node.xyz[0] = rotated[0];
        node.xyz[1] = rotated[1];
        node.xyz[2] = rotated[2];
        kept_nodes.push_back(node);
      }
    }

    std::unordered_map<std::size_t, std::size_t> node_id_to_index;
    node_id_to_index.reserve(kept_nodes.size());
    for (std::size_t i = 0; i < kept_nodes.size(); ++i) {
      node_id_to_index[kept_nodes[i].id] = i;
    }

    SnapResult snap;
    if (snap_enabled_) {
      if (result.L <= 0.0) {
        report.issues.push_back({ExitCode::kUsageError, "E_USAGE: align_axes:snap requires a valid L from kabsch", 0});
        return report;
      }

      // Snap only after the Kabsch rotation has been applied to every exported
      // node.  Global refit uses per-axis L values, and non-uniform global
      // axes route to per-plane clustering instead of phase+k*L.
      for (int axis = 0; axis < 3; ++axis) {
        snap.axis[axis] = (global_refit_ && refit.nonuniform[axis])
                              ? SnapAxisVariable(kept_nodes, axis, refit.L[axis])
                              : SnapAxis(kept_nodes, axis, refit.L[axis]);
      }

      const char axis_names[3] = {'x', 'y', 'z'};
      for (int axis = 0; axis < 3; ++axis) {
        const double snap_sep = 0.25 * refit.L[axis];
        if (snap.axis[axis].max_plane_spread >= snap_sep) {
          report.issues.push_back({ExitCode::kUsageError,
                                   std::string("E_USAGE: align_axes:snap axis ") + axis_names[axis] +
                                       " not separable: intra-plane spread exceeds 0.25*L (mesh is not rectilinear at this L)",
                                   0});
          return report;
        }
      }

      for (const auto& element : mesh->elements) {
        for (int axis = 0; axis < 3; ++axis) {
          std::set<long> planes;
          for (std::size_t v = 0; v < element.node_ids.size(); ++v) {
            const auto it = node_id_to_index.find(element.node_ids[v]);
            if (it != node_id_to_index.end()) {
              planes.insert(snap.axis[axis].k[it->second]);
            }
          }
          if (planes.size() < 2) {
            ++snap.cells_collapsed;
          } else if (planes.size() > 2) {
            ++snap.cells_nonunit;
          }
        }
      }
      if (snap.cells_collapsed > 0) {
        report.issues.push_back({ExitCode::kUsageError,
                                 "E_USAGE: align_axes:snap would collapse " + std::to_string(snap.cells_collapsed) +
                                     " cell-faces to zero volume (wrong L or non-grid input)",
                                 0});
        return report;
      }
      if (snap.cells_nonunit > 0) {
        std::cerr << "W_QUALITY: align_axes:snap " << snap.cells_nonunit << " cells span >2 planes on an axis\n";
      }

      for (std::size_t i = 0; i < kept_nodes.size(); ++i) {
        for (int axis = 0; axis < 3; ++axis) {
          kept_nodes[i].xyz[axis] = snap.axis[axis].snapped.empty()
                                         ? snap.axis[axis].phase + static_cast<double>(snap.axis[axis].k[i]) * refit.L[axis]
                                         : snap.axis[axis].snapped[i];
        }
      }
      snap.ok = true;
    }

    mesh->nodes = std::move(kept_nodes);
    mesh->node_id_to_index = std::move(node_id_to_index);

    const std::streamsize old_export_precision = std::cout.precision();
    const auto old_export_flags = std::cout.flags();
    std::cout << std::fixed << std::setprecision(6);
    if (snap_enabled_) {
      const double max_residual = std::max(snap.axis[0].max_residual, std::max(snap.axis[1].max_residual, snap.axis[2].max_residual));
      std::cout << "mesh.align_axes.snap=on\n";
      std::cout << "mesh.align_axes.snap.L=" << result.L << "\n";
      std::cout << "mesh.align_axes.snap.L.x=" << refit.L[0] << "\n";
      std::cout << "mesh.align_axes.snap.L.y=" << refit.L[1] << "\n";
      std::cout << "mesh.align_axes.snap.L.z=" << refit.L[2] << "\n";
      std::cout << "mesh.align_axes.snap.phase.x=" << snap.axis[0].phase << "\n";
      std::cout << "mesh.align_axes.snap.phase.y=" << snap.axis[1].phase << "\n";
      std::cout << "mesh.align_axes.snap.phase.z=" << snap.axis[2].phase << "\n";
      std::cout << "mesh.align_axes.snap.unique_planes.x=" << snap.axis[0].unique_planes << "\n";
      std::cout << "mesh.align_axes.snap.unique_planes.y=" << snap.axis[1].unique_planes << "\n";
      std::cout << "mesh.align_axes.snap.unique_planes.z=" << snap.axis[2].unique_planes << "\n";
      std::cout << "mesh.align_axes.snap.max_residual.x=" << snap.axis[0].max_residual << "\n";
      std::cout << "mesh.align_axes.snap.max_residual.y=" << snap.axis[1].max_residual << "\n";
      std::cout << "mesh.align_axes.snap.max_residual.z=" << snap.axis[2].max_residual << "\n";
      std::cout << "mesh.align_axes.snap.max_residual=" << max_residual << "\n";
      std::cout << "mesh.align_axes.snap.max_plane_spread.x=" << snap.axis[0].max_plane_spread << "\n";
      std::cout << "mesh.align_axes.snap.max_plane_spread.y=" << snap.axis[1].max_plane_spread << "\n";
      std::cout << "mesh.align_axes.snap.max_plane_spread.z=" << snap.axis[2].max_plane_spread << "\n";
      std::cout << "mesh.align_axes.snap.cells_collapsed=" << snap.cells_collapsed << "\n";
      std::cout << "mesh.align_axes.snap.cells_nonunit=" << snap.cells_nonunit << "\n";
    } else {
      std::cout << "mesh.align_axes.snap=off\n";
    }
    std::cout << "mesh.align_axes.elements.exported=" << mesh->elements.size() << "\n";
    std::cout << "mesh.align_axes.elements.dropped=" << (original_elements - mesh->elements.size()) << "\n";
    std::cout << "mesh.align_axes.nodes.exported=" << mesh->nodes.size() << "\n";
    std::cout << "mesh.align_axes.nodes.dropped=" << (original_nodes - mesh->nodes.size()) << "\n";
    std::cout.flags(old_export_flags);
    std::cout.precision(old_export_precision);
    return report;
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

    SymmetricSpectralSystem3x3 eigen = JacobiEigenDecomposition(covariance);
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
    std::cout << "mesh.align_axes.snap=off\n";
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
    std::cerr << "meshpp apply: starting operation " << spec << "\n";
    auto apply_report = op->Apply(mesh);
    if (!apply_report.ok()) {
      return apply_report;
    }
  }
  return report;
}

}  // namespace meshpp
