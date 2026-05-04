#pragma once

#include <array>
#include <cstddef>
#include <istream>
#include <string>
#include <unordered_map>
#include <vector>

namespace mesh2stl {

enum class ExitCode {
  kSuccess = 0,
  kUsageError = 2,
  kParseError = 3,
  kUnsupported = 4,
  kTopologyError = 5,
  kIoError = 6,
};

struct Node {
  std::size_t id;
  std::array<double, 3> xyz;
};

struct HexElement {
  std::size_t id;
  std::array<std::size_t, 8> node_ids;
};

struct MeshData {
  std::vector<Node> nodes;
  std::vector<HexElement> elements;
  std::unordered_map<std::size_t, std::size_t> node_id_to_index;
};

struct ValidationIssue {
  ExitCode code;
  std::string message;
  std::size_t line_number;
};

struct ValidationReport {
  std::vector<ValidationIssue> issues;

  bool ok() const { return issues.empty(); }
  ExitCode first_error_code() const;
};

class MeshReader {
 public:
  ValidationReport Read(std::istream& input, MeshData* out_mesh) const;
};

ValidationReport ValidateReferences(const MeshData& mesh);

}  // namespace mesh2stl
