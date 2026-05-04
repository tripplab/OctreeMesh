#include "mesh_reader.h"

#include <algorithm>
#include <cctype>
#include <sstream>

namespace mesh2stl {
namespace {

std::string Trim(const std::string& s) {
  std::size_t start = 0;
  while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start]))) {
    ++start;
  }
  std::size_t end = s.size();
  while (end > start && std::isspace(static_cast<unsigned char>(s[end - 1]))) {
    --end;
  }
  return s.substr(start, end - start);
}

bool StartsWith(const std::string& line, const std::string& prefix) {
  return line.rfind(prefix, 0) == 0;
}

}  // namespace

ExitCode ValidationReport::first_error_code() const {
  if (issues.empty()) {
    return ExitCode::kSuccess;
  }
  return issues.front().code;
}

ValidationReport MeshReader::Read(std::istream& input, MeshData* out_mesh) const {
  ValidationReport report;
  out_mesh->nodes.clear();
  out_mesh->elements.clear();
  out_mesh->node_id_to_index.clear();

  bool in_coordinates = false;
  bool in_elements = false;
  bool saw_mesh_header = false;

  std::string line;
  std::size_t line_number = 0;
  while (std::getline(input, line)) {
    ++line_number;
    const std::string t = Trim(line);
    if (t.empty() || StartsWith(t, "#")) {
      continue;
    }

    if (StartsWith(t, "MESH")) {
      saw_mesh_header = true;
      if (t.find("ElemType Hexahedra") == std::string::npos) {
        report.issues.push_back({ExitCode::kUnsupported,
                                 "E_UNSUPPORTED: only ElemType Hexahedra is supported in v1",
                                 line_number});
        return report;
      }
      continue;
    }

    if (t == "Coordinates") {
      in_coordinates = true;
      in_elements = false;
      continue;
    }

    if (t == "End Coordinates") {
      in_coordinates = false;
      continue;
    }

    if (t == "Elements") {
      in_elements = true;
      in_coordinates = false;
      continue;
    }

    if (t == "End Elements") {
      in_elements = false;
      continue;
    }

    if (in_coordinates) {
      std::istringstream iss(t);
      Node node;
      if (!(iss >> node.id >> node.xyz[0] >> node.xyz[1] >> node.xyz[2])) {
        report.issues.push_back({ExitCode::kParseError,
                                 "E_PARSE: invalid Coordinates row; expected '<id> <x> <y> <z>'",
                                 line_number});
        return report;
      }
      if (out_mesh->node_id_to_index.count(node.id) > 0) {
        report.issues.push_back({ExitCode::kParseError,
                                 "E_PARSE: duplicate node id in Coordinates",
                                 line_number});
        return report;
      }
      out_mesh->node_id_to_index[node.id] = out_mesh->nodes.size();
      out_mesh->nodes.push_back(node);
      continue;
    }

    if (in_elements) {
      std::istringstream iss(t);
      HexElement elem;
      if (!(iss >> elem.id >> elem.node_ids[0] >> elem.node_ids[1] >> elem.node_ids[2] >> elem.node_ids[3] >>
            elem.node_ids[4] >> elem.node_ids[5] >> elem.node_ids[6] >> elem.node_ids[7])) {
        report.issues.push_back({ExitCode::kParseError,
                                 "E_PARSE: invalid Elements row; expected '<id> n1..n8'",
                                 line_number});
        return report;
      }
      out_mesh->elements.push_back(elem);
      continue;
    }

    if (StartsWith(t, "End Coordinate")) {
      report.issues.push_back({ExitCode::kParseError,
                               "E_PARSE: malformed section terminator; expected 'End Coordinates'",
                               line_number});
      return report;
    }
  }

  if (!saw_mesh_header) {
    report.issues.push_back({ExitCode::kParseError, "E_PARSE: missing MESH header", line_number});
    return report;
  }
  if (in_coordinates || in_elements) {
    report.issues.push_back({ExitCode::kParseError,
                             "E_PARSE: unterminated Coordinates or Elements section",
                             line_number});
    return report;
  }
  if (out_mesh->nodes.empty() || out_mesh->elements.empty()) {
    report.issues.push_back({ExitCode::kParseError,
                             "E_PARSE: mesh must contain at least one node and one element",
                             line_number});
    return report;
  }

  return report;
}

ValidationReport ValidateReferences(const MeshData& mesh) {
  ValidationReport report;
  for (const auto& element : mesh.elements) {
    for (std::size_t i = 0; i < element.node_ids.size(); ++i) {
      if (mesh.node_id_to_index.count(element.node_ids[i]) == 0) {
        report.issues.push_back({
            ExitCode::kParseError,
            "E_PARSE: element references undefined node id", 0});
        return report;
      }
    }
  }
  return report;
}

}  // namespace mesh2stl
