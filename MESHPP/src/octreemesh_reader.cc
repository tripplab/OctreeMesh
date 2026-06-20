#include "octreemesh_reader.h"

#include <cctype>
#include <sstream>
#include <string>

namespace meshpp {
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

std::string StripOctreeComment(const std::string& line) {
  const std::size_t comment = line.find(';');
  if (comment == std::string::npos) {
    return Trim(line);
  }
  return Trim(line.substr(0, comment));
}

void EmitProgress(const PostMshReadProgress* progress, int percent) {
  if (progress == nullptr || progress->output == nullptr) {
    return;
  }
  *progress->output << "meshpp read: progress=" << percent << "%\n";
}

bool ReadSingleSize(const std::string& line, std::size_t* value) {
  std::istringstream iss(line);
  long long parsed = 0;
  if (!(iss >> parsed) || parsed < 0) {
    return false;
  }
  std::string extra;
  if (iss >> extra) {
    return false;
  }
  *value = static_cast<std::size_t>(parsed);
  return true;
}

bool ReadSingleInt(const std::string& line, int* value) {
  std::istringstream iss(line);
  if (!(iss >> *value)) {
    return false;
  }
  std::string extra;
  return !(iss >> extra);
}

}  // namespace

ValidationReport OctreeMeshReader::Read(std::istream& input, MeshData* out_mesh, const PostMshReadProgress* progress) const {
  ValidationReport report;
  out_mesh->nodes.clear();
  out_mesh->elements.clear();
  out_mesh->node_id_to_index.clear();

  enum class State {
    kSeekingNodes,
    kReadingDimension,
    kReadingNodeCount,
    kReadingNodes,
    kSeekingMesh,
    kReadingElementType,
    kReadingNodesPerElement,
    kReadingElementCount,
    kReadingElements,
    kDone,
  };

  State state = State::kSeekingNodes;
  std::size_t expected_nodes = 0;
  std::size_t expected_elements = 0;
  std::size_t nodes_read = 0;
  std::size_t elements_read = 0;

  std::string line;
  std::size_t line_number = 0;
  int last_reported_percent = -1;
  if (progress != nullptr && progress->output != nullptr) {
    EmitProgress(progress, 0);
    last_reported_percent = 0;
  }

  while (std::getline(input, line)) {
    ++line_number;
    if (progress != nullptr && progress->output != nullptr && progress->total_bytes > 0) {
      const std::streamoff pos = input.tellg();
      if (pos >= 0) {
        int percent = static_cast<int>((pos * 100) / progress->total_bytes);
        if (percent > 100) {
          percent = 100;
        }
        const int rounded_percent = (percent / 5) * 5;
        if (rounded_percent > last_reported_percent) {
          EmitProgress(progress, rounded_percent);
          last_reported_percent = rounded_percent;
        }
      }
    }

    // OctreeMesh files use ';' for both full-line and inline comments. Strip
    // comments before interpreting section headers or numeric rows.
    const std::string t = StripOctreeComment(line);
    if (t.empty()) {
      continue;
    }

    switch (state) {
      case State::kSeekingNodes:
        if (t != "{Nodes}") {
          report.issues.push_back({ExitCode::kParseError, "E_PARSE: missing {Nodes} section", line_number});
          return report;
        }
        state = State::kReadingDimension;
        break;

      case State::kReadingDimension: {
        int dimension = 0;
        if (!ReadSingleInt(t, &dimension)) {
          report.issues.push_back({ExitCode::kParseError, "E_PARSE: invalid OctreeMesh dimension row; expected integer dimension", line_number});
          return report;
        }
        if (dimension != 3) {
          report.issues.push_back({ExitCode::kUnsupported, "E_UNSUPPORTED: only 3D OctreeMesh files are supported", line_number});
          return report;
        }
        state = State::kReadingNodeCount;
        break;
      }

      case State::kReadingNodeCount:
        if (!ReadSingleSize(t, &expected_nodes)) {
          report.issues.push_back({ExitCode::kParseError, "E_PARSE: invalid OctreeMesh node count", line_number});
          return report;
        }
        if (expected_nodes == 0) {
          report.issues.push_back({ExitCode::kParseError, "E_PARSE: OctreeMesh file must contain at least one node", line_number});
          return report;
        }
        out_mesh->nodes.reserve(expected_nodes);
        out_mesh->node_id_to_index.reserve(expected_nodes);
        state = State::kReadingNodes;
        break;

      case State::kReadingNodes: {
        if (t == "{Mesh}") {
          report.issues.push_back({ExitCode::kParseError, "E_PARSE: OctreeMesh file ended before declared node count", line_number});
          return report;
        }
        std::istringstream iss(t);
        Node node;
        node.id = nodes_read + 1;
        if (!(iss >> node.xyz[0] >> node.xyz[1] >> node.xyz[2])) {
          report.issues.push_back({ExitCode::kParseError, "E_PARSE: invalid OctreeMesh node row; expected '<x> <y> <z>'", line_number});
          return report;
        }
        std::string extra;
        if (iss >> extra) {
          report.issues.push_back({ExitCode::kParseError, "E_PARSE: invalid OctreeMesh node row; expected exactly three coordinates", line_number});
          return report;
        }
        out_mesh->node_id_to_index[node.id] = out_mesh->nodes.size();
        out_mesh->nodes.push_back(node);
        ++nodes_read;
        if (nodes_read == expected_nodes) {
          state = State::kSeekingMesh;
        }
        break;
      }

      case State::kSeekingMesh:
        if (t != "{Mesh}") {
          report.issues.push_back({ExitCode::kParseError, "E_PARSE: missing {Mesh} section", line_number});
          return report;
        }
        state = State::kReadingElementType;
        break;

      case State::kReadingElementType: {
        int element_type = 0;
        if (!ReadSingleInt(t, &element_type)) {
          report.issues.push_back({ExitCode::kParseError, "E_PARSE: invalid OctreeMesh element type row; expected integer element type", line_number});
          return report;
        }
        if (element_type != 5) {
          report.issues.push_back({ExitCode::kUnsupported, "E_UNSUPPORTED: only OctreeMesh element type 5 (Hexahedra) is supported", line_number});
          return report;
        }
        state = State::kReadingNodesPerElement;
        break;
      }

      case State::kReadingNodesPerElement: {
        int nodes_per_element = 0;
        if (!ReadSingleInt(t, &nodes_per_element)) {
          report.issues.push_back({ExitCode::kParseError, "E_PARSE: invalid OctreeMesh nodes-per-element row; expected integer count", line_number});
          return report;
        }
        if (nodes_per_element != 8) {
          report.issues.push_back({ExitCode::kUnsupported, "E_UNSUPPORTED: only 8-node OctreeMesh hexahedra are supported", line_number});
          return report;
        }
        state = State::kReadingElementCount;
        break;
      }

      case State::kReadingElementCount:
        if (!ReadSingleSize(t, &expected_elements)) {
          report.issues.push_back({ExitCode::kParseError, "E_PARSE: invalid OctreeMesh element count", line_number});
          return report;
        }
        if (expected_elements == 0) {
          report.issues.push_back({ExitCode::kParseError, "E_PARSE: OctreeMesh file must contain at least one element", line_number});
          return report;
        }
        out_mesh->elements.reserve(expected_elements);
        state = State::kReadingElements;
        break;

      case State::kReadingElements: {
        std::istringstream iss(t);
        std::size_t material_id = 0;
        HexElement elem;
        elem.id = elements_read + 1;
        if (!(iss >> material_id >> elem.node_ids[0] >> elem.node_ids[1] >> elem.node_ids[2] >> elem.node_ids[3] >> elem.node_ids[4] >> elem.node_ids[5] >> elem.node_ids[6] >> elem.node_ids[7])) {
          report.issues.push_back({ExitCode::kParseError, "E_PARSE: invalid OctreeMesh element row; expected '<material> n1..n8'", line_number});
          return report;
        }
        std::string extra;
        if (iss >> extra) {
          report.issues.push_back({ExitCode::kParseError, "E_PARSE: invalid OctreeMesh element row; expected exactly one material id and eight node ids", line_number});
          return report;
        }
        (void)material_id;  // Geometry-only meshpp_apply pipelines intentionally ignore material ids.
        out_mesh->elements.push_back(elem);
        ++elements_read;
        if (elements_read == expected_elements) {
          state = State::kDone;
        }
        break;
      }

      case State::kDone:
        report.issues.push_back({ExitCode::kParseError, "E_PARSE: unexpected content after declared OctreeMesh elements", line_number});
        return report;
    }
  }

  if (nodes_read < expected_nodes) {
    report.issues.push_back({ExitCode::kParseError, "E_PARSE: OctreeMesh file ended before declared node count", line_number});
    return report;
  }
  if (elements_read < expected_elements) {
    report.issues.push_back({ExitCode::kParseError, "E_PARSE: OctreeMesh file ended before declared element count", line_number});
    return report;
  }
  if (state != State::kDone) {
    report.issues.push_back({ExitCode::kParseError, "E_PARSE: incomplete OctreeMesh file", line_number});
    return report;
  }

  if (progress != nullptr && progress->output != nullptr && last_reported_percent < 100) {
    EmitProgress(progress, 100);
  }
  return report;
}

}  // namespace meshpp
