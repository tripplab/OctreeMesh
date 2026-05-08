#include "post_msh_writer.h"

#include <ostream>

namespace meshpp {

ValidationReport WritePostMsh(const MeshData& mesh, const PostMshWriteOptions& options, std::ostream& output) {
  ValidationReport report;
  if (!output.good()) {
    report.issues.push_back({ExitCode::kIoError, "E_IO: output stream is not writable", 0});
    return report;
  }

  output << "MESH \"" << options.mesh_name << "\" dimension " << options.dimension
         << " ElemType " << options.element_type << " Nnode " << options.nodes_per_element << "\n\n";
  output << "Coordinates\n";
  for (const auto& node : mesh.nodes) {
    output << node.id << " " << node.xyz[0] << " " << node.xyz[1] << " " << node.xyz[2] << "\n";
  }
  output << "End Coordinates\n";

  output << "Elements\n";
  for (const auto& element : mesh.elements) {
    output << element.id;
    for (const auto node_id : element.node_ids) {
      output << " " << node_id;
    }
    output << "\n";
  }
  output << "End Elements\n";

  if (!output.good()) {
    report.issues.push_back({ExitCode::kIoError, "E_IO: failed while writing .post.msh output", 0});
  }
  return report;
}

}  // namespace meshpp
