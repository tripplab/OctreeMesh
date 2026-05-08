#pragma once

#include "post_msh_reader.h"

#include <iosfwd>
#include <string>

namespace meshpp {

struct PostMshWriteOptions {
  std::string mesh_name = "mesh";
  int dimension = 3;
  std::string element_type = "Hexahedra";
  int nodes_per_element = 8;
};

ValidationReport WritePostMsh(const MeshData& mesh, const PostMshWriteOptions& options, std::ostream& output);

}  // namespace meshpp
