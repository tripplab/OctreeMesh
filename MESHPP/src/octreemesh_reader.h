#pragma once

#include "post_msh_reader.h"

#include <istream>

namespace meshpp {

// Reads OctreeMesh geometry data files into the shared meshpp MeshData model.
//
// The current meshpp_apply pipeline only transforms geometry/topology, so the
// material column in each OctreeMesh element row is parsed for validation and
// then intentionally ignored.
class OctreeMeshReader {
 public:
  ValidationReport Read(std::istream& input, MeshData* out_mesh, const PostMshReadProgress* progress = nullptr) const;
};

}  // namespace meshpp
