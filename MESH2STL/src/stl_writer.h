#pragma once

#include <ostream>
#include <string>

#include "mesh_reader.h"
#include "topology.h"

namespace mesh2stl {

struct StlWriteStats {
  std::size_t written_triangles = 0;
  std::size_t skipped_degenerate_triangles = 0;
};

ValidationReport WriteAsciiStl(const MeshData& mesh,
                               const TopologyResult& topology,
                               const std::string& solid_name,
                               std::ostream& out,
                               StlWriteStats* out_stats);

}  // namespace mesh2stl
