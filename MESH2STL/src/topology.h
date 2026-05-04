#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include "mesh_reader.h"

namespace mesh2stl {

struct QuadFace {
  std::array<std::size_t, 4> node_ids;
};

struct TopologyStats {
  std::size_t total_quad_faces = 0;
  std::size_t boundary_quad_faces = 0;
  std::size_t interior_quad_faces = 0;
  std::size_t non_manifold_quad_faces = 0;
  std::size_t output_quad_faces = 0;
  std::size_t output_triangles = 0;
};

enum class ExtractionMode {
  kSurfaceOnly,
  kAllFaces,
};

struct TopologyResult {
  TopologyStats stats;
  std::vector<QuadFace> output_faces;
  ValidationReport report;
};

TopologyResult AnalyzeTopology(const MeshData& mesh, ExtractionMode mode);

}  // namespace mesh2stl
