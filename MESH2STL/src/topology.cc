#include "topology.h"

#include <algorithm>
#include <unordered_map>

namespace mesh2stl {
namespace {

using FaceKey = std::array<std::size_t, 4>;

constexpr int kHexFaceMap[6][4] = {
    {0, 1, 2, 3},
    {4, 5, 6, 7},
    {0, 1, 5, 4},
    {1, 2, 6, 5},
    {2, 3, 7, 6},
    {3, 0, 4, 7},
};

FaceKey CanonicalizeFace(const std::array<std::size_t, 4>& nodes) {
  FaceKey key = nodes;
  std::sort(key.begin(), key.end());
  return key;
}

struct FaceWithOwner {
  QuadFace face;
  FaceKey key;
};

}  // namespace

TopologyResult AnalyzeTopology(const MeshData& mesh, ExtractionMode mode) {
  TopologyResult result;

  struct FaceKeyHash {
    std::size_t operator()(const FaceKey& key) const {
      std::size_t h = 1469598103934665603ull;
      for (std::size_t v : key) {
        h ^= v + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
      }
      return h;
    }
  };

  std::vector<FaceWithOwner> all_faces;
  all_faces.reserve(mesh.elements.size() * 6);
  std::unordered_map<FaceKey, std::size_t, FaceKeyHash> owner_counts;
  owner_counts.reserve(mesh.elements.size() * 6);

  for (const auto& element : mesh.elements) {
    for (int i_face = 0; i_face < 6; ++i_face) {
      QuadFace face{{
          element.node_ids[kHexFaceMap[i_face][0]],
          element.node_ids[kHexFaceMap[i_face][1]],
          element.node_ids[kHexFaceMap[i_face][2]],
          element.node_ids[kHexFaceMap[i_face][3]],
      }};
      FaceKey key = CanonicalizeFace(face.node_ids);
      owner_counts[key] += 1;
      all_faces.push_back({face, key});
      result.stats.total_quad_faces += 1;
    }
  }

  for (const auto& it : owner_counts) {
    if (it.second == 1) {
      result.stats.boundary_quad_faces += 1;
    } else if (it.second == 2) {
      result.stats.interior_quad_faces += 1;
    } else {
      result.stats.non_manifold_quad_faces += 1;
    }
  }

  if (mode == ExtractionMode::kSurfaceOnly) {
    if (result.stats.non_manifold_quad_faces > 0) {
      result.report.issues.push_back(
          {ExitCode::kTopologyError, "E_TOPOLOGY: non-manifold faces detected in --surface-only mode", 0});
      return result;
    }
    for (const auto& fw : all_faces) {
      if (owner_counts[fw.key] == 1) {
        result.output_faces.push_back(fw.face);
      }
    }
  } else {
    for (const auto& fw : all_faces) {
      result.output_faces.push_back(fw.face);
    }
  }

  result.stats.output_quad_faces = result.output_faces.size();
  result.stats.output_triangles = result.stats.output_quad_faces * 2;
  return result;
}

}  // namespace mesh2stl
