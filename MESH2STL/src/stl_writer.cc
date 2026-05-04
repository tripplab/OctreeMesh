#include "stl_writer.h"

#include <cmath>
#include <iomanip>
#include <unordered_map>

namespace mesh2stl {
namespace {

struct Vec3 {
  double x;
  double y;
  double z;
};

Vec3 Sub(const Vec3& a, const Vec3& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
Vec3 Cross(const Vec3& a, const Vec3& b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double Norm(const Vec3& v) { return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z); }

}  // namespace

ValidationReport WriteAsciiStl(const MeshData& mesh,
                               const TopologyResult& topology,
                               const std::string& solid_name,
                               std::ostream& out,
                               StlWriteStats* out_stats) {
  ValidationReport report;
  out_stats->written_triangles = 0;
  out_stats->skipped_degenerate_triangles = 0;

  std::unordered_map<std::size_t, Vec3> coords;
  coords.reserve(mesh.nodes.size());
  for (const auto& node : mesh.nodes) {
    coords[node.id] = {node.xyz[0], node.xyz[1], node.xyz[2]};
  }

  out << std::fixed << std::setprecision(6);
  out << "solid " << solid_name << "\n";
  for (const auto& face : topology.output_faces) {
    const std::array<std::array<std::size_t, 3>, 2> tris = {{
        {{face.node_ids[0], face.node_ids[1], face.node_ids[2]}},
        {{face.node_ids[0], face.node_ids[2], face.node_ids[3]}},
    }};

    for (const auto& tri : tris) {
      const Vec3 p0 = coords[tri[0]];
      const Vec3 p1 = coords[tri[1]];
      const Vec3 p2 = coords[tri[2]];

      const Vec3 e1 = Sub(p1, p0);
      const Vec3 e2 = Sub(p2, p0);
      Vec3 n = Cross(e1, e2);
      const double len = Norm(n);
      if (len <= 1e-14) {
        out_stats->skipped_degenerate_triangles += 1;
        continue;
      }
      n = {n.x / len, n.y / len, n.z / len};

      out << "  facet normal " << n.x << " " << n.y << " " << n.z << "\n";
      out << "    outer loop\n";
      out << "      vertex " << p0.x << " " << p0.y << " " << p0.z << "\n";
      out << "      vertex " << p1.x << " " << p1.y << " " << p1.z << "\n";
      out << "      vertex " << p2.x << " " << p2.y << " " << p2.z << "\n";
      out << "    endloop\n";
      out << "  endfacet\n";
      out_stats->written_triangles += 1;
    }
  }
  out << "endsolid " << solid_name << "\n";
  return report;
}

}  // namespace mesh2stl
