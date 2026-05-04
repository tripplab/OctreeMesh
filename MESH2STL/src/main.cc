#include "mesh_reader.h"
#include "topology.h"
#include "stl_writer.h"

#include <fstream>
#include <iostream>
#include <string>

using mesh2stl::ExitCode;

namespace {

void PrintUsage() {
  std::cout << "Usage: mesh2stl <input_mesh.gid> <output_mesh.stl> [options]\n"
               "Options:\n"
               "  --surface-only   Export boundary faces only (default; not yet implemented)\n"
               "  --all-faces      Export all faces (not yet implemented)\n"
               "  --validate       Run parser and reference validation checks\n"
               "  --units <name>   Log units metadata only\n"
               "  --help           Show this message\n";
}

int ToInt(ExitCode code) { return static_cast<int>(code); }

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    PrintUsage();
    return ToInt(ExitCode::kUsageError);
  }

  std::string input_path = argv[1];
  std::string output_path = argv[2];
  bool validate = false;
  std::string units;
  bool saw_surface_only = false;
  bool saw_all_faces = false;

  for (int i = 3; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--help") {
      PrintUsage();
      return ToInt(ExitCode::kSuccess);
    }
    if (arg == "--validate") {
      validate = true;
      continue;
    }
    if (arg == "--units") {
      if (i + 1 >= argc) {
        std::cerr << "E_USAGE: --units requires a value\n";
        return ToInt(ExitCode::kUsageError);
      }
      units = argv[++i];
      continue;
    }
    if (arg == "--surface-only") {
      saw_surface_only = true;
      continue;
    }
    if (arg == "--all-faces") {
      saw_all_faces = true;
      continue;
    }
    std::cerr << "E_USAGE: unknown option: " << arg << "\n";
    return ToInt(ExitCode::kUsageError);
  }


  if (saw_surface_only && saw_all_faces) {
    std::cerr << "E_USAGE: --surface-only and --all-faces are mutually exclusive\n";
    return ToInt(ExitCode::kUsageError);
  }

  std::ifstream input(input_path.c_str());
  if (!input) {
    std::cerr << "E_IO: cannot open input file: " << input_path << "\n";
    return ToInt(ExitCode::kIoError);
  }

  mesh2stl::MeshData mesh;
  mesh2stl::MeshReader reader;
  mesh2stl::ValidationReport read_report = reader.Read(input, &mesh);
  if (!read_report.ok()) {
    const auto& issue = read_report.issues.front();
    std::cerr << issue.message;
    if (issue.line_number > 0) {
      std::cerr << " (line " << issue.line_number << ")";
    }
    std::cerr << "\n";
    return ToInt(issue.code);
  }

  mesh2stl::ValidationReport ref_report = mesh2stl::ValidateReferences(mesh);
  if (!ref_report.ok()) {
    std::cerr << ref_report.issues.front().message << "\n";
    return ToInt(ref_report.issues.front().code);
  }

  mesh2stl::ExtractionMode mode = mesh2stl::ExtractionMode::kSurfaceOnly;
  if (saw_all_faces) {
    mode = mesh2stl::ExtractionMode::kAllFaces;
  }

  mesh2stl::TopologyResult topology = mesh2stl::AnalyzeTopology(mesh, mode);
  if (!topology.report.ok()) {
    std::cerr << topology.report.issues.front().message << "\n";
    return ToInt(topology.report.issues.front().code);
  }


  std::ofstream output(output_path.c_str());
  if (!output) {
    std::cerr << "E_IO: cannot open output file: " << output_path << "\n";
    return ToInt(ExitCode::kIoError);
  }

  mesh2stl::StlWriteStats write_stats;
  mesh2stl::ValidationReport write_report =
      mesh2stl::WriteAsciiStl(mesh, topology, "mesh2stl", output, &write_stats);
  if (!write_report.ok()) {
    std::cerr << write_report.issues.front().message << "\n";
    return ToInt(write_report.issues.front().code);
  }

  if (validate) {
    std::cout << "Validation report:\n";
    std::cout << "  input: " << input_path << "\n";
    std::cout << "  output: " << output_path << "\n";
    std::cout << "  units: " << (units.empty() ? "(not set)" : units) << "\n";
    std::cout << "  nodes: " << mesh.nodes.size() << "\n";
    std::cout << "  elements: " << mesh.elements.size() << "\n";
    std::cout << "  mode: " << (mode == mesh2stl::ExtractionMode::kSurfaceOnly ? "surface-only" : "all-faces") << "\n";
    std::cout << "  total_quad_faces: " << topology.stats.total_quad_faces << "\n";
    std::cout << "  boundary_quad_faces: " << topology.stats.boundary_quad_faces << "\n";
    std::cout << "  interior_quad_faces: " << topology.stats.interior_quad_faces << "\n";
    std::cout << "  non_manifold_quad_faces: " << topology.stats.non_manifold_quad_faces << "\n";
    std::cout << "  output_quad_faces: " << topology.stats.output_quad_faces << "\n";
    std::cout << "  output_triangles: " << topology.stats.output_triangles << "\n";
    std::cout << "  written_triangles: " << write_stats.written_triangles << "\n";
    std::cout << "  skipped_degenerate_triangles: " << write_stats.skipped_degenerate_triangles << "\n";
    std::cout << "  result: OK\n";
  } else {
    std::cout << "STL written successfully. Use --validate for details.\n";
  }

  return ToInt(ExitCode::kSuccess);
}
