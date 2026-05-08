#include "post_msh_reader.h"
#include "post_msh_writer.h"
#include "perf.h"

#include <fstream>
#include <iostream>

using meshpp::ExitCode;

namespace {
int ToInt(ExitCode code) { return static_cast<int>(code); }
}

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cout << "Usage: meshpp_roundtrip <input.post.msh> <output.post.msh> [--validate] [--stats]\n";
    return ToInt(ExitCode::kUsageError);
  }

  const std::string input_path = argv[1];
  const std::string output_path = argv[2];
  bool validate = false;
  bool stats = false;
  for (int i = 3; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--validate") {
      validate = true;
    } else if (arg == "--stats") {
      stats = true;
    } else {
      std::cerr << "E_USAGE: unknown option: " << arg << "\n";
      return ToInt(ExitCode::kUsageError);
    }
  }

  std::ifstream input(input_path.c_str());
  static char input_buffer[1 << 20];
  input.rdbuf()->pubsetbuf(input_buffer, sizeof(input_buffer));
  if (!input) {
    std::cerr << "E_IO: cannot open input file: " << input_path << "\n";
    return ToInt(ExitCode::kIoError);
  }

  meshpp::PerfStats perf;
  meshpp::MeshData mesh;
  meshpp::PostMshReader reader;
  {
    meshpp::ScopedTimer timer(&perf.read_ms);
    auto report = reader.Read(input, &mesh);
    if (!report.ok()) {
      std::cerr << report.issues.front().message << "\n";
      return ToInt(report.issues.front().code);
    }
  }

  {
    meshpp::ScopedTimer timer(&perf.validate_ms);
    auto ref = meshpp::ValidateReferences(mesh);
    if (!ref.ok()) {
      std::cerr << ref.issues.front().message << "\n";
      return ToInt(ref.issues.front().code);
    }
  }

  std::ofstream output(output_path.c_str());
  static char output_buffer[1 << 20];
  output.rdbuf()->pubsetbuf(output_buffer, sizeof(output_buffer));
  if (!output) {
    std::cerr << "E_IO: cannot open output file: " << output_path << "\n";
    return ToInt(ExitCode::kIoError);
  }
  {
    meshpp::ScopedTimer timer(&perf.write_ms);
    meshpp::PostMshWriteOptions options;
    options.mesh_name = "roundtrip";
    auto wr = meshpp::WritePostMsh(mesh, options, output);
    if (!wr.ok()) {
      std::cerr << wr.issues.front().message << "\n";
      return ToInt(wr.issues.front().code);
    }
  }

  if (validate) {
    std::cout << "nodes: " << mesh.nodes.size() << "\n";
    std::cout << "elements: " << mesh.elements.size() << "\n";
  }

  if (stats) {
    std::cout << "stats.read_ms=" << perf.read_ms << "\n";
    std::cout << "stats.validate_ms=" << perf.validate_ms << "\n";
    std::cout << "stats.write_ms=" << perf.write_ms << "\n";
  }

  return ToInt(ExitCode::kSuccess);
}
