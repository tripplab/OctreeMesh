#include "operations.h"
#include "post_msh_reader.h"
#include "post_msh_writer.h"
#include "perf.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using meshpp::ExitCode;

namespace {
int ToInt(ExitCode code) { return static_cast<int>(code); }
}

int main(int argc, char** argv) {
  if (argc < 6) {
    std::cout << "Usage: meshpp_apply --in <input.post.msh> --out <output.post.msh> --op <spec> [--op <spec> ...] [--mesh_stats] [--perf_stats]\n";
    return ToInt(ExitCode::kUsageError);
  }

  std::string input_path;
  std::string output_path;
  std::vector<std::string> ops;
  bool perf_stats = false;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--in" && i + 1 < argc) {
      input_path = argv[++i];
    } else if (arg == "--out" && i + 1 < argc) {
      output_path = argv[++i];
    } else if (arg == "--op" && i + 1 < argc) {
      ops.push_back(argv[++i]);
    } else if (arg == "--mesh_stats") {
      ops.push_back("mesh_stats");
    } else if (arg == "--perf_stats") {
      perf_stats = true;
    } else {
      std::cerr << "E_USAGE: unknown or incomplete option: " << arg << "\n";
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
    auto read_report = reader.Read(input, &mesh);
    if (!read_report.ok()) {
      std::cerr << read_report.issues.front().message << "\n";
      return ToInt(read_report.issues.front().code);
    }
  }

  {
    meshpp::ScopedTimer timer(&perf.validate_ms);
    auto ref_report = meshpp::ValidateReferences(mesh);
    if (!ref_report.ok()) {
      std::cerr << ref_report.issues.front().message << "\n";
      return ToInt(ref_report.issues.front().code);
    }
  }

  {
    meshpp::ScopedTimer timer(&perf.operations_ms);
    auto op_report = meshpp::ApplyOperationPipeline(ops, &mesh);
    if (!op_report.ok()) {
      std::cerr << op_report.issues.front().message << "\n";
      return ToInt(op_report.issues.front().code);
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
    options.mesh_name = "meshpp_apply";
    auto write_report = meshpp::WritePostMsh(mesh, options, output);
    if (!write_report.ok()) {
      std::cerr << write_report.issues.front().message << "\n";
      return ToInt(write_report.issues.front().code);
    }
  }

  if (perf_stats) {
    std::cout << "stats.read_ms=" << perf.read_ms << "\n";
    std::cout << "stats.validate_ms=" << perf.validate_ms << "\n";
    std::cout << "stats.operations_ms=" << perf.operations_ms << "\n";
    std::cout << "stats.write_ms=" << perf.write_ms << "\n";
    std::cout << "stats.nodes=" << mesh.nodes.size() << "\n";
    std::cout << "stats.elements=" << mesh.elements.size() << "\n";
  }

  std::cout << "meshpp apply: OK\n";
  return ToInt(ExitCode::kSuccess);
}
