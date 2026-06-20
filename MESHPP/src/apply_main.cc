#include "operations.h"
#include "octreemesh_reader.h"
#include "post_msh_reader.h"
#include "post_msh_writer.h"
#include "perf.h"

#include <cctype>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using meshpp::ExitCode;

namespace {
int ToInt(ExitCode code) { return static_cast<int>(code); }

enum class InputFormat {
  kAuto,
  kGid,
  kOctree,
};

std::streamoff StreamSize(std::ifstream* input) {
  const std::streampos original = input->tellg();
  input->seekg(0, std::ios::end);
  const std::streampos end = input->tellg();
  input->seekg(original);
  if (end == std::streampos(-1)) {
    return 0;
  }
  return static_cast<std::streamoff>(end);
}

std::string Trim(const std::string& s) {
  std::size_t start = 0;
  while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start]))) {
    ++start;
  }
  std::size_t end = s.size();
  while (end > start && std::isspace(static_cast<unsigned char>(s[end - 1]))) {
    --end;
  }
  return s.substr(start, end - start);
}

std::string StripAutodetectComment(const std::string& line) {
  std::string candidate = line;
  const std::size_t octree_comment = candidate.find(';');
  if (octree_comment != std::string::npos) {
    candidate = candidate.substr(0, octree_comment);
  }
  return Trim(candidate);
}

bool StartsWith(const std::string& line, const std::string& prefix) { return line.rfind(prefix, 0) == 0; }

bool ParseInputFormat(const std::string& value, InputFormat* format) {
  if (value == "auto") {
    *format = InputFormat::kAuto;
    return true;
  }
  if (value == "gid") {
    *format = InputFormat::kGid;
    return true;
  }
  if (value == "octree") {
    *format = InputFormat::kOctree;
    return true;
  }
  return false;
}

const char* InputFormatName(InputFormat format) {
  switch (format) {
    case InputFormat::kAuto:
      return "auto";
    case InputFormat::kGid:
      return "gid";
    case InputFormat::kOctree:
      return "octree";
  }
  return "unknown";
}

meshpp::ValidationReport DetectInputFormat(std::ifstream* input, InputFormat* detected) {
  meshpp::ValidationReport report;
  input->clear();
  input->seekg(0);

  std::string line;
  std::size_t line_number = 0;
  while (std::getline(*input, line)) {
    ++line_number;
    const std::string t = StripAutodetectComment(line);
    if (t.empty() || StartsWith(t, "#")) {
      continue;
    }
    if (StartsWith(t, "MESH")) {
      *detected = InputFormat::kGid;
      input->clear();
      input->seekg(0);
      return report;
    }
    if (t == "{Nodes}") {
      *detected = InputFormat::kOctree;
      input->clear();
      input->seekg(0);
      return report;
    }
    report.issues.push_back({ExitCode::kParseError, "E_PARSE: cannot auto-detect input format; expected GiD MESH header or OctreeMesh {Nodes} section", line_number});
    input->clear();
    input->seekg(0);
    return report;
  }

  report.issues.push_back({ExitCode::kParseError, "E_PARSE: cannot auto-detect input format from empty input", line_number});
  input->clear();
  input->seekg(0);
  return report;
}

void PrintHelp() {
  std::cout << "meshpp_apply - Apply operation pipeline to mesh files and write GiD .post.msh output\n\n";
  std::cout << "Usage:\n";
  std::cout << "  meshpp_apply --in <input> [--in_format auto|gid|octree] --out <output.post.msh> --op <spec> [--op <spec> ...] [--mesh_stats] [--perf_stats]\n";
  std::cout << "  meshpp_apply -h\n\n";
  std::cout << "Input formats:\n";
  std::cout << "  --in_format auto                Auto-detect GiD .post.msh or OctreeMesh data input (default).\n";
  std::cout << "  --in_format gid                 Read GiD ASCII .post.msh input.\n";
  std::cout << "  --in_format octree              Read OctreeMesh data input with {Nodes}/{Mesh} sections.\n";
  std::cout << "                                  OctreeMesh material ids are ignored; output remains GiD .post.msh.\n\n";
  std::cout << "Available operations (for --op <spec>):\n";
  std::cout << "  scale:<factor>                  Multiply node coordinates by <factor>.\n";
  std::cout << "  translate:<dx>,<dy>,<dz>        Add offsets to node coordinates.\n";
  std::cout << "  rotate:<ux>,<uy>,<uz>,<degrees> Rotate node coordinates around an arbitrary axis.\n";
  std::cout << "  align_axes:pca                  Align PCA principal directions with X/Y/Z axes.\n";
  std::cout << "  align_axes:simple               Export the first ten elements after aligning first-element edge axes to X/Y/Z.\n";
  std::cout << "  align_axes:kabsch               Export first ten elements after Kabsch corner alignment.\n";
  std::cout << "  align_axes:kabsch:global        Seed Kabsch, then refit one global rotation from all hex edge directions.\n";
  std::cout << "  align_axes:kabsch:global:snap   Global Kabsch refit, then snap rotated nodes to the fitted rectilinear lattice.\n";
  std::cout << "  align_axes:kabsch:global:snap:cube\n";
  std::cout << "                                  Snap with one center-anchored cubic cell edge L* preserving bbox volume.\n";
  std::cout << "  align_axes:kabsch:snap          Kabsch-align, then snap rotated nodes to the fitted uniform rectilinear lattice.\n";
  std::cout << "  mesh_stats[:format=text]        Print mesh metrics to stdout.\n";
  std::cout << "  octet:(+|-)x(+|-)y(+|-)z        Keep elements in one octant only.\n";
  std::cout << "  cylinder:<radius>               Keep elements with all nodes in r<=radius and z>0.\n\n";
  std::cout << "Convenience flags:\n";
  std::cout << "  --mesh_stats                    Same as adding --op mesh_stats.\n";
  std::cout << "  --perf_stats                    Print timing and mesh size counters.\n\n";
  std::cout << "Examples:\n";
  std::cout << "  meshpp_apply --in in.post.msh --out out.post.msh --op scale:2.0\n";
  std::cout << "  meshpp_apply --in octreemesh.dat --in_format octree --out out.post.msh --op scale:2.0\n";
  std::cout << "  meshpp_apply --in in.post.msh --out out.post.msh --op translate:1,2,3 --op mesh_stats\n";
  std::cout << "  meshpp_apply --in in.post.msh --out rot.post.msh --op rotate:0,0,1,90\n";
  std::cout << "  meshpp_apply --in in.post.msh --out aligned.post.msh --op align_axes:pca\n";
  std::cout << "  meshpp_apply --in in.post.msh --out first10.post.msh --op align_axes:simple\n";
  std::cout << "  meshpp_apply --in in.post.msh --out aligned.post.msh --op align_axes:kabsch\n";
  std::cout << "  meshpp_apply --in in.post.msh --out snapped.post.msh --op align_axes:kabsch:snap\n";
  std::cout << "  meshpp_apply --in in.post.msh --out global.post.msh --op align_axes:kabsch:global:snap\n";
  std::cout << "  meshpp_apply --in in.post.msh --out cubes.post.msh --op align_axes:kabsch:global:snap:cube\n";
  std::cout << "  meshpp_apply --in in.post.msh --out octant.post.msh --op octet:+x+y-z\n";
  std::cout << "  meshpp_apply --in in.post.msh --out cyl.post.msh --op cylinder:35 --mesh_stats --perf_stats\n";
}
}

int main(int argc, char** argv) {
  if (argc == 2 && std::string(argv[1]) == "-h") {
    PrintHelp();
    return ToInt(ExitCode::kSuccess);
  }

  if (argc < 6) {
    PrintHelp();
    return ToInt(ExitCode::kUsageError);
  }

  std::string input_path;
  std::string output_path;
  std::vector<std::string> ops;
  bool perf_stats = false;
  InputFormat input_format = InputFormat::kAuto;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--in" && i + 1 < argc) {
      input_path = argv[++i];
    } else if (arg == "--out" && i + 1 < argc) {
      output_path = argv[++i];
    } else if (arg == "--in_format" && i + 1 < argc) {
      const std::string value = argv[++i];
      if (!ParseInputFormat(value, &input_format)) {
        std::cerr << "E_USAGE: unsupported --in_format: " << value << " (expected auto, gid, or octree)\n";
        return ToInt(ExitCode::kUsageError);
      }
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
  {
    InputFormat resolved_format = input_format;
    if (resolved_format == InputFormat::kAuto) {
      auto detect_report = DetectInputFormat(&input, &resolved_format);
      if (!detect_report.ok()) {
        std::cerr << detect_report.issues.front().message << "\n";
        return ToInt(detect_report.issues.front().code);
      }
    }

    std::cerr << "meshpp apply: starting read/parse input " << input_path << " (format=" << InputFormatName(resolved_format) << ")\n";
    meshpp::PostMshReadProgress progress;
    progress.output = &std::cerr;
    progress.total_bytes = StreamSize(&input);
    meshpp::ScopedTimer timer(&perf.read_ms);
    meshpp::ValidationReport read_report;
    if (resolved_format == InputFormat::kOctree) {
      meshpp::OctreeMeshReader reader;
      read_report = reader.Read(input, &mesh, &progress);
    } else {
      meshpp::PostMshReader reader;
      read_report = reader.Read(input, &mesh, &progress);
    }
    if (!read_report.ok()) {
      std::cerr << read_report.issues.front().message << "\n";
      return ToInt(read_report.issues.front().code);
    }
  }

  {
    std::cerr << "meshpp apply: starting reference validation\n";
    meshpp::ScopedTimer timer(&perf.validate_ms);
    auto ref_report = meshpp::ValidateReferences(mesh);
    if (!ref_report.ok()) {
      std::cerr << ref_report.issues.front().message << "\n";
      return ToInt(ref_report.issues.front().code);
    }
  }

  {
    std::cerr << "meshpp apply: starting operation pipeline (" << ops.size() << " operations)\n";
    meshpp::ScopedTimer timer(&perf.operations_ms);
    auto op_report = meshpp::ApplyOperationPipeline(ops, &mesh);
    if (!op_report.ok()) {
      std::cerr << op_report.issues.front().message << "\n";
      return ToInt(op_report.issues.front().code);
    }
  }

  std::ofstream output(output_path.c_str());
  std::cerr << "meshpp apply: starting output open " << output_path << "\n";
  static char output_buffer[1 << 20];
  output.rdbuf()->pubsetbuf(output_buffer, sizeof(output_buffer));
  if (!output) {
    std::cerr << "E_IO: cannot open output file: " << output_path << "\n";
    return ToInt(ExitCode::kIoError);
  }

  {
    std::cerr << "meshpp apply: starting write output\n";
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
