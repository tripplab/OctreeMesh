#pragma once

#include "post_msh_reader.h"

#include <memory>
#include <string>
#include <vector>

namespace meshpp {

class MeshOperation {
 public:
  virtual ~MeshOperation() = default;
  virtual const char* Name() const = 0;
  virtual ValidationReport Configure(const std::string& spec) = 0;
  virtual ValidationReport Apply(MeshData* mesh) const = 0;
};

std::unique_ptr<MeshOperation> CreateOperation(const std::string& name);
ValidationReport ApplyOperationPipeline(const std::vector<std::string>& specs, MeshData* mesh);

}  // namespace meshpp
