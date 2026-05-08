#pragma once

#include <chrono>

namespace meshpp {

struct PerfStats {
  double read_ms = 0;
  double validate_ms = 0;
  double operations_ms = 0;
  double write_ms = 0;
};

class ScopedTimer {
 public:
  explicit ScopedTimer(double* target_ms) : target_ms_(target_ms), start_(std::chrono::steady_clock::now()) {}
  ~ScopedTimer() {
    if (target_ms_ != nullptr) {
      const auto end = std::chrono::steady_clock::now();
      *target_ms_ = std::chrono::duration<double, std::milli>(end - start_).count();
    }
  }

 private:
  double* target_ms_;
  std::chrono::steady_clock::time_point start_;
};

}  // namespace meshpp
