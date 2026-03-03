# GTS-CPT Refactoring Plan: C++20 Modernization and Dependency Reduction

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Refactor GTS-CPT to depend only on GTS library, adopt C++20 best practices, and establish clear separation of concerns.

**Architecture:** 
- Core library (`libgts-cpt`) with only GTS dependency
- Optional export module (`libgts-cpt-export`) with Silo dependency
- CLI application (`gts-cpt`) that uses both
- C++20 modules for dependency isolation

**Tech Stack:** C++20, GTS library, GLib, Silo (optional), CMake build system

---

## Overview

### Current State
```
gts-cpt/
├── src/main.cpp         # CLI + GTS I/O + algorithm (tightly coupled)
├── src/cpt.cpp          # Core CPT algorithm
├── src/export.cpp       # Silo/VTK export (forced dependency)
├── src/gtstools.cpp     # GTS utilities
├── include/array3d.h    # DEAD CODE (unused)
└── include/*.h          # Mixed headers
```

### Target Architecture
```
gts-cpt/
├── src/
│   ├── gts-cpt/
│   │   ├── cpt.cppm         # C++20 module: core CPT algorithm
│   │   ├── signed_distance.cppm  # Module: SignedDistance class
│   │   ├── math_utils.cppm  # Module: utility functions
│   │   └── gts_cpt.cppm    # Main module interface
│   ├── gts-cpt-export/     # OPTIONAL export module
│   │   ├── export.cppm      # Silo/VTK output (separate library)
│   │   └── silo_export.cpp
│   ├── gts-cpt-cli/         # CLI application
│   │   └── main.cpp
│   └── gts_tools.cpp       # Independent utilities
├── include/
│   └── (deprecated - use modules)
├── tests/
│   ├── unit/
│   └── integration/
└── CMakeLists.txt          # Modern CMake with optional components
```

### Dependencies After Refactoring
```
Required:
  - GTS library (triangulated surfaces)
  - GLib 2.0 (GTS dependency)

Optional:
  - Silo (for VisIt export)
  
Removed:
  - array3d.h (dead code)
```

---

## Phase 1: Foundation (Remove Dead Code + Prepare)

### Task 1.1: Remove Dead Code

**Files:**
- Delete: `include/array3d.h`

**Step 1: Verify array3d.h is unused**

```bash
grep -r "array3d" src/ include/ tests/
# Expected: No matches (array3d.h is legacy Blitz++ template, never used)
```

**Step 2: Remove the file**

```bash
git rm include/array3d.h
git commit -m "refactor: remove unused array3d.h legacy template"
```

### Task 1.2: Audit Include Dependencies

**Files:**
- Analyze: All source files

**Step 1: Map current dependencies**

Create `docs/dependency-map.md`:
```markdown
| File | External Dependencies | Internal Dependencies |
|------|----------------------|---------------------|
| main.cpp | GTS, Silo, getopt | gtscpt.h, cpt.h, export.h |
| cpt.cpp | GTS, math.h | gtscpt.h, cpt.h |
| export.cpp | GTS, Silo | gtscpt.h, cpt.h, export.h |
| gtstools.cpp | GTS | gtstools.h |
```

**Step 2: Identify minimal core files**

Core algorithm (GTS-only):
- `cpt.cpp` / `cpt.h`
- `gtscpt.h` (SignedDistance class)

Optional export (Silo):
- `export.cpp` / `export.h`

**Step 3: Verify no Blitz++ or array3d dependencies remain**

```bash
grep -rn "array3d\|blitz" src/ include/ 
# Expected: No matches
```

### Task 1.3: Update Build System for Component Separation

**Files:**
- Create: `CMakeLists.txt` (replace Makefile)
- Create: `cmake/FindGTS.cmake`

**Step 1: Create CMakeLists.txt**

```cmake
cmake_minimum_required(VERSION 3.20)
project(gts-cpt VERSION 2.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Enable warnings as errors
add_compile_options(-Wall -Wextra -Werror -pedantic)

# Find required dependencies
find_package(PkgConfig REQUIRED)
pkg_check_modules(GTS REQUIRED gts)
pkg_check_modules(GLIB REQUIRED glib-2.0)

# Optional: Silo library
option(GTS_CPT_ENABLE_SILO "Enable Silo export support" ON)
if(GTS_CPT_ENABLE_SILO)
  find_package(SILO)
  if(NOT SILO_FOUND)
    message(WARNING "Silo not found, export functionality disabled")
    set(GTS_CPT_ENABLE_SILO OFF)
  endif()
endif()

# Core library (GTS-only dependency)
add_library(gts-cpt-core STATIC
  src/cpt.cpp
  src/gtstools.cpp
)
target_include_directories(gts-cpt-core PUBLIC
  ${GTS_INCLUDE_DIRS}
  ${GLIB_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(gts-cpt-core PUBLIC ${GTS_LIBRARIES} ${GLIB_LIBRARIES})

# Optional export library
if(GTS_CPT_ENABLE_SILO)
  add_library(gts-cpt-export STATIC
    src/export.cpp
  )
  target_include_directories(gts-cpt-export PUBLIC ${SILO_INCLUDE_DIRS})
  target_link_libraries(gts-cpt-export PUBLIC gts-cpt-core ${SILO_LIBRARIES})
endif()

# CLI application
add_executable(gts-cpt src/main.cpp)
if(GTS_CPT_ENABLE_SILO)
  target_link_libraries(gts-cpt PRIVATE gts-cpt-core gts-cpt-export)
else()
  target_link_libraries(gts-cpt PRIVATE gts-cpt-core)
endif()

# Tests
option(BUILD_TESTING "Build tests" ON)
if(BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests)
endif()
```

**Step 2: Create tests/CMakeLists.txt**

```cmake
add_executable(test_sign test_sign.cpp)
add_executable(test_safe_mesh_size test_safe_mesh_size.cpp)
add_executable(test_cpt_math test_cpt_math.cpp)

foreach(test_name IN ITEMS test_sign test_safe_mesh_size test_cpt_math)
  target_include_directories(${test_name} PRIVATE ${CMAKE_SOURCE_DIR}/include)
  target_link_libraries(${test_name} PRIVATE m)
  add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()
```

**Step 3: Test build**

```bash
mkdir build && cd build
cmake .. -DGTS_CPT_ENABLE_SILO=ON
make -j$(nproc)
```

**Step 4: Verify tests**

```bash
cd build && ctest --output-on-failure
# Expected: All tests pass
```

**Step 5: Keep old Makefile for reference**

```bash
git mv Makefile Makefile.legacy
git add CMakeLists.txt cmake/ tests/CMakeLists.txt
git commit -m "build: add CMake build system with component separation

- Create gts-cpt-core library (GTS-only dependency)
- Create gts-cpt-export library (optional Silo dependency)
- Modern CMake 3.20 with C++20
- Preserve old Makefile as Makefile.legacy"
```

---

## Phase 2: Core Library Refactoring

### Task 2.1: Extract SignedDistance Class

**Files:**
- Create: `src/gts-cpt/signed_distance.hpp`
- Create: `src/gts-cpt/signed_distance.cpp`
- Modify: `include/gtscpt.h`

**Step 1: Create signed_distance.hpp**

```cpp
// src/gts-cpt/signed_distance.hpp
#pragma once

#include <vector>
#include <memory>
#include <limits>
#include <stdexcept>

#include <gts.h>

namespace gts_cpt {

class SignedDistance {
public:
    std::vector<gdouble> value;
    gdouble distCut = 0.0;
    gdouble distMax = 0.0;
    gdouble sigma = 0.0;
    GtsSurface* pSurface = nullptr;
    
    GtsVector coordMin = {std::numeric_limits<gdouble>::max(),
                          std::numeric_limits<gdouble>::max(),
                          std::numeric_limits<gdouble>::max()};
    GtsVector coordMax = {std::numeric_limits<gdouble>::lowest(),
                          std::numeric_limits<gdouble>::lowest(),
                          std::numeric_limits<gdouble>::lowest()};
    GtsVector delta = {0.0, 0.0, 0.0};
    size_t size[3] = {0, 0, 0};
    
    SignedDistance() = default;
    ~SignedDistance();
    
    // Delete copy
    SignedDistance(const SignedDistance&) = delete;
    SignedDistance& operator=(const SignedDistance&) = delete;
    
    // Move operations
    SignedDistance(SignedDistance&& other) noexcept;
    SignedDistance& operator=(SignedDistance&& other) noexcept;
    
    // Safe mesh size calculation
    void allocate_mesh(size_t sx, size_t sy, size_t sz);
    
private:
    static size_t safe_mesh_size(size_t sx, size_t sy, size_t sz);
};

// Utility namespace
namespace util {
    [[nodiscard]] constexpr double sign(double x) noexcept {
        return std::copysign(1.0, x);
    }
}

} // namespace gts_cpt
```

**Step 2: Create signed_distance.cpp**

```cpp
// src/gts-cpt/signed_distance.cpp
#include "signed_distance.hpp"

namespace gts_cpt {

SignedDistance::~SignedDistance() {
    if (pSurface) {
        gts_object_destroy(GTS_OBJECT(pSurface));
        pSurface = nullptr;
    }
}

SignedDistance::SignedDistance(SignedDistance&& other) noexcept
    : value(std::move(other.value))
    , distCut(other.distCut)
    , distMax(other.distMax)
    , sigma(other.sigma)
    , pSurface(other.pSurface)
{
    for (int i = 0; i < 3; ++i) {
        coordMin[i] = other.coordMin[i];
        coordMax[i] = other.coordMax[i];
        delta[i] = other.delta[i];
        size[i] = other.size[i];
    }
    other.pSurface = nullptr;
}

SignedDistance& SignedDistance::operator=(SignedDistance&& other) noexcept {
    if (this != &other) {
        value = std::move(other.value);
        distCut = other.distCut;
        distMax = other.distMax;
        sigma = other.sigma;
        if (pSurface) gts_object_destroy(GTS_OBJECT(pSurface));
        pSurface = other.pSurface;
        other.pSurface = nullptr;
        for (int i = 0; i < 3; ++i) {
            coordMin[i] = other.coordMin[i];
            coordMax[i] = other.coordMax[i];
            delta[i] = other.delta[i];
            size[i] = other.size[i];
        }
    }
    return *this;
}

size_t SignedDistance::safe_mesh_size(size_t sx, size_t sy, size_t sz) {
    if (sx == 0 || sy == 0 || sz == 0) return 0;
    constexpr size_t max_size = std::numeric_limits<size_t>::max();
    if (sx > max_size / sy) {
        throw std::overflow_error("Mesh size overflow: X * Y exceeds limits");
    }
    size_t sxy = sx * sy;
    if (sxy > max_size / sz) {
        throw std::overflow_error("Mesh size overflow: total size exceeds limits");
    }
    return sxy * sz;
}

void SignedDistance::allocate_mesh(size_t sx, size_t sy, size_t sz) {
    size_t mesh_size = safe_mesh_size(sx, sy, sz);
    value.resize(mesh_size);
    size[0] = sx;
    size[1] = sy;
    size[2] = sz;
}

} // namespace gts_cpt
```

**Step 3: Update gtscpt.h to include new header**

```cpp
// include/gtscpt.h
#pragma once

// Forward declarations for C library
extern "C" {
#include <gts.h>
}

// Include refactored C++ component
#include "signed_distance.hpp"

// Backward compatibility
using SignedDistance = gts_cpt::SignedDistance;
using SizeVector = size_t[3];
```

**Step 4: Verify compilation**

```bash
cd build && make -j$(nproc)
# Expected: Compiles successfully
```

**Step 5: Run tests**

```bash
cd build && ctest
# Expected: All tests pass
```

**Step 6: Commit**

```bash
git add src/gts-cpt/ include/gtscpt.h
git commit -m "refactor: extract SignedDistance class into separate module

- Create signed_distance.hpp/cpp with proper C++20 namespace
- Keep backward compatibility with old header
- Safe mesh allocation with overflow checks"
```

### Task 2.2: Extract Math Utilities

**Files:**
- Create: `src/gts-cpt/math_utils.hpp`
- Create: `src/gts-cpt/math_utils.cpp`
- Modify: `src/cpt.cpp`, `include/cpt.h`

**Step 1: Create math_utils.hpp**

```cpp
// src/gts-cpt/math_utils.hpp
#pragma once

#include <cmath>
#include <numbers>

namespace gts_cpt {
namespace math {

// Safe sign function - replaces dangerous SIGN macro
[[nodiscard]] constexpr double sign(double x) noexcept {
    return std::copysign(1.0, x);
}

// Safe mesh size with overflow protection
[[nodiscard]] size_t safe_mesh_size(size_t sx, size_t sy, size_t sz);

// Utility functions (constexpr where possible)
[[nodiscard]] constexpr double min(double a, double b) noexcept {
    return (a < b) ? a : b;
}

[[nodiscard]] constexpr double max(double a, double b) noexcept {
    return (a > b) ? a : b;
}

[[nodiscard]] double vector_angle(const double v1[3], const double v2[3]) noexcept;

[[nodiscard]] constexpr double get_dist_cut(double size_sup, double delta_max) noexcept {
    constexpr double SQRT_3 = std::numbers::sqrt3;
    return SQRT_3 * 0.5 * size_sup * delta_max;
}

[[nodiscard]] constexpr double get_dist_max(double dist_cut, double dist_extra) noexcept {
    return dist_cut + dist_extra;
}

} // namespace math
} // namespace gts_cpt
```

**Step 2: Create math_utils.cpp**

```cpp
// src/gts-cpt/math_utils.cpp
#include "math_utils.hpp"
#include <stdexcept>
#include <limits>

namespace gts_cpt {
namespace math {

size_t safe_mesh_size(size_t sx, size_t sy, size_t sz) {
    if (sx == 0 || sy == 0 || sz == 0) return 0;
    
    constexpr size_t max_size = std::numeric_limits<size_t>::max();
    
    if (sx > max_size / sy) {
        throw std::overflow_error("Mesh size overflow: X * Y exceeds limits");
    }
    size_t sxy = sx * sy;
    
    if (sxy > max_size / sz) {
        throw std::overflow_error("Mesh size overflow: total size exceeds limits");
    }
    
    return sxy * sz;
}

double vector_angle(const double v1[3], const double v2[3]) noexcept {
    const double pvx = v1[1]*v2[2] - v1[2]*v2[1];
    const double pvy = v1[2]*v2[0] - v1[0]*v2[2];
    const double pvz = v1[0]*v2[1] - v1[1]*v2[0];
    
    return std::atan2(
        std::sqrt(pvx*pvx + pvy*pvy + pvz*pvz),
        v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
    );
}

} // namespace math
} // namespace gts_cpt
```

**Step 3: Update cpt.h for compatibility**

```cpp
// include/cpt.h - backward compatibility layer
#pragma once

#include "math_utils.hpp"

// Backward compatibility - delegate to math namespace
inline gdouble cpt_min(gdouble a, gdouble b) noexcept {
    return gts_cpt::math::min(a, b);
}

inline gdouble cpt_max(gdouble a, gdouble b) noexcept {
    return gts_cpt::math::max(a, b);
}

inline gdouble cpt_vector_angle(const GtsVector v1, const GtsVector v2) noexcept {
    return gts_cpt::math::vector_angle(v1, v2);
}

// ... rest of CPT function declarations
```

**Step 4: Verify and commit**

```bash
cd build && make && ctest
git add src/gts-cpt/math_utils.*
git commit -m "refactor: extract math utilities into separate module

- Create math_utils.hpp with constexpr functions
- Use std::numbers::sqrt3 (C++20)
- Maintain backward compatibility layer"
```

### Task 2.3: Add C++20 Concepts for Type Safety

**Files:**
- Create: `src/gts-cpt/concepts.hpp`

**Step 1: Create concepts.hpp**

```cpp
// src/gts-cpt/concepts.hpp
#pragma once

#include <concepts>
#include <type_traits>

namespace gts_cpt {

// Concept for GTS compatible floating point
template<typename T>
concept GTSFloat = std::floating_point<T> && (sizeof(T) == 4 || sizeof(T) == 8);

// Concept for mesh size type
template<typename T>
concept MeshSize = std::integral<T> && std::is_unsigned_v<T>;

// Concept for callable distance function
template<typename F>
concept DistanceFunction = requires(F f, void* point) {
    { f(point) } -> std::convertible_to<double>;
};

} // namespace gts_cpt
```

**Step 2: Use concepts in signed_distance.hpp**

```cpp
// Add to signed_distance.hpp
#include "concepts.hpp"

namespace gts_cpt {
    template<GTSFloat T = gdouble>
    // ... rest of class
}
```

**Step 3: Commit**

```bash
git add src/gts-cpt/concepts.hpp
git commit -m "feat: add C++20 concepts for type safety

- GTSFloat concept for GTS-compatible floating types
- MeshSize concept for unsigned integral sizes
- DistanceFunction concept for callable validation"
```

---

## Phase 3: Export Module Separation

### Task 3.1: Create Optional Export Library

**Files:**
- Create: `src/gts-cpt-export/export.cppm`
- Move: `src/export.cpp` → `src/gts-cpt-export/silo_export.cpp`

**Step 1: Create export module interface**

```cpp
// src/gts-cpt-export/export.cppm
export module gts_cpt.export;

export {
#include "silo_export.hpp"
}
```

**Step 2: Create silo_export.hpp**

```cpp
// src/gts-cpt-export/silo_export.hpp
#pragma once

#include "gts-cpt/signed_distance.hpp"
#include <cstddef>

#ifdef GTS_CPT_HAS_SILO
#include <silo.h>

namespace gts_cpt::export {

// Write to Silo format for VisIt visualization
bool write_silo(const SignedDistance& sd, const char* filename);

// Write to Silo with zone-centered data
bool write_silo_zones(const SignedDistance& sd, const char* filename);

} // namespace gts_cpt::export
#else
// Stub implementation when Silo is not available
namespace gts_cpt::export {
inline bool write_silo(const SignedDistance&, const char*) {
    return false; // Silo support not compiled
}
inline bool write_silo_zones(const SignedDistance&, const char*) {
    return false;
}
} // namespace gts_cpt::export
#endif

#endif // GTS_CPT_EXPORT_HPP
```

**Step 3: Update CMakeLists.txt**

```cmake
# Optional export library
if(GTS_CPT_ENABLE_SILO)
  add_library(gts-cpt-export STATIC
    src/gts-cpt-export/silo_export.cpp
  )
  target_compile_definitions(gts-cpt-export PUBLIC GTS_CPT_HAS_SILO)
  target_include_directories(gts-cpt-export PUBLIC ${SILO_INCLUDE_DIRS})
  target_link_libraries(gts-cpt-export PUBLIC gts-cpt-core ${SILO_LIBRARIES})
endif()
```

**Step 4: Update main.cpp to conditionally use export**

```cpp
// src/main.cpp
#include "gts-cpt/signed_distance.hpp"

#ifdef GTS_CPT_HAS_SILO
#include "gts-cpt-export/silo_export.hpp"
#endif

// ... in main()
#ifdef GTS_CPT_HAS_SILO
    gts_cpt::export::write_silo(*pSignedDistance, "out/eulerianmesh.silo");
#else
    std::cerr << "Warning: Silo export not available\n";
#endif
```

**Step 5: Verify optional build**

```bash
# Build without Silo
cd build
cmake .. -DGTS_CPT_ENABLE_SILO=OFF
make -j$(nproc)
ctest

# Build with Silo
cmake .. -DGTS_CPT_ENABLE_SILO=ON
make -j$(nproc)
ctest
```

**Step 6: Commit**

```bash
git add src/gts-cpt-export/ CMakeLists.txt src/main.cpp
git commit -m "refactor: create optional export module

- Silo dependency is now optional
- Export functionality compiled conditionally
- Main application gracefully handles missing export support"
```

---

## Phase 4: CLI Refactoring

### Task 4.1: Separate CLI from Core Logic

**Files:**
- Create: `src/gts-cpt-cli/cli_parser.hpp`
- Create: `src/gts-cpt-cli/cli_parser.cpp`
- Modify: `src/main.cpp`

**Step 1: Create cli_parser.hpp**

```cpp
// src/gts-cpt-cli/cli_parser.hpp
#pragma once

#include "gts-cpt/signed_distance.hpp"
#include <string>
#include <optional>
#include <array>

namespace gts_cpt::cli {

struct CliOptions {
    // Mesh bounds
    std::array<double, 3> begin = {0.0, 0.0, 0.0};
    std::array<double, 3> end = {0.0, 0.0, 0.0};
    
    // Mesh resolution
    std::array<long, 3> size = {0, 0, 0};
    
    // Flags
    bool verbose = false;
    bool normalize = false;
    
    // Input file (optional, stdin if empty)
    std::optional<std::string> input_file;
    
    // Output directory
    std::string output_dir = "out";
    
    // Parse from command line
    static std::optional<CliOptions> parse(int argc, char* argv[]);
    
    // Validate options
    bool validate() const;
    
    // Print usage
    static void print_usage(const char* program_name);
};

// Process input file and compute signed distance
int process(const CliOptions& options);

} // namespace gts_cpt::cli
```

**Step 2: Create cli_parser.cpp**

```cpp
// src/gts-cpt-cli/cli_parser.cpp
#include "cli_parser.hpp"
#include <getopt.h>
#include <iostream>
#include <fstream>

namespace gts_cpt::cli {

namespace {
    constexpr long MIN_MESH_SIZE = 1;
}

std::optional<CliOptions> CliOptions::parse(int argc, char* argv[]) {
    CliOptions opts;
    
    static struct option long_options[] = {
        {"begin-x",   required_argument, nullptr, 'A'},
        {"begin-y",   required_argument, nullptr, 'B'},
        {"begin-z",   required_argument, nullptr, 'C'},
        {"end-x",     required_argument, nullptr, 'D'},
        {"end-y",     required_argument, nullptr, 'E'},
        {"end-z",     required_argument, nullptr, 'F'},
        {"size-x",    required_argument, nullptr, 'G'},
        {"size-y",    required_argument, nullptr, 'H'},
        {"size-z",    required_argument, nullptr, 'I'},
        {"normalize", no_argument,       nullptr, 'o'},
        {"verbose",   no_argument,       nullptr, 'v'},
        {"help",      no_argument,       nullptr, 'h'},
        {nullptr,     0,                 nullptr, 0}
    };
    
    int c;
    while ((c = getopt_long(argc, argv, "A:B:C:D:E:F:G:H:I:os:vh", long_options, nullptr)) != EOF) {
        switch (c) {
            case 'A': opts.begin[0] = std::stod(optarg); break;
            case 'B': opts.begin[1] = std::stod(optarg); break;
            case 'C': opts.begin[2] = std::stod(optarg); break;
            case 'D': opts.end[0] = std::stod(optarg); break;
            case 'E': opts.end[1] = std::stod(optarg); break;
            case 'F': opts.end[2] = std::stod(optarg); break;
            case 'G':
                opts.size[0] = std::stol(optarg);
                if (opts.size[0] <= 0) {
                    std::cerr << "Error: size-x must be positive\n";
                    return std::nullopt;
                }
                break;
            case 'H':
                opts.size[1] = std::stol(optarg);
                if (opts.size[1] <= 0) {
                    std::cerr << "Error: size-y must be positive\n";
                    return std::nullopt;
                }
                break;
            case 'I':
                opts.size[2] = std::stol(optarg);
                if (opts.size[2] <= 0) {
                    std::cerr << "Error: size-z must be positive\n";
                    return std::nullopt;
                }
                break;
            case 'o': opts.normalize = true; break;
            case 'v': opts.verbose = true; break;
            case 'h':
                print_usage(argv[0]);
                return std::nullopt;
            default:
                return std::nullopt;
        }
    }
    
    return opts;
}

bool CliOptions::validate() const {
    // Validate mesh bounds
    if (begin[0] >= end[0] || begin[1] >= end[1] || begin[2] >= end[2]) {
        std::cerr << "Error: Invalid mesh bounds\n";
        return false;
    }
    
    // Validate sizes
    if (size[0] <= 0 || size[1] <= 0 || size[2] <= 0) {
        std::cerr << "Error: Grid size must be positive\n";
        return false;
    }
    
    return true;
}

void CliOptions::print_usage(const char* program_name) {
    std::cout << "Usage: " << program_name << " [OPTIONS] < input.gts\n"
              << "CPT using the GTS library.\n\n"
              << "Options:\n"
              << "  --begin-x VALUE    Minimum X coord\n"
              << "  --begin-y VALUE    Minimum Y coord\n"
              << "  --begin-z VALUE    Minimum Z coord\n"
              << "  --end-x VALUE      Maximum X coord\n"
              << "  --end-y VALUE      Maximum Y coord\n"
              << "  --end-z VALUE      Maximum Z coord\n"
              << "  --size-x N         Grid cells in X\n"
              << "  --size-y N         Grid cells in Y\n"
              << "  --size-z N         Grid cells in Z\n"
              << "  --normalize        Normalize to unit cube\n"
              << "  --verbose          Print progress\n"
              << "  --help             Show this help\n";
}

} // namespace gts_cpt::cli
```

**Step 3: Refactor main.cpp to use CLI**

```cpp
// src/main.cpp (simplified)
#include "cli_parser.hpp"
#include "gts-cpt/signed_distance.hpp"

#ifdef GTS_CPT_HAS_SILO
#include "gts-cpt-export/silo_export.hpp"
#endif

#include <iostream>

int main(int argc, char* argv[]) {
    auto opts = gts_cpt::cli::CliOptions::parse(argc, argv);
    if (!opts) {
        return 1;
    }
    
    if (!opts->validate()) {
        return 1;
    }
    
    return gts_cpt::cli::process(*opts);
}
```

**Step 4: Commit**

```bash
git add src/gts-cpt-cli/
git commit -m "refactor: separate CLI parsing from core logic

- Create cli_parser module with validation
- Type-safe option handling
- Clean separation of concerns"
```

---

## Phase 5: Testing Infrastructure

### Task 5.1: Add Unit Test Framework

**Files:**
- Create: `tests/CMakeLists.txt` (already done in Phase 1)
- Add: `tests/unit/CMakeLists.txt`
- Add: Test utilities

**Step 1: Install testing framework**

Use CTest (built into CMake) with minimalist approach:

```cmake
# tests/CMakeLists.txt
enable_testing()

# Unit tests
add_subdirectory(unit)

# Integration tests (optional, require GTS files)
option(BUILD_INTEGRATION_TESTS "Build integration tests" OFF)
if(BUILD_INTEGRATION_TESTS)
  add_subdirectory(integration)
endif()
```

**Step 2: Create catch.hpp wrapper**

```cpp
// tests/catch.hpp
#pragma once

// Use C++20 features
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
```

**Step 3: Update existing tests to use Catch2**

```cpp
// tests/unit/test_sign.cpp (updated)
#include <catch2/catch_test_macros.hpp>
#include "gts-cpt/math_utils.hpp"

TEST_CASE("sign function handles edge cases", "[math]") {
    using gts_cpt::math::sign;
    
    SECTION("positive values") {
        REQUIRE(sign(5.0) == Catch::Approx(1.0));
        REQUIRE(sign(0.001) == Catch::Approx(1.0));
    }
    
    SECTION("negative values") {
        REQUIRE(sign(-5.0) == Catch::Approx(-1.0));
        REQUIRE(sign(-0.001) == Catch::Approx(-1.0));
    }
    
    SECTION("zero") {
        REQUIRE(sign(0.0) == Catch::Approx(1.0)); // copysign(1.0, 0.0) = 1.0
    }
    
    SECTION("subnormal") {
        double subnormal = std::numeric_limits<double>::denorm_min();
        REQUIRE(sign(subnormal) == Catch::Approx(1.0));
        REQUIRE(sign(-subnormal) == Catch::Approx(-1.0));
    }
}
```

**Step 4: Commit**

```bash
git add tests/
git commit -m "test: add Catch2 testing framework

- Replace raw asserts with Catch2 test macros
- Add unit test organization
- Enable integration test option"
```

---

## Phase 6: Documentation and Final Cleanup

### Task 6.1: Update README.md

**Files:**
- Update: `README.md`

**Add sections:**
- Build options (with/without Silo)
- CMake configuration
- Testing
- Module architecture

### Task 6.2: Add API Documentation

**Files:**
- Create: `docs/API.md`

**Document:**
- All public APIs in C++ namespaces
- Usage examples
- Migration guide from old API

### Task 6.3: Remove Deprecated Code

**Files:**
- Remove: `Makefile.legacy`
- Remove: Old style headers from `include/`
- Keep: `include/gtscpt.h` as compatibility wrapper only

---

## Verification Checklist

After completing all phases:

```bash
# Build without Silo
cmake .. -DGTS_CPT_ENABLE_SILO=OFF
make clean && make -j$(nproc)

# Run tests
ctest --output-on-failure

# Build with Silo
cmake .. -DGTS_CPT_ENABLE_SILO=ON
make clean && make -j$(nproc)

# Run tests
ctest --output-on-failure

# Build with C++20 modules (if compiler supports)
cmake .. -DGTS_CPT_ENABLE_MODULES=ON
make -j$(nproc)

# Verify no warnings
make 2>&1 | grep -i warning
# Expected: No warnings
```

---

## Rollback Plan

If issues arise:
1. Keep `Makefile.legacy` for backward compatibility
2. Keep old headers in `include/` with deprecation warnings
3. All phases are reversible via git revert

---

## Summary

This refactoring plan will:
1. Remove all dead code (array3d.h)
2. Separate core (GTS-only) from optional (Silo) dependencies
3. Modernize to C++20 with concepts and modules
4. Establish clean module boundaries
5. Add comprehensive test coverage
6. Enable CMake build system with optional components

**Estimated time:** 3-5 days for experienced developer
**Risk level:** Medium (incremental phases, each independently testable)