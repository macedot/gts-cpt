# GTS-CPT C++ RAII Encapsulation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Encapsulate GTS C library calls with modern C++ RAII patterns, eliminating raw pointers for owned resources.

**Architecture:** Create C++ wrapper classes that manage GTS object lifecycles via RAII. Use `unique_ptr` for ownership, provide clean C++ interfaces that hide C implementation details.

**Tech Stack:** C++17, GTS C library, GLib, RAII patterns, smart pointers

---

## Overview

### Current State
- Raw pointers to GTS objects (`GtsSurface*`, `GtsVertex*`, etc.)
- Manual memory management with `gts_object_destroy()`
- C-style callbacks with `gpointer` void pointers
- No resource safety on exceptions

### Target State
```cpp
// Before (C-style)
GtsSurface* surface = gts_surface_new(...);
// ... use surface ...
gts_object_destroy(GTS_OBJECT(surface));  // Easy to forget

// After (C++ RAII)
auto surface = gts::make_surface(...);
surface->add_triangle(v1, v2, v3);
// Automatic cleanup via destructor
```

---

## Phase 1: GTS Resource Wrappers

### Task 1.1: Create `gts_wrapper.hpp` Base

**Files:**
- Create: `include/gts-cpp/gts_wrapper.hpp`

**Step 1: Analyze current GTS usage**

```bash
grep -h "gts_surface_new\|gts_object_destroy\|GtsSurface\*\|GtsVertex\*\|GtsTriangle\*" src/*.cpp include/*.h
```

Expected: Identify all GTS object types used.

**Step 2: Create base RAII wrapper**

Create `include/gts-cpp/gts_wrapper.hpp`:

```cpp
#pragma once

#include <gts.h>
#include <memory>
#include <utility>
#include <type_traits>

namespace gts {

// Custom deleter for GTS objects
struct GtsDeleter {
    void operator()(GtsObject* obj) const noexcept {
        if (obj) {
            gts_object_destroy(obj);
        }
    }
};

// Type alias for unique_ptr with GTS deleter
template<typename T>
using gts_unique_ptr = std::unique_ptr<T, GtsDeleter>;

// Base deleter for containers
struct GtsContainerDeleter {
    void operator()(GtsContainer* container) const noexcept {
        if (container) {
            gts_object_destroy(GTS_OBJECT(container));
        }
    }
};

// Surface deleter
struct SurfaceDeleter {
    void operator()(GtsSurface* surface) const noexcept {
        if (surface) {
            gts_object_destroy(GTS_OBJECT(surface));
        }
    }
};

// Safe GTS pointer
template<typename T>
using GtsPtr = std::unique_ptr<T, std::add_pointer_t<void(T*)>>;

// Create unique_ptr with GTS deleter
template<typename T>
auto make_gts_ptr(T* ptr) {
    return std::unique_ptr<T, GtsDeleter>(ptr);
}

} // namespace gts
```

**Step 3: Verify compilation**

```bash
g++ -std=c++17 -I./include -c include/gts-cpp/gts_wrapper.hpp -o /dev/null
```

Expected: No errors (header-only).

**Step 4: Commit**

```bash
git add include/gts-cpp/
git commit -m "feat: add GTS RAII wrapper base

- Custom deleter for GTS objects
- unique_ptr alias with automatic cleanup
- Foundation for type-safe GTS wrappers"
```

### Task 1.2: Create Surface Wrapper Class

**Files:**
- Create: `include/gts-cpp/surface.hpp`
- Create: `src/gts-cpp/surface.cpp`
- Modify: `Makefile`

**Step 1: Create Surface class header**

Create `include/gts-cpp/surface.hpp`:

```cpp
#pragma once

#include "gts_wrapper.hpp"
#include <gts.h>
#include <memory>
#include <stdexcept>

namespace gts {

class Surface {
public:
    // Factory function - creates empty surface
    static std::unique_ptr<Surface> create();
    
    // Factory function - reads from file
    static std::unique_ptr<Surface> from_file(const char* filename);
    
    // Non-copyable
    Surface(const Surface&) = delete;
    Surface& operator=(const Surface&) = delete;
    
    // Movable
    Surface(Surface&& other) noexcept;
    Surface& operator=(Surface&& other) noexcept;
    
    // Destructor
    ~Surface();
    
    // Access underlying GTS surface (for C API calls)
    GtsSurface* get() noexcept { return surface_.get(); }
    const GtsSurface* get() const noexcept { return surface_.get(); }
    
    // Type-safe operations
    void add_triangle(GtsVertex* v1, GtsVertex* v2, GtsVertex* v3);
    void read(const char* filename);
    void write_vtk(const char* filename) const;
    
    // Surface properties
    bool is_closed() const;
    gdouble volume() const;
    gdouble area() const;
    
private:
    Surface() = default;
    gts_unique_ptr<GtsSurface> surface_;
};

// RAII helper for reading surface
std::unique_ptr<Surface> read_surface(const char* filename);

// RAII helper for writing surface  
void write_surface_vtk(const Surface& surface, const char* filename);

} // namespace gts
```

**Step 2: Implement Surface class**

Create `src/gts-cpp/surface.cpp`:

```cpp
#include "gts-cpp/surface.hpp"
#include <fstream>

namespace gts {

std::unique_ptr<Surface> Surface::create() {
    auto surf = std::unique_ptr<Surface>(new Surface());
    surf->surface_ = make_gts_ptr(
        gts_surface_new(
            gts_surface_class(),
            gts_face_class(),
            gts_edge_class(),
            gts_vertex_class()
        )
    );
    if (!surf->surface_) {
        throw std::runtime_error("Failed to create GTS surface");
    }
    return surf;
}

std::unique_ptr<Surface> Surface::from_file(const char* filename) {
    auto surf = create();
    surf->read(filename);
    return surf;
}

Surface::Surface(Surface&& other) noexcept
    : surface_(std::move(other.surface_))
{
}

Surface& Surface::operator=(Surface&& other) noexcept {
    if (this != &other) {
        surface_ = std::move(other.surface_);
    }
    return *this;
}

Surface::~Surface() = default;

void Surface::add_triangle(GtsVertex* v1, GtsVertex* v2, GtsVertex* v3) {
    // Create edge and face
    GtsEdge* e1 = GTS_EDGE(gts_edge_new(gts_edge_class(), v1, v2));
    GtsEdge* e2 = GTS_EDGE(gts_edge_new(gts_edge_class(), v2, v3));
    GtsEdge* e3 = GTS_EDGE(gts_edge_new(gts_edge_class(), v3, v1));
    GtsFace* face = GTS_FACE(gts_face_new(gts_face_class(), e1, e2, e3));
    gts_surface_add_face(surface_.get(), face);
}

void Surface::read(const char* filename) {
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        throw std::runtime_error("Failed to open file: " + std::string(filename));
    }
    
    GtsFile* gf = gts_file_new(fp);
    if (gts_surface_read(surface_.get(), gf) != 0) {
        fclose(fp);
        gts_file_destroy(gf);
        throw std::runtime_error("Failed to read surface from file");
    }
    
    fclose(fp);
    gts_file_destroy(gf);
}

void Surface::write_vtk(const char* filename) const {
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        throw std::runtime_error("Failed to open file for writing: " + std::string(filename));
    }
    // VTK writing implementation
    // ... 
}

bool Surface::is_closed() const {
    return gts_surface_is_closed(surface_.get()) != FALSE;
}

gdouble Surface::volume() const {
    return gts_surface_volume(surface_.get());
}

gdouble Surface::area() const {
    return gts_surface_area(surface_.get());
}

// Helper functions
std::unique_ptr<Surface> read_surface(const char* filename) {
    return Surface::from_file(filename);
}

void write_surface_vtk(const Surface& surface, const char* filename) {
    surface.write_vtk(filename);
}

} // namespace gts
```

**Step 3: Update Makefile**

Add to `Makefile`:
```makefile
SRC += src/gts-cpp/surface.cpp
INCLPATH += -I./include/gts-cpp
```

**Step 4: Verify build**

```bash
make clean && make
```

Expected: Compiles successfully.

**Step 5: Commit**

```bash
git add include/gts-cpp/surface.* src/gts-cpp/
git commit -m "feat: add Surface RAII wrapper

- Surface class with factory methods
- Automatic cleanup via destructor
- Type-safe interface for GTS surface operations
- Exception-safe resource management"
```

### Task 1.3: Create Vertex and Edge Wrappers

**Files:**
- Create: `include/gts-cpp/vertex.hpp`
- Create: `include/gts-cpp/edge.hpp`

**Step 1: Create Vertex wrapper**

Create `include/gts-cpp/vertex.hpp`:

```cpp
#pragma once

#include "gts_wrapper.hpp"
#include <gts.h>
#include <memory>

namespace gts {

class Vertex {
public:
    // Create vertex at coordinates
    static std::unique_ptr<Vertex> create(gdouble x, gdouble y, gdouble z);
    
    // Non-copyable
    Vertex(const Vertex&) = delete;
    Vertex& operator=(const Vertex&) = delete;
    
    // Access
    GtsVertex* get() noexcept { return vertex_.get(); }
    const GtsVertex* get() const noexcept { return vertex_.get(); }
    
    // Coordinates
    gdouble x() const noexcept;
    gdouble y() const noexcept;
    gdouble z() const noexcept;
    
private:
    Vertex() = default;
    gts_unique_ptr<GtsVertex> vertex_;
};

} // namespace gts
```

**Step 2: Create Edge wrapper (similar pattern)**

Similar structure for Edge class.

**Step 3: Verify and commit**

```bash
make && git add -A && git commit -m "feat: add Vertex and Edge RAII wrappers"
```

---

## Phase 2: SignedDistance Class Refactoring

### Task 2.1: Update SignedDistance to Use RAII Wrappers

**Files:**
- Modify: `include/gtscpt.h`
- Modify: `src/main.cpp`

**Step 1: Update SignedDistance to use gts::Surface**

Modify `include/gtscpt.h`:

```cpp
#include "gts-cpp/surface.hpp"

class SignedDistance {
public:
    std::vector<gdouble> value;
    gdouble distCut = 0.0;
    gdouble distMax = 0.0;
    gdouble sigma = 0.0;
    
    // Use gts::Surface wrapper instead of raw pointer
    std::unique_ptr<gts::Surface> surface;
    
    GtsVector coordMin = {std::numeric_limits<gdouble>::max(), ...};
    GtsVector coordMax = {std::numeric_limits<gdouble>::lowest(), ...};
    GtsVector delta = {0.0, 0.0, 0.0};
    SizeVector size = {0, 0, 0};
    
    // ... rest same
};
```

**Step 2: Update main.cpp usage**

```cpp
// Before
pSignedDistance->pSurface = gts_surface_new(...);

// After
pSignedDistance->surface = gts::Surface::create();
```

**Step 3: Remove manual cleanup**

```cpp
// Remove these lines:
gts_object_destroy(GTS_OBJECT(pSignedDistance->pSurface));
```

**Step 4: Verify all tests pass**

```bash
cd tests && make check
```

**Step 5: Commit**

```bash
git add -A && git commit -m "refactor: use gts::Surface in SignedDistance

- Replace raw GtsSurface* with gts::Surface wrapper
- Automatic cleanup via RAII
- Remove manual gts_object_destroy calls"
```

---

## Phase 3: Exception Safety

### Task 3.1: Add Exception-Safe Surface Loading

**Files:**
- Modify: `src/main.cpp`

**Step 1: Wrap surface loading in try-catch**

```cpp
try {
    pSignedDistance->surface = gts::Surface::from_file(input_file);
} catch (const std::exception& e) {
    std::cerr << "Error loading surface: " << e.what() << std::endl;
    return 1;
}
```

**Step 2: Ensure cleanup on early returns**

With RAII, all cleanup happens automatically. Verify:

```cpp
// No more manual cleanup needed - RAII handles it
try {
    // ... all operations ...
} catch (...) {
    // unique_ptr automatically cleans up
    throw;
}
```

**Step 3: Test exception paths**

Verify that exceptions don't leak resources.

---

## Phase 4: C++ Style Callbacks

### Task 4.1: Create C++ Callback Wrapper

**Files:**
- Create: `include/gts-cpp/callbacks.hpp`

**Step 1: Create type-safe callback wrapper**

```cpp
#pragma once

#include <gts.h>
#include <functional>
#include <type_traits>

namespace gts {

// Type-erased callback data
template<typename Func>
struct CallbackData {
    Func func;
    
    template<typename F>
    CallbackData(F&& f) : func(std::forward<F>(f)) {}
};

// Helper to call C++ callback from C
template<typename Func>
void callback_adapter(GtsObject* obj, gpointer data) {
    auto* cb_data = static_cast<CallbackData<Func>*>(data);
    cb_data->func(obj);
}

// RAII wrapper for callback registration
template<typename Func>
class CallbackScope {
    CallbackData<Func> data_;
    
public:
    template<typename F>
    CallbackScope(F&& func) : data_(std::forward<F>(func)) {}
    
    gpointer data() noexcept { return &data_; }
    gpointer data() const noexcept { return const_cast<gpointer>(&data_); }
};

} // namespace gts
```

**Step 2: Use in CPT functions**

```cpp
// Before (C-style)
void cpt_vertex(GtsVertex* v, gpointer data);

// After (C++ wrapper)
gts::CallbackScope<std::function<void(GtsVertex*)>> callback([](GtsVertex* v) {
    // C++ code here
});

gts_surface_foreach_vertex(surface, gts::callback_adapter, callback.data());
```

---

## Phase 5: Integration

### Task 5.1: Update CPT Functions to Use Wrappers

**Files:**
- Modify: `src/cpt.cpp`

**Step 1: Update function signatures**

```cpp
// Before
void cpt_face(GtsTriangle *pTriangle, GtsVector normal, SignedDistance* pSignedDistance, gdouble inside_outside);

// After - use gts::Surface&
void cpt_face(const gts::Surface& surface, /* ... */);
```

**Step 2: Update internal usage**

Replace raw pointer access with wrapper methods where applicable.

---

## Verification

After all phases:

```bash
# Build
make clean && make

# Run tests
cd tests && make check

# Check for memory leaks (if valgrind available)
valgrind --leak-check=full ./debug/gts-cpt --verbose ...

# Verify no manual deletes/destroys
grep -rn "delete\|destroy\|free" src/*.cpp include/*.h
# Expected: Only in RAII destructors, not manual calls
```

---

## Summary

This plan transforms the codebase from C-style manual resource management to modern C++ RAII:

1. **GTS Wrapper Base** - Custom deleters and unique_ptr aliases
2. **Surface/Vertex/Edge Wrappers** - Type-safe RAII classes
3. **SignedDistance Update** - Use gts::Surface instead of raw pointer
4. **Exception Safety** - Automatic cleanup on exceptions
5. **C++ Callbacks** - Type-safe callback wrappers

**Key Benefits:**
- No more raw pointers for owned resources
- Automatic cleanup via destructors
- Exception-safe resource management
- Cleaner API hiding C implementation

**Estimated Time:** 1-2 sessions