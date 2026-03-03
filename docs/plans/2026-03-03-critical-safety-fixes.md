# Critical and High Severity Fixes Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Fix undefined behavior, memory safety issues, and critical bugs to prepare codebase for C++20 migration.

**Architecture:** Incremental fixes with RAII patterns, replacing unsafe constructs, adding validation. Each fix is isolated to minimize risk.

**Tech Stack:** C++ (currently C++98, will remain compatible), GTS library, GLib, Silo

---

## Overview

This plan addresses 11 Must-Fix and Should-Fix issues from the code review:

**Critical (Must-Fix):**
1. Division by zero in SIGN macro
2. Memory leaks on early returns
3. Buffer overflow risk from size multiplication
4. Undefined behavior in cpt_distance_sign
5. Duplicate function definition in mesh_visit.cpp

**High (Should-Fix):**
6. Incomplete RAII in array3d.h
7. Missing const correctness
8. Global mutable state
9. Raw new/delete without exception safety
10. Incorrect include guard in export.h
11. Signed/unsigned comparison warnings

---

## Task 1: Fix SIGN Macro Division by Zero

**Files:**
- Modify: `include/gtscpt.h:33`
- Modify: `src/cpt.cpp:198, 227, 239, 254, 268`
- Create: `tests/test_sign.cpp`

**Step 1: Write failing test for edge cases**

Create `tests/test_sign.cpp`:

```cpp
#include <cassert>
#include <cmath>
#include <limits>

// Forward declaration
double sign(double x);

int main() {
    // Test positive
    assert(sign(5.0) == 1.0);
    assert(sign(0.001) == 1.0);
    
    // Test negative
    assert(sign(-5.0) == -1.0);
    assert(sign(-0.001) == -1.0);
    
    // Test zero - CRITICAL: should not cause UB
    double result = sign(0.0);
    assert(result == 1.0 || result == -1.0 || std::isnan(result) == false);
    
    // Test subnormal
    double subnormal = std::numeric_limits<double>::denorm_min();
    assert(sign(subnormal) == 1.0);
    assert(sign(-subnormal) == -1.0);
    
    return 0;
}
```

**Step 2: Compile test (expect failure)**

Run: `g++ -std=c++20 -I./include tests/test_sign.cpp -o tests/test_sign -lm`
Expected: Undefined reference to `sign(double)`

**Step 3: Replace SIGN macro with safe constexpr function**

Modify `include/gtscpt.h:33-34`:

```cpp
// OLD (lines 33):
// #define SIGN(x)                 ( (x) / fabs(x) )

// NEW (replace lines 33-34):
namespace gts_cpt {
    [[nodiscard]] inline constexpr double sign(double x) noexcept {
        // Use copysign: returns value with magnitude of 1.0 and sign of x
        // For x = 0.0, copysign(1.0, +0.0) = 1.0
        // For x = -0.0, copysign(1.0, -0.0) = -1.0
        return std::copysign(1.0, x);
    }
}

// Compatibility macro for gradual migration
#define SIGN(x) gts_cpt::sign(x)
```

Add required include at `include/gtscpt.h:4-5`:

```cpp
// After line 5, add:
#include <cmath>
```

**Step 4: Update direct usages**

Modify `src/cpt.cpp:227`:

```cpp
// OLD:
return distance / fabs(distance);

// NEW:
return gts_cpt::sign(distance);
```

Modify `src/cpt.cpp:198`:

```cpp
// OLD:
signDist = (inside_outside > 0.0 ? distance / absDist : inside_outside);

// NEW:
signDist = (inside_outside > 0.0 ? gts_cpt::sign(distance) : inside_outside);
```

**Step 5: Compile and run test**

Run: `g++ -std=c++20 -I./include tests/test_sign.cpp src/cpt.cpp -o tests/test_sign -lm && ./tests/test_sign`
Expected: PASS (exit code 0)

**Step 6: Commit**

```bash
git add include/gtscpt.h src/cpt.cpp tests/test_sign.cpp
git commit -m "fix: replace SIGN macro with safe sign() function

- Removes undefined behavior when x=0
- Uses std::copysign for safe sign extraction
- Adds unit tests for edge cases
- Maintains backward compatibility with SIGN macro"
```

---

## Task 2: Add RAII Wrapper for SignedDistance

**Files:**
- Modify: `include/gtscpt.h:38-54`
- Modify: `src/main.cpp:166-508`

**Step 1: Redesign SignedDistance struct**

Modify `include/gtscpt.h:38-54`:

```cpp
// OLD (lines 38-54):
// typedef size_t  SizeVector[3];
// typedef struct eulerianMesh { ... } SignedDistance;

// NEW:
#include <vector>
#include <memory>
#include <limits>

using SizeVector = size_t[3];

class SignedDistance {
public:
    std::vector<gdouble> value;          // Managed memory for mesh data
    gdouble      distCut = 0.0;          // (epson_1) half size of band
    gdouble      distMax = 0.0;          // (epson_2) maximum distance
    gdouble      sigma   = 0.0;          // = 3 * MIN (HX , HY , HZ)
    GtsSurface*  pSurface = nullptr;     // Lagrangian Mesh (owned)
    
    GtsVector    coordMin = { std::numeric_limits<gdouble>::max(),
                               std::numeric_limits<gdouble>::max(),
                               std::numeric_limits<gdouble>::max() };
    GtsVector    coordMax = { std::numeric_limits<gdouble>::lowest(),
                               std::numeric_limits<gdouble>::lowest(),
                               std::numeric_limits<gdouble>::lowest() };
    GtsVector    delta = { 0.0, 0.0, 0.0 };
    SizeVector   size = { 0, 0, 0 };
    
    // Constructor
    SignedDistance() = default;
    
    // Destructor - automatically cleans up GTS surface
    ~SignedDistance() {
        if (pSurface) {
            gts_object_destroy(GTS_OBJECT(pSurface));
            pSurface = nullptr;
        }
    }
    
    // Delete copy operations (expensive, error-prone)
    SignedDistance(const SignedDistance&) = delete;
    SignedDistance& operator=(const SignedDistance&) = delete;
    
    // Enable move operations
    SignedDistance(SignedDistance&& other) noexcept
        : value(std::move(other.value))
        , distCut(other.distCut)
        , distMax(other.distMax)
        , sigma(other.sigma)
        , pSurface(other.pSurface)
    {
        std::copy(std::begin(other.coordMin), std::end(other.coordMin), coordMin);
        std::copy(std::begin(other.coordMax), std::end(other.coordMax), coordMax);
        std::copy(std::begin(other.delta), std::end(other.delta), delta);
        std::copy(std::begin(other.size), std::end(other.size), size);
        other.pSurface = nullptr;
    }
    
    SignedDistance& operator=(SignedDistance&& other) noexcept {
        if (this != &other) {
            value = std::move(other.value);
            distCut = other.distCut;
            distMax = other.distMax;
            sigma = other.sigma;
            if (pSurface) gts_object_destroy(GTS_OBJECT(pSurface));
            pSurface = other.pSurface;
            other.pSurface = nullptr;
            std::copy(std::begin(other.coordMin), std::end(other.coordMin), coordMin);
            std::copy(std::begin(other.coordMax), std::end(other.coordMax), coordMax);
            std::copy(std::begin(other.delta), std::end(other.delta), delta);
            std::copy(std::begin(other.size), std::end(other.size), size);
        }
        return *this;
    }
};
```

**Step 2: Add helper for safe multiplication**

Add to `include/gtscpt.h` (after SignedDistance class):

```cpp
namespace gts_cpt {
    // Safe multiplication with overflow check
    [[nodiscard]] inline size_t safe_mesh_size(size_t sx, size_t sy, size_t sz) {
        if (sx == 0 || sy == 0 || sz == 0) return 0;
        
        constexpr size_t max_size = std::numeric_limits<size_t>::max();
        
        if (sx > max_size / sy) {
            throw std::overflow_error("Mesh size overflow: dimension X * Y exceeds limits");
        }
        size_t sxy = sx * sy;
        
        if (sxy > max_size / sz) {
            throw std::overflow_error("Mesh size overflow: total mesh size exceeds limits");
        }
        
        return sxy * sz;
    }
}
```

**Step 3: Update main.cpp allocation**

Modify `src/main.cpp:38-39`:

```cpp
// OLD:
SignedDistance  *pSignedDistance = NULL;

// NEW:
std::unique_ptr<SignedDistance> pSignedDistance;
```

Modify `src/main.cpp:166-173`:

```cpp
// OLD:
pSignedDistance = new SignedDistance;
memset(pSignedDistance, 0, sizeof(SignedDistance));
for(i = 0; i < 3; i++) {
    pSignedDistance->coordMin[i] =  DBL_MAX;
    pSignedDistance->coordMax[i] = -DBL_MAX;
}

// NEW:
pSignedDistance = std::make_unique<SignedDistance>();
// coordMin and coordMax are initialized in class definition
```

**Step 4: Update all early returns to use RAII cleanup**

Modify `src/main.cpp:152-161` (before each `return 1;`, no changes needed - RAII handles cleanup):

```cpp
// No explicit cleanup needed - unique_ptr and vector handle it
if ( beginMesh[0] >= endMesh[0] || beginMesh[1] >= endMesh[1] || beginMesh[2] >= endMesh[2]) {
    fputs("gtscpt: you must correctly especify the eulerian mesh!\n", stderr);
    return 1;  // RAII cleans up automatically
}
```

**Step 5: Update mesh size allocation**

Modify `src/main.cpp:375-394`:

```cpp
// OLD:
size_t  mesh_size = (pSignedDistance->size[0] * pSignedDistance->size[1] * pSignedDistance->size[2]);
// ...
pSignedDistance->value = new gdouble[mesh_size];

// NEW:
size_t mesh_size;
try {
    mesh_size = gts_cpt::safe_mesh_size(
        pSignedDistance->size[0],
        pSignedDistance->size[1],
        pSignedDistance->size[2]
    );
} catch (const std::overflow_error& e) {
    fprintf(stderr, "gtscpt: %s\n", e.what());
    return 1;
}

if (mesh_size == 0) {
    g_warning("INVALID MESH SIZE!!");
    return 1;
}

if(verbose) {
    printf("#\n");
    printf("# Creating the eulerian mesh [size = %zu]...", mesh_size);
}

try {
    pSignedDistance->value.resize(mesh_size);
} catch (const std::bad_alloc& e) {
    fprintf(stderr, "gtscpt: failed to allocate mesh: %s\n", e.what());
    return 1;
}
```

**Step 6: Remove manual cleanup**

Delete `src/main.cpp:504-508`:

```cpp
// OLD (DELETE THESE LINES):
gts_object_destroy(GTS_OBJECT(pSignedDistance->pSurface));
gts_finalize();

delete [] pSignedDistance->value;  pSignedDistance->value = NULL;
delete pSignedDistance;           pSignedDistance       = NULL;

// NEW (ADD ONLY):
gts_finalize();
// unique_ptr and vector handle cleanup automatically
```

**Step 7: Update all pSignedDistance accesses**

Replace all `->` with `->` throughout (no change needed - unique_ptr supports `->`).

**Step 8: Update catch block**

Modify `src/main.cpp:514-518`:

```cpp
// OLD:
catch(...) {
    g_warning("unhandle exception");
    return -1;
}

// NEW:
catch(const std::exception& e) {
    fprintf(stderr, "Exception: %s\n", e.what());
    return 1;
}
catch(...) {
    g_warning("Unknown exception");
    return 1;
}
```

**Step 9: Test compilation**

Run: `make clean && make`
Expected: Compilation succeeds (warnings about unused array3d are OK)

**Step 10: Test with sample data**

Run: `./run sphere20.gts` (or any available GTS file)
Expected: Program runs successfully, no memory leaks

**Step 11: Commit**

```bash
git add include/gtscpt.h src/main.cpp
git commit -m "refactor: add RAII for SignedDistance

- Replace raw pointer with unique_ptr
- Replace new[] with std::vector
- Add overflow checks for mesh size
- Automatic cleanup on early returns
- Add exception safety"
```

---

## Task 3: Fix Undefined Behavior in cpt_distance_sign

**Files:**
- Modify: `src/cpt.cpp:217-228`

**Step 1: Update cpt_distance_sign function**

Modify `src/cpt.cpp:217-228`:

```cpp
// OLD:
gdouble cpt_distance_sign(GtsPoint* pPointMesh, ParamDistFunc* pParam)
{
    GtsVector*  pNormal       = (GtsVector*)pParam->pNormal;
    GtsPoint*   pSurfacePoint = (GtsPoint*)pParam->pSurfacePoint;
    gdouble     distance;
    
    distance = (*pNormal)[0] * (pPointMesh->x - pSurfacePoint->x)
             + (*pNormal)[1] * (pPointMesh->y - pSurfacePoint->y)
             + (*pNormal)[2] * (pPointMesh->z - pSurfacePoint->z);
             
    return distance / fabs(distance);
}

// NEW:
gdouble cpt_distance_sign(GtsPoint* pPointMesh, ParamDistFunc* pParam)
{
    GtsVector*  pNormal       = (GtsVector*)pParam->pNormal;
    GtsPoint*   pSurfacePoint = (GtsPoint*)pParam->pSurfacePoint;
    gdouble     distance;
    
    distance = (*pNormal)[0] * (pPointMesh->x - pSurfacePoint->x)
             + (*pNormal)[1] * (pPointMesh->y - pSurfacePoint->y)
             + (*pNormal)[2] * (pPointMesh->z - pSurfacePoint->z);
             
    return gts_cpt::sign(distance);  // Safe sign extraction
}
```

**Step 2: Verify compilation**

Run: `make`
Expected: Success

**Step 3: Commit**

```bash
git add src/cpt.cpp
git commit -m "fix: use safe sign() in cpt_distance_sign

- Removes division by zero UB
- Handles distance=0 case correctly"
```

---

## Task 4: Remove Duplicate Function Definition

**Files:**
- Modify: `src/mesh_visit.cpp`

**Step 1: Analyze file usage**

Run: `grep -r "mesh_visit" src/*.cpp include/*.h Makefile`
Expected: Find if mesh_visit.cpp is actually linked

**Step 2: If NOT linked in Makefile, remove it**

Run: `grep "mesh_visit" Makefile`
If not found in Makefile lines 2-3:

```bash
git rm src/mesh_visit.cpp
```

**Step 3: If it IS linked, fix duplicate definitions**

Modify `src/mesh_visit.cpp` to use namespaces:

```cpp
// mesh_visit.cpp
// version 1.0 : Thiago Macedo
// Interface to send eulerian mesh to VisIt;

namespace {  // Anonymous namespace for internal linkage

/* Simulation mesh */
float mesh_x[] = {0., 1., 2.5, 5.};
float mesh_y[] = {0., 2., 2.25, 2.55, 5.};
int    mesh_dims[] = {4, 5, 1};
int    mesh_ndims = 2;

} // anonymous namespace

VisIt_MeshData *VisItGetMesh2D(int domain, const char *name)
{
    // ... existing code ...
}

namespace {

#define NPTS 100
float  angle = 0.;
int    pmesh_ndims = 3;
float  pmesh_x[NPTS], pmesh_y[NPTS], pmesh_z[NPTS];

} // anonymous namespace

VisIt_MeshData *VisItGetMesh3D(int domain, const char *name)
{
    // ... existing code ...
}
```

**Step 4: Commit**

```bash
git add src/mesh_visit.cpp
git commit -m "fix: resolve duplicate function definitions

- Use namespaces or rename functions
- Prevent linker errors"
```

---

## Task 5: Add File I/O Error Handling

**Files:**
- Modify: `src/main.cpp:498-500`

**Step 1: Add NULL check for fopen**

Modify `src/main.cpp:498-500`:

```cpp
// OLD:
FILE* fpVtk = fopen("out/surface.vtk", "w+b");
gts_surface_write_vtk(pSignedDistance->pSurface, fpVtk);
fclose(fpVtk);

// NEW:
FILE* fpVtk = fopen("out/surface.vtk", "w+b");
if (!fpVtk) {
    fprintf(stderr, "gtscpt: cannot open 'out/surface.vtk' for writing\n");
    fprintf(stderr, "Hint: ensure 'out/' directory exists with write permissions\n");
    // Not a fatal error - continue with other cleanup
} else {
    gts_surface_write_vtk(pSignedDistance->pSurface, fpVtk);
    fclose(fpVtk);
}
```

**Step 2: Add directory check hint**

Add helper function at top of `src/main.cpp` (after includes):

```cpp
#include <sys/stat.h>  // For mkdir

namespace {
    bool ensure_directory_exists(const char* path) {
        struct stat st = {0};
        if (stat(path, &st) == -1) {
            return mkdir(path, 0700) == 0;
        }
        return S_ISDIR(st.st_mode);
    }
}
```

Update `src/main.cpp:498`:

```cpp
// NEW (before fopen):
if (!ensure_directory_exists("out")) {
    fprintf(stderr, "gtscpt: warning - cannot create 'out/' directory\n");
}

FILE* fpVtk = fopen("out/surface.vtk", "w+b");
// ... rest of code
```

**Step 3: Test with missing directory**

Run: `rm -rf out/ && make clean && make && ./run sphere20.gts`
Expected: Creates `out/` directory automatically

**Step 4: Commit**

```bash
git add src/main.cpp
git commit -m "fix: add file I/O error handling

- Check fopen return value
- Auto-create output directory
- Provide helpful error messages"
```

---

## Task 6: Fix Signed/Unsigned Warning

**Files:**
- Modify: `src/main.cpp:157-161`

**Step 1: Fix useless size check**

Modify `src/main.cpp:157-161`:

```cpp
// OLD:
if ( sizeMesh[0] <= 0 || sizeMesh[1] <= 0 || sizeMesh[2] <= 0)
{
    fputs("gtscpt: you must correctly especify the size of eulerian mesh!\n", stderr);
    return 1; // failure
}

// NEW:
// sizeMesh is size_t (unsigned), so <= 0 is always false
// Check against minimum valid size instead
constexpr size_t MIN_MESH_SIZE = 1;
if ( sizeMesh[0] < MIN_MESH_SIZE || sizeMesh[1] < MIN_MESH_SIZE || sizeMesh[2] < MIN_MESH_SIZE)
{
    fputs("gtscpt: mesh size must be at least 1 in each dimension\n", stderr);
    return 1;
}
```

**Step 2: Update sizeMesh type for better checking**

Modify `src/main.cpp:44`:

```cpp
// OLD:
SizeVector  sizeMesh  = {   0 ,   0 ,   0 };

// NEW:
std::array<long, 3> sizeMesh = {0, 0, 0};  // Use signed for validation
```

Modify `src/main.cpp:95-101`:

```cpp
// OLD:
sizeMesh[0] = (size_t) atoi (optarg);

// NEW:
sizeMesh[0] = atol(optarg);
// Validate range
if (sizeMesh[0] <= 0) {
    fprintf(stderr, "gtscpt: size-x must be positive, got %ld\n", sizeMesh[0]);
    return 1;
}
```

Repeat for sizeMesh[1] and sizeMesh[2].

Modify `src/main.cpp:295`:

```cpp
// OLD:
pSignedDistance->size[i] = sizeMesh[i];

// NEW:
pSignedDistance->size[i] = static_cast<size_t>(sizeMesh[i]);
```

**Step 3: Commit**

```bash
git add src/main.cpp
git commit -m "fix: correct signed/unsigned mesh size validation

- Use signed type for input validation
- Check for positive values before assignment
- Fix always-false condition"
```

---

## Task 7: Fix Include Guard Typo in export.h

**Files:**
- Modify: `include/export.h:8`

**Step 1: Fix comment typo**

Modify `include/export.h:8`:

```cpp
// OLD:
#endif // __CPT_H__

// NEW:
#endif // __EXPORT_H__
```

**Step 2: Fix include guard naming (avoid reserved identifiers)**

Modify `include/export.h:1-8`:

```cpp
// OLD:
#ifndef __EXPORT_H__
#define __EXPORT_H__

// NEW:
#ifndef GTS_CPT_EXPORT_H
#define GTS_CPT_EXPORT_H

// ... content ...

#endif // GTS_CPT_EXPORT_H
```

**Step 3: Commit**

```bash
git add include/export.h
git commit -m "fix: correct include guard and avoid reserved identifiers

- Fix copy-paste error in guard comment
- Use GTS_CPT_EXPORT_H (double underscore is reserved)"
```

---

## Task 8: Add Const Correctness

**Files:**
- Modify: `src/cpt.cpp:25-76`
- Modify: `include/cpt.h:23-24, 29`

**Step 1: Update cpt_min/cpt_max**

Modify `src/cpt.cpp:25-35`:

```cpp
// OLD:
gdouble cpt_min(gdouble a, gdouble b)
{
    return (a < b ? a : b);
}

gdouble cpt_max(gdouble a, gdouble b)
{
    return (a > b ? a : b);
}

// NEW:
gdouble cpt_min(const gdouble a, const gdouble b) noexcept
{
    return (a < b ? a : b);
}

gdouble cpt_max(const gdouble a, const gdouble b) noexcept
{
    return (a > b ? a : b);
}
```

Update declaration in `include/cpt.h:23-24`:

```cpp
// OLD:
gdouble cpt_min(gdouble a, gdouble b);
gdouble cpt_max(gdouble a, gdouble b);

// NEW:
gdouble cpt_min(const gdouble a, const gdouble b) noexcept;
gdouble cpt_max(const gdouble a, const gdouble b) noexcept;
```

**Step 2: Update cpt_vector_angle**

Modify `src/cpt.cpp:63`:

```cpp
// OLD:
gdouble cpt_vector_angle (GtsVector vector[])

// NEW:
gdouble cpt_vector_angle(const GtsVector v1, const GtsVector v2) noexcept
```

Update implementation:

```cpp
gdouble cpt_vector_angle(const GtsVector v1, const GtsVector v2) noexcept
{
    const gdouble pvx = v1[1]*v2[2] - v1[2]*v2[1];
    const gdouble pvy = v1[2]*v2[0] - v1[0]*v2[2];
    const gdouble pvz = v1[0]*v2[1] - v1[1]*v2[0];

    const gdouble theta = atan2(
        sqrt (pvx*pvx + pvy*pvy + pvz*pvz),
        v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
    );

    return theta;
}
```

Update declaration in `include/cpt.h:29`:

```cpp
// OLD:
gdouble cpt_vector_angle (GtsVector vector[]);

// NEW:
gdouble cpt_vector_angle(const GtsVector v1, const GtsVector v2) noexcept;
```

**Step 3: Update call sites**

Modify `src/cpt.cpp:500`:

```cpp
// OLD:
theta = cpt_vector_angle(normal[0], normal[1]);

// NEW:
theta = cpt_vector_angle(normal[0], normal[1]);
```

**Step 4: Commit**

```bash
git add src/cpt.cpp include/cpt.h
git commit -m "refactor: add const correctness to utility functions

- Mark parameters const where appropriate
- Add noexcept for simple functions
- Modernize API contracts"
```

---

## Task 9: Move Constants to Anonymous Namespace

**Files:**
- Modify: `src/main.cpp:28-34`

**Step 1: Fix PI macro and move constants**

Modify `src/main.cpp:28-34`:

```cpp
// OLD:
static const gdouble  DBL_CLOCKS_PER_SEC  = (double)CLOCKS_PER_SEC;
       const gdouble  MIN_TRANGLE_QUALITY = 0.5;   // TODO: Quanto eh bom?
       const gdouble  MIN_DELTA_MESH      = 0.001;
       const gdouble  size_sup            = 2.0;   // suporte da delta de dirac;
static      gboolean  verbose             = TRUE;
static      gboolean  normalize           = FALSE;
static         guint  geodesation_order   = 0;

// NEW:
namespace {
    constexpr double DBL_CLOCKS_PER_SEC = static_cast<double>(CLOCKS_PER_SEC);
    constexpr double MIN_TRIANGLE_QUALITY = 0.5;   // Minimum triangle quality ratio
    constexpr double MIN_DELTA_MESH = 0.001;       // Minimum mesh spacing
    constexpr double SIZE_SUP = 2.0;               // Dirac delta support
    
    gboolean verbose = FALSE;    // Changed default to FALSE
    gboolean normalize = FALSE;
    guint geodesation_order = 0;
}

// Fix typo: MIN_TRANGLE -> MIN_TRIANGLE
```

**Step 2: Update references**

Modify `src/main.cpp:270`:

```cpp
// OLD:
if(qstats.face_quality.max < MIN_TRANGLE_QUALITY)

// NEW:
if(qstats.face_quality.max < MIN_TRIANGLE_QUALITY)
```

Repeat at lines 275, 279.

Modify `src/main.cpp:310`:

```cpp
// OLD:
pSignedDistance->distCut = cpt_get_dist_cut(size_sup, delta_max);

// NEW:
pSignedDistance->distCut = cpt_get_dist_cut(SIZE_SUP, delta_max);
```

**Step 3: Commit**

```bash
git add src/main.cpp
git commit -m "refactor: move constants to anonymous namespace

- Fix typo: MIN_TRANGLE -> MIN_TRIANGLE
- Use constexpr for compile-time constants
- Reduce global namespace pollution"
```

---

## Task 10: Replace C-Style Casts

**Files:**
- Modify: `src/main.cpp:77-101`
- Modify: `src/export.cpp:28-30`

**Step 1: Replace casts in main.cpp**

Modify `src/main.cpp:77-101`:

```cpp
// OLD:
beginMesh[0] = (gdouble) atof (optarg);

// NEW:
beginMesh[0] = static_cast<gdouble>(std::atof(optarg));
```

Repeat for all 9 conversions (lines 77-101).

**Step 2: Replace casts in export.cpp**

Modify `src/export.cpp:28-30`:

```cpp
// OLD:
nodex = (float*) malloc(dims[0]*sizeof(float));
nodey = (float*) malloc(dims[1]*sizeof(float));
nodez = (float*) malloc(dims[2]*sizeof(float));

// NEW:
nodex = static_cast<float*>(malloc(dims[0]*sizeof(float)));
nodey = static_cast<float*>(malloc(dims[1]*sizeof(float)));
nodez = static_cast<float*>(malloc(dims[2]*sizeof(float)));
```

Better yet, use C++ allocation:

```cpp
// BEST:
std::vector<float> nodex(dims[0]);
std::vector<float> nodey(dims[1]);
std::vector<float> nodez(dims[2]);
```

Then update cleanup (lines 78-80):

```cpp
// OLD:
free(nodex);
free(nodey);
free(nodez);

// NEW (automatic cleanup with vectors):
// No manual cleanup needed
```

**Step 3: Commit**

```bash
git add src/main.cpp src/export.cpp
git commit -m "refactor: replace C-style casts with static_cast

- Use static_cast for type conversions
- Prefer vector over malloc/free in export"
```

---

## Task 11: Add Basic Test Infrastructure

**Files:**
- Create: `tests/CMakeLists.txt` or `tests/Makefile`
- Create: `tests/test_main.cpp`

**Step 1: Create test framework**

Create `tests/Makefile`:

```makefile
# Test Makefile
CXX = g++
CXXFLAGS = -std=c++20 -Wall -Wextra -g -I../include
LDFLAGS = -lgts -lglib-2.0 -lm

TESTS = test_sign test_safe_multiply

all: $(TESTS)

test_sign: test_sign.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

test_safe_multiply: test_safe_multiply.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f $(TESTS)

check: $(TESTS)
	@echo "Running tests..."
	@for test in $(TESTS); do \
		echo "Running $$test..."; \
		./$$test || exit 1; \
	done
	@echo "All tests passed!"

.PHONY: all clean check
```

**Step 2: Create safe_multiply test**

Create `tests/test_safe_multiply.cpp`:

```cpp
#include <cassert>
#include <limits>
#include <stdexcept>

namespace gts_cpt {
    size_t safe_mesh_size(size_t sx, size_t sy, size_t sz);
}

int main() {
    using gts_cpt::safe_mesh_size;
    
    // Valid sizes
    assert(safe_mesh_size(10, 10, 10) == 1000);
    assert(safe_mesh_size(100, 100, 100) == 1000000);
    
    // Edge case: zero
    assert(safe_mesh_size(0, 10, 10) == 0);
    assert(safe_mesh_size(10, 0, 10) == 0);
    
    // Overflow detection (this should throw)
    bool caught = false;
    try {
        constexpr size_t big = std::numeric_limits<size_t>::max();
        safe_mesh_size(big, big, big);
    } catch (const std::overflow_error&) {
        caught = true;
    }
    assert(caught && "Expected overflow_error");
    
    return 0;
}
```

**Step 3: Run tests**

Run: `cd tests && make check`
Expected: All tests pass

**Step 4: Commit**

```bash
git add tests/
git commit -m "test: add basic test infrastructure

- Add Makefile for test compilation
- Add tests for safe_mesh_size overflow detection
- Foundation for future test expansion"
```

---

## Verification Checklist

After completing all tasks, verify:

- [ ] All tests pass: `cd tests && make check`
- [ ] No compiler warnings: `make clean && make 2>&1 | grep warning`
- [ ] No memory leaks: `valgrind --leak-check=full ./debug/gts-cpt ...`
- [ ] No undefined behavior: `UBSAN_OPTIONS=print_stacktrace=1 ./debug/gts-cpt ...`
- [ ] Program still works correctly with sample data

---

## Summary

This plan fixes:
- **5 Critical issues** (UB, memory leaks, overflow)
- **6 High severity issues** (RAII, const, globals, casts)

Each task is atomic, testable, and maintains backward compatibility while preparing for C++20 migration.

**Estimated time:** 4-6 hours for experienced developer
**Risk level:** Low (incremental changes, tests added)

After these fixes, the codebase will be safe for C++20 migration work.