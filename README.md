# GTS-CPT: Signed Distance Function using Closest Point Transform

An implementation of a signed distance function based on Mauch's Closest Point Transform (CPT) employing the GNU Triangulated Surface Library.

## Overview

The signed distance function between an arbitrary point in 3D space and a given closed surface returns the minimum distance from that point to the collection of triangles representing the surface. By convention, the sign is positive if the point is outside and negative if the point is inside the region determined by the surface.

In Computational Fluid Dynamics, the signed distance function is useful to:
- Locate the separation interface between two fluids (zero level-set surface)
- Determine instantaneous mass density and viscosity for each Eulerian grid point

## Requirements

- **C++17** compatible compiler (GCC 8+, Clang 7+)
- **GTS Library** - GNU Triangulated Surface Library
- **GLib 2.0** - Required by GTS
- **Silo** - For mesh output (optional, for VisIt visualization)

### Installing Dependencies (Arch Linux)

```bash
sudo pacman -S gts silo
```

### Installing Dependencies (Ubuntu/Debian)

```bash
sudo apt install libgts-dev silo-bin
```

## Building

```bash
make clean && make
```

The binary will be created at `debug/gts-cpt`.

### Build Options

The Makefile uses these default flags:
- `-std=c++17` - C++17 standard
- `-Wall -pedantic` - Enable warnings
- `-O3 -funroll-loops -ftree-vectorize` - Optimizations

## Usage

### Basic Syntax

```bash
./debug/gts-cpt [OPTIONS] < input_file.gts
```

### Required Options

All mesh boundary and size parameters must be specified:

| Option | Description |
|--------|-------------|
| `--begin-x VALUE` | Minimum X coordinate of Eulerian mesh |
| `--begin-y VALUE` | Minimum Y coordinate of Eulerian mesh |
| `--begin-z VALUE` | Minimum Z coordinate of Eulerian mesh |
| `--end-x VALUE` | Maximum X coordinate of Eulerian mesh |
| `--end-y VALUE` | Maximum Y coordinate of Eulerian mesh |
| `--end-z VALUE` | Maximum Z coordinate of Eulerian mesh |
| `--size-x N` | Number of cells along X axis |
| `--size-y N` | Number of cells along Y axis |
| `--size-z N` | Number of cells along Z axis |

### Optional Flags

| Flag | Description |
|------|-------------|
| `--normalize` | Fit surface in a unit cube centered at origin |
| `--verbose` | Print statistics and progress |
| `--help` | Display help message |

### Example: Unit Sphere

```bash
# Using a sphere with geodesic order 4
./debug/gts-cpt --verbose \
    --begin-x -1.5 --begin-y -1.5 --begin-z -1.5 \
    --end-x 1.5 --end-y 1.5 --end-z 1.5 \
    --size-x 100 --size-y 100 --size-z 100 \
    < sphere20.gts
```

### Example: Custom Mesh

```bash
# Process a GTS file with custom extents
./debug/gts-cpt --verbose --normalize \
    --begin-x -2.0 --begin-y -2.0 --begin-z -2.0 \
    --end-x 2.0 --end-y 2.0 --end-z 2.0 \
    --size-x 50 --size-y 50 --size-z 50 \
    < my_mesh.gts
```

### Using the Run Script

```bash
# The 'run' script simplifies execution
./run sphere20.gts

# Or with STL files (auto-converted)
./run my_model.stl
```

## Output

The program creates output in the `out/` directory:

| File | Format | Description |
|------|--------|-------------|
| `eulerianmesh.silo` | Silo | Distance function for VisIt visualization |
| `surface.vtk` | VTK | Surface mesh for visualization |

### Visualizing Results

Using VisIt:
```bash
visit -o out/eulerianmesh.silo
```

Using ParaView:
```bash
paraview out/surface.vtk
```

## Testing

The test suite covers critical algorithms and safety checks:

```bash
cd tests && make check
```

### Test Coverage

| Test File | Coverage |
|-----------|----------|
| `test_sign.cpp` | Sign function (zero, bounds, subnormals) |
| `test_safe_mesh_size.cpp` | Overflow detection, zero handling |
| `test_cpt_math.cpp` | min/max, vector angle, distance thresholds |

## Project Structure

```
gts-cpt/
├── include/
│   ├── gtscpt.h      # Main header, SignedDistance class
│   ├── cpt.h         # CPT function declarations
│   ├── export.h      # Export function declarations
│   ├── gtstools.h    # GTS utility declarations
│   ├── array3d.h     # Legacy array template (unused)
│   └── gts-cpp/      # C++ RAII wrappers
│       ├── gts_wrapper.hpp   # Template deleters for GTS objects
│       ├── surface.hpp       # gts::Surface RAII wrapper
│       ├── vertex.hpp        # gts::Vertex, gts::Edge wrappers
│       └── callbacks.hpp     # Callback utilities
├── src/
│   ├── main.cpp      # Entry point, CLI parsing
│   ├── cpt.cpp       # CPT algorithm implementation
│   ├── export.cpp    # Silo/VTK export functions
│   ├── gtstools.cpp  # GTS utility functions
│   └── gts-cpp/      # C++ wrapper implementations
│       └── surface.cpp
├── tests/
│   ├── Makefile
│   ├── test_sign.cpp
│   ├── test_safe_mesh_size.cpp
│   └── test_cpt_math.cpp
├── docs/plans/
│   └── *.md          # Design and implementation plans
├── Makefile
├── run               # Execution helper script
└── README.md
```

## Technical Notes

### Recent Improvements

This codebase has been modernized from C++98 to C++17 with the following fixes:

1. **Memory Safety**
   - RAII wrapper for `SignedDistance` class
   - `gts::Surface` wrapper with `std::unique_ptr` ownership
   - Template `GtsDeleter<T>` for automatic GTS object cleanup
   - `std::unique_ptr` and `std::vector` for automatic cleanup
   - Overflow-safe mesh size allocation

2. **Undefined Behavior Fixes**
   - Safe `sign()` function using `std::copysign` (no division by zero)
   - Division by zero protection in `cpt_point_angle()` with epsilon check
   - File I/O error handling with NULL checks
   - Proper exception handling

3. **Code Quality**
   - Const correctness for utility functions
   - `noexcept` for simple functions
   - Proper include guards
   - Removed dead code and declarations
   - Clean compile with `-Werror`

### Supported Input Formats

- **GTS** - Native GTS triangulated surface format
- **STL** - STL format (converted via `stl2gts` in run script)

### Performance

For surfaces where the number of Lagrangian points (triangle vertices) is less than N*log(N) of Eulerian grid points, the CPT algorithm provides:
- Linear complexity in geometric elements
- O(1) proportion constant per element

## Sample Data

Download GTS sample files from:
- http://gts.sourceforge.net/samples.html

Note: The program requires **closed manifold** surfaces. Non-manifold surfaces may produce incorrect results.

## References

If you use this implementation, please cite:

1. **Ceniceros and Roma 2005**: "A study of the long-time dynamics of a three-dimensional drop rising through a matching-density medium"
   - Journal of Computational Physics, Vol. 205, Issue 2, May 2005, Pages 391-400
   - DOI: [10.1016/j.jcp.2004.11.013](https://doi.org/10.1016/j.jcp.2004.11.013)

2. **Ceniceros et al. 2009**: "A fast, robust, and accurate adaptive method for the simulation of surface tension flows"
   - Commun. Comput. Phys., Vol. 8, 2010, Pages 51-94
   - DOI: [10.4208/cicp.050509.141009a](https://doi.org/10.4208/cicp.050509.141009a)

3. **Mauch's Algorithm**: Closest Point Transform
   - http://www.cacr.caltech.edu/~sean/projects/cpt/html3/index.html

## Authors

- **Thiago P. Macedo** - Original implementation and maintenance
- **Prof. Alexandre M. Roma** - Supervision
- Instituto de Matemática e Estatística, Universidade de São Paulo (IME-USP)

## License

See [LICENSE](LICENSE) file for details.

## Bug Reports

Please report issues to:
- Thiago P. Macedo (tmacedo@usp.br)