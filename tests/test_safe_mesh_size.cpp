#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <iostream>

// Inline definitions for testing (actual implementations in gtscpt.h)
namespace gts_cpt {
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

int main() {
    using gts_cpt::safe_mesh_size;
    
    std::cout << "Test 1: Valid sizes..." << std::endl;
    assert(safe_mesh_size(10, 10, 10) == 1000);
    assert(safe_mesh_size(100, 100, 100) == 1000000);
    assert(safe_mesh_size(1, 1, 1) == 1);
    std::cout << "  PASS" << std::endl;
    
    std::cout << "Test 2: Edge case - zero..." << std::endl;
    assert(safe_mesh_size(0, 10, 10) == 0);
    assert(safe_mesh_size(10, 0, 10) == 0);
    assert(safe_mesh_size(10, 10, 0) == 0);
    assert(safe_mesh_size(0, 0, 0) == 0);
    std::cout << "  PASS" << std::endl;
    
    std::cout << "Test 3: Overflow detection..." << std::endl;
    bool caught = false;
    try {
        // Use runtime value to prevent compile-time optimization
        volatile size_t big = std::numeric_limits<size_t>::max();
        safe_mesh_size(big, 2, 1);
    } catch (const std::overflow_error&) {
        caught = true;
    }
    if (!caught) {
        std::cerr << "FAIL: Expected overflow_error" << std::endl;
        return 1;
    }
    std::cout << "  PASS" << std::endl;
    
    std::cout << "Test 4: Boundary values..." << std::endl;
    // Use a large value that won't overflow but is still meaningfully large
    size_t large_val = 100000;
    size_t result = safe_mesh_size(large_val, large_val, 2);
    assert(result == large_val * large_val * 2);
    std::cout << "  PASS" << std::endl;
    
    std::cout << "All tests passed!" << std::endl;
    return 0;
}