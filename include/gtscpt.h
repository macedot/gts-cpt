#ifndef __GTSCPT_H__
#define __GTSCPT_H__

// Standard C++ library
#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <cmath>
#include <vector>
#include <memory>
#include <limits>
#include <stdexcept>
// Blitz related
//#include <blitz/numinquire.h>
//#include <blitz/array.h>

#ifdef __cplusplus
extern "C" {
#endif

// Standard C library
#include <stdlib.h>
#include <locale.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <memory.h>

// GTS related
#include "gts.h"
#include "gtstools.h"

#ifdef __cplusplus
}
#endif

// Safe sign function - replaces dangerous SIGN macro
namespace gts_cpt {
    [[nodiscard]] inline constexpr double sign(double x) noexcept {
        // Use copysign: returns value with magnitude of 1.0 and sign of x
        // For x = 0.0, copysign(1.0, +0.0) = 1.0
        // For x = -0.0, copysign(1.0, -0.0) = -1.0
        return std::copysign(1.0, x);
    }
    
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

// Compatibility macro for gradual migration
#define SIGN(x) gts_cpt::sign(x)


using SizeVector = size_t[3];

// RAII class for signed distance mesh data
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
        for (int i = 0; i < 3; i++) {
            coordMin[i] = other.coordMin[i];
            coordMax[i] = other.coordMax[i];
            delta[i] = other.delta[i];
            size[i] = other.size[i];
        }
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
            for (int i = 0; i < 3; i++) {
                coordMin[i] = other.coordMin[i];
                coordMax[i] = other.coordMax[i];
                delta[i] = other.delta[i];
                size[i] = other.size[i];
            }
        }
        return *this;
    }
};


#endif // __GTSCPT_H__
