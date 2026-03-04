#pragma once

#include "gts_wrapper.hpp"
#include <gts.h>
#include <memory>
#include <stdexcept>
#include <string>

namespace gts {

/**
 * @class Surface
 * @brief RAII wrapper for GtsSurface
 * 
 * Provides automatic resource management for GTS surface objects.
 * The surface is automatically destroyed when the wrapper goes out of scope.
 */
class Surface {
public:
    /**
     * @brief Factory function - creates an empty surface
     * @return unique_ptr to new Surface
     * @throws std::runtime_error if surface creation fails
     */
    static std::unique_ptr<Surface> create();
    
    /**
     * @brief Factory function - reads surface from file
     * @param filename Path to GTS file
     * @return unique_ptr to new Surface
     * @throws std::runtime_error if file cannot be opened or parsed
     */
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
    
    // Operator for implicit conversion
    operator GtsSurface*() noexcept { return surface_.get(); }
    operator const GtsSurface*() const noexcept { return surface_.get(); }
    
    // Type-safe operations
    void add_triangle(GtsVertex* v1, GtsVertex* v2, GtsVertex* v3);
    void read(const char* filename);
    void read(const std::string& filename);
    
    // Surface properties
    bool is_closed() const;
    gdouble volume() const;
    gdouble area() const;
    
    // Check if surface is valid
    bool is_valid() const noexcept { return surface_ != nullptr; }
    
private:
    Surface() = default;
    gts_unique_ptr<GtsSurface> surface_;
};

/**
 * @brief RAII helper for reading surface from file
 * @param filename Path to GTS file
 * @return unique_ptr to Surface
 * @throws std::runtime_error on failure
 */
std::unique_ptr<Surface> read_surface(const char* filename);

} // namespace gts
