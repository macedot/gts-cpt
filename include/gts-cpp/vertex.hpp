#pragma once

#include "gts_wrapper.hpp"
#include <gts.h>
#include <memory>

namespace gts {

/**
 * @class Vertex
 * @brief RAII wrapper for GtsVertex
 * 
 * Provides automatic resource management for GTS vertex objects.
 */
class Vertex {
public:
    /**
     * @brief Create vertex at coordinates
     * @param x X coordinate
     * @param y Y coordinate
     * @param z Z coordinate
     * @return unique_ptr to new Vertex
     */
    static std::unique_ptr<Vertex> create(gdouble x, gdouble y, gdouble z);
    
    // Non-copyable
    Vertex(const Vertex&) = delete;
    Vertex& operator=(const Vertex&) = delete;
    
    // Movable
    Vertex(Vertex&& other) noexcept;
    Vertex& operator=(Vertex&& other) noexcept;
    
    // Destructor
    ~Vertex();
    
    // Access
    GtsVertex* get() noexcept { return vertex_.get(); }
    const GtsVertex* get() const noexcept { return vertex_.get(); }
    
    // Operator for implicit conversion
    operator GtsVertex*() noexcept { return vertex_.get(); }
    operator const GtsVertex*() const noexcept { return vertex_.get(); }
    
    // Coordinates
    gdouble x() const noexcept;
    gdouble y() const noexcept;
    gdouble z() const noexcept;
    
    // Check if valid
    bool is_valid() const noexcept { return vertex_ != nullptr; }
    
private:
    Vertex() = default;
    gts_unique_ptr<GtsVertex> vertex_;
};

/**
 * @class Edge
 * @brief RAII wrapper for GtsEdge
 * 
 * Provides automatic resource management for GTS edge objects.
 */
class Edge {
public:
    /**
     * @brief Create edge between two vertices
     * @param v1 First vertex
     * @param v2 Second vertex
     * @return unique_ptr to new Edge
     */
    static std::unique_ptr<Edge> create(GtsVertex* v1, GtsVertex* v2);
    
    // Non-copyable
    Edge(const Edge&) = delete;
    Edge& operator=(const Edge&) = delete;
    
    // Movable
    Edge(Edge&& other) noexcept;
    Edge& operator=(Edge&& other) noexcept;
    
    // Destructor
    ~Edge();
    
    // Access
    GtsEdge* get() noexcept { return edge_.get(); }
    const GtsEdge* get() const noexcept { return edge_.get(); }
    
    // Operator for implicit conversion
    operator GtsEdge*() noexcept { return edge_.get(); }
    operator const GtsEdge*() const noexcept { return edge_.get(); }
    
    // Check if valid
    bool is_valid() const noexcept { return edge_ != nullptr; }
    
private:
    Edge() = default;
    gts_unique_ptr<GtsEdge> edge_;
};

} // namespace gts
