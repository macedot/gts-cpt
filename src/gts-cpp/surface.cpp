#include "gts-cpp/surface.hpp"
#include <cstdio>
#include <stdexcept>

namespace gts {

std::unique_ptr<Surface> Surface::create() {
    auto surf = std::unique_ptr<Surface>(new Surface());
    surf->surface_ = adopt(
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
    // Create edges and face
    // Note: GTS takes ownership of edges and face when added to surface.
    // We use RAII wrappers to ensure cleanup if add_face fails.
    auto e1 = adopt(GTS_EDGE(gts_edge_new(gts_edge_class(), v1, v2)));
    auto e2 = adopt(GTS_EDGE(gts_edge_new(gts_edge_class(), v2, v3)));
    auto e3 = adopt(GTS_EDGE(gts_edge_new(gts_edge_class(), v3, v1)));
    auto face = adopt(GTS_EDGE(gts_face_new(gts_face_class(), e1.get(), e2.get(), e3.get())));
    
    // Add to surface - GTS takes ownership, release our pointers
    gts_surface_add_face(surface_.get(), GTS_FACE(face.get()));
    e1.release();
    e2.release();
    e3.release();
    face.release();
}

void Surface::read(const char* filename) {
    FILE* fp = std::fopen(filename, "r");
    if (!fp) {
        throw std::runtime_error(std::string("Failed to open file: ") + filename);
    }
    
    GtsFile* gf = gts_file_new(fp);
    if (gts_surface_read(surface_.get(), gf) != 0) {
        gts_file_destroy(gf);   // Destroy GtsFile first
        std::fclose(fp);        // Then close the file
        throw std::runtime_error(std::string("Failed to read surface from file: ") + filename);
    }
    
    gts_file_destroy(gf);
    std::fclose(fp);
}

void Surface::read(const std::string& filename) {
    read(filename.c_str());
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

} // namespace gts
