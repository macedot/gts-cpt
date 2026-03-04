#pragma once

#include <gts.h>
#include <memory>
#include <utility>
#include <type_traits>

/**
 * @file gts_wrapper.hpp
 * @brief RAII wrappers for GTS C library objects
 * 
 * This header provides C++ wrappers for GTS (GNU Triangulated Surface) library
 * objects, ensuring automatic resource management via RAII pattern.
 */

namespace gts {

/**
 * @brief Custom deleter for GTS objects
 * 
 * GTS objects are reference-counted and must be destroyed using
 * gts_object_destroy(). This deleter provides automatic cleanup
 * for std::unique_ptr.
 * 
 * @tparam T GTS object type (must inherit from GtsObject)
 */
template<typename T>
struct GtsDeleter {
    /**
     * @brief Destroy a GTS object
     * @param obj Pointer to GTS object (can be nullptr)
     */
    void operator()(T* obj) const noexcept {
        if (obj != nullptr) {
            gts_object_destroy(GTS_OBJECT(obj));
        }
    }
};

/**
 * @brief Type alias for unique_ptr with GTS deleter
 * 
 * Use this for any GTS object that should be automatically
 * destroyed when the unique_ptr goes out of scope.
 * 
 * @code
 * auto surface = gts::make_unique<GtsSurface>(surface_ptr);
 * // Automatically calls gts_object_destroy when surface goes out of scope
 * @endcode
 */
template<typename T>
using gts_unique_ptr = std::unique_ptr<T, GtsDeleter<T>>;

/**
 * @brief Create a unique_ptr with GTS deleter
 * 
 * Helper function to create a gts_unique_ptr from a raw GTS pointer.
 * The unique_ptr will automatically call gts_object_destroy() on destruction.
 * 
 * @tparam T GTS object type
 * @param ptr Raw pointer to GTS object (can be nullptr)
 * @return unique_ptr with GTS deleter
 * 
 * @code
 * GtsSurface* raw = gts_surface_new(...);
 * auto surf = gts::make_unique(raw);
 * @endcode
 */
template<typename T>
gts_unique_ptr<T> make_unique(T* ptr) noexcept {
    return gts_unique_ptr<T>(ptr);
}

/**
 * @brief No-op deleter for non-owning pointers
 * 
 * Use this when you want to wrap a pointer in unique_ptr
 * but do NOT want it to be deleted.
 */
struct NoOpDeleter {
    template<typename T>
    void operator()(T*) const noexcept {}
};

/**
 * @brief Create a unique_ptr from existing object without taking ownership
 * 
 * WARNING: This is dangerous! The returned unique_ptr will NOT call
 * any destructor when it goes out of scope. This is useful only for
 * interfacing with APIs that expect unique_ptr but shouldn't own the object.
 * 
 * Prefer using raw pointers for non-owning references in modern C++.
 * 
 * @tparam T GTS object type  
 * @param ptr Raw pointer (will NOT be destroyed)
 * @return unique_ptr with no-op deleter
 */
template<typename T>
std::unique_ptr<T, NoOpDeleter> make_observer(T* ptr) noexcept {
    return std::unique_ptr<T, NoOpDeleter>(ptr);
}

/**
 * @brief Move a GTS object into unique_ptr ownership
 * 
 * Use this when you receive a GTS object from a function that
 * transfers ownership to you.
 * 
 * @tparam T GTS object type
 * @param ptr Raw pointer (ownership transferred to unique_ptr)
 * @return unique_ptr with GTS deleter
 */
template<typename T>
gts_unique_ptr<T> adopt(T* ptr) noexcept {
    return gts_unique_ptr<T>(ptr);
}

/**
 * @brief Check if a GTS object pointer is valid
 * 
 * @param ptr GTS pointer to check
 * @return true if pointer is non-null
 */
inline bool is_valid(GtsObject* ptr) noexcept {
    return ptr != nullptr;
}

} // namespace gts
