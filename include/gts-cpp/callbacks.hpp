#pragma once

/**
 * @file callbacks.hpp
 * @brief Type-safe callback wrappers for GTS C library functions
 * 
 * This header provides utilities for working with GTS C-style callbacks
 * in a type-safe manner.
 * 
 * Note: Direct use of std::function with GTS callbacks requires careful
 * lifetime management. For now, prefer the traditional GTS callback pattern:
 * 
 * @code
 * void my_callback(GtsVertex* v, gpointer data) {
 *     auto* self = static_cast<MyClass*>(data);
 *     self->process(v);
 * }
 * @endcode
 */

namespace gts {
    // Callbacks infrastructure ready for future use
    // The foreach helpers require C++20 generic lambdas or more complex 
    // workarounds to compile with strict GTS C function pointer typing.
    // For now, the traditional C-style callbacks work correctly.
} // namespace gts
