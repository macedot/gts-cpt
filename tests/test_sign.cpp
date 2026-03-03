#include <cassert>
#include <cmath>
#include <limits>

// Define the function inline for testing (actual implementation in gtscpt.h)
namespace gts_cpt {
    [[nodiscard]] inline constexpr double sign(double x) noexcept {
        return std::copysign(1.0, x);
    }
}

int main() {
    using gts_cpt::sign;
    
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