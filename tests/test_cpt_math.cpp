#include <cassert>
#include <cmath>
#include <limits>

// Standalone implementations for testing (no GTS dependency)
namespace {

double cpt_min(const double a, const double b) noexcept {
    return (a < b ? a : b);
}

double cpt_max(const double a, const double b) noexcept {
    return (a > b ? a : b);
}

double cpt_vector_angle(const double v1[3], const double v2[3]) noexcept {
    const double pvx = v1[1]*v2[2] - v1[2]*v2[1];
    const double pvy = v1[2]*v2[0] - v1[0]*v2[2];
    const double pvz = v1[0]*v2[1] - v1[1]*v2[0];

    const double theta = atan2(
        sqrt(pvx*pvx + pvy*pvy + pvz*pvz),
        v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
    );

    return theta;
}

double cpt_get_dist_cut(double size_sup, double delta_max) noexcept {
    constexpr double SQRT_3 = 1.7320508075688772;
    return SQRT_3 * 0.5 * size_sup * delta_max;
}

double cpt_get_dist_max(double dist_cut, double dist_extra) noexcept {
    constexpr double c = 1.0;
    return dist_cut + c * dist_extra;
}

}

int main() {
    // ========================================
    // Test cpt_min
    // ========================================
    assert(cpt_min(1.0, 2.0) == 1.0);
    assert(cpt_min(2.0, 1.0) == 1.0);
    assert(cpt_min(-1.0, 1.0) == -1.0);
    assert(cpt_min(0.0, 0.0) == 0.0);
    
    // NaN handling (min should return the non-NaN or NaN behavior)
    double nan = std::numeric_limits<double>::quiet_NaN();
    assert(std::isnan(cpt_min(nan, 1.0)) || cpt_min(nan, 1.0) == 1.0);
    
    // ========================================
    // Test cpt_max
    // ========================================
    assert(cpt_max(1.0, 2.0) == 2.0);
    assert(cpt_max(2.0, 1.0) == 2.0);
    assert(cpt_max(-1.0, 1.0) == 1.0);
    assert(cpt_max(0.0, 0.0) == 0.0);
    
    // ========================================
    // Test cpt_vector_angle
    // ========================================
    // Parallel vectors (same direction) -> 0 radians
    double v1[3] = {1.0, 0.0, 0.0};
    double v2[3] = {1.0, 0.0, 0.0};
    double angle = cpt_vector_angle(v1, v2);
    assert(std::abs(angle - 0.0) < 1e-10);
    
    // Perpendicular vectors -> pi/2 radians
    double v3[3] = {0.0, 1.0, 0.0};
    angle = cpt_vector_angle(v1, v3);
    assert(std::abs(angle - M_PI/2) < 1e-10);
    
    // Opposite vectors -> pi radians
    double v4[3] = {-1.0, 0.0, 0.0};
    angle = cpt_vector_angle(v1, v4);
    assert(std::abs(angle - M_PI) < 1e-10);
    
    // Arbitrary vectors
    double v5[3] = {1.0, 1.0, 0.0};  // 45 degrees in XY plane
    angle = cpt_vector_angle(v1, v5);
    assert(std::abs(angle - M_PI/4) < 1e-10);
    
    // Normalized vectors
    double v6[3] = {0.5, 0.5, 0.7071067811865476};  // normalized (1,1,1)/sqrt(3)
    double v7[3] = {0.5, 0.5, 0.7071067811865476};
    angle = cpt_vector_angle(v6, v7);
    assert(std::abs(angle - 0.0) < 1e-10);
    
    // ========================================
    // Test cpt_get_dist_cut
    // ========================================
    // Formula: sqrt(3) * 0.5 * size_sup * delta_max
    double dist_cut = cpt_get_dist_cut(2.0, 1.0);
    double expected = 1.7320508075688772 * 0.5 * 2.0 * 1.0;
    assert(std::abs(dist_cut - expected) < 1e-10);
    
    dist_cut = cpt_get_dist_cut(4.0, 0.5);
    expected = 1.7320508075688772 * 0.5 * 4.0 * 0.5;
    assert(std::abs(dist_cut - expected) < 1e-10);
    
    // Zero inputs
    assert(cpt_get_dist_cut(0.0, 1.0) == 0.0);
    assert(cpt_get_dist_cut(1.0, 0.0) == 0.0);
    
    // ========================================
    // Test cpt_get_dist_max
    // ========================================
    // Formula: dist_cut + 1.0 * dist_extra
    double dist_max = cpt_get_dist_max(1.0, 1.0);
    assert(std::abs(dist_max - 2.0) < 1e-10);
    
    dist_max = cpt_get_dist_max(0.5, 0.25);
    assert(std::abs(dist_max - 0.75) < 1e-10);
    
    dist_max = cpt_get_dist_max(0.0, 0.0);
    assert(dist_max == 0.0);
    
    return 0;
}