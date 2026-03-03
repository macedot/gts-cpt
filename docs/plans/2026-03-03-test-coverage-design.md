# Test Coverage Design

**Date:** 2026-03-03

## Goal

Add comprehensive tests for critical algorithms and safety-critical functions in the GTS-CPT codebase.

## Scope

Focus on **5 critical components** that prevent undefined behavior and ensure algorithmic correctness:

### Components to Test

| Component | File | Risk if Buggy | Test Priority |
|-----------|------|---------------|---------------|
| `gts_cpt::sign()` | gtscpt.h | Division by zero UB | Critical |
| `gts_cpt::safe_mesh_size()` | gtscpt.h | Overflow crash | Critical |
| `cpt_min/cpt_max` | cpt.cpp | Wrong distances | High |
| `cpt_vector_angle` | cpt.cpp | Edge angle errors | High |
| `cpt_get_dist_cut/max` | cpt.cpp | Wrong thresholds | Medium |

## Test Files

1. **test_sign.cpp** - Already exists, needs enhancement
2. **test_safe_mesh_size.cpp** - NEW: Overflow detection tests
3. **test_cpt_math.cpp** - NEW: Math utility tests

## Test Strategy

### Principles
- **No external dependencies** - Tests run without GTS/Silo installed
- **Edge case focus** - Zero, overflow, NaN, boundary conditions
- **Fast execution** - All tests complete in <1 second
- **Descriptive assertions** - Clear test names and failure messages

### Coverage Goals
- All safety-critical functions: 100% branch coverage
- Math utilities: happy path + edge cases
- RAII: construction, destruction, move semantics

## Implementation Plan

### Phase 1: Enhance Existing Tests
- Expand test_sign.cpp with NaN and infinity cases

### Phase 2: Add Critical Safety Tests
- Create test_safe_mesh_size.cpp for overflow detection
- Add test for SignedDistance RAII behavior

### Phase 3: Add Math Utility Tests
- Create test_cpt_math.cpp for min/max/angle functions

## Success Criteria

- [ ] All tests pass with `make check`
- [ ] No undefined behavior (run with UBSAN)
- [ ] No memory leaks (run with ASan if possible)
- [ ] Total test execution < 1 second