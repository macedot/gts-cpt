# Test Makefile
CXX = g++
CXXFLAGS = -std=c++20 -Wall -Wextra -g -I../include
LDFLAGS = -lm

TESTS = test_sign test_safe_mesh_size test_cpt_math

all: $(TESTS)

test_sign: test_sign.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

test_safe_mesh_size: test_safe_mesh_size.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

test_cpt_math: test_cpt_math.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

clean:
	rm -f $(TESTS)

check: $(TESTS)
	@echo "Running tests..."
	@for test in $(TESTS); do \
		echo "Running $$test..."; \
		./$$test || exit 1; \
	done
	@echo "✅ All tests passed!"

.PHONY: all clean check