.PHONY: all check fmt fmt-check clippy test test-python build build-release clean prepush header install-python help

# Default target
all: check test

# Format code
fmt:
	cargo fmt --all

# Check formatting without modifying
fmt-check:
	cargo fmt --all -- --check

# Run clippy lints
clippy:
	PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 cargo clippy --workspace -- -D warnings

# Run all Rust tests
test:
	cargo test -p mordred-core -p mordred-ffi -p mordred-cli

# Run Python tests (requires maturin develop first)
test-python:
	cd mordred-py && maturin develop --release && pytest tests/

# Run all checks (formatting + clippy)
check: fmt-check clippy

# Build debug
build:
	PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 cargo build --workspace

# Build release
build-release:
	PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 cargo build --workspace --release

# Pre-push: run everything before pushing
prepush: fmt-check clippy test

# Clean build artifacts
clean:
	cargo clean

# Generate C header
header:
	cargo build -p mordred-ffi
	@echo "Header generated at mordred-ffi/include/mordred.h"

# Install Python package in development mode
install-python:
	cd mordred-py && maturin develop --release

# Show help
help:
	@echo "Available targets:"
	@echo "  all            - Run checks and tests (default)"
	@echo "  fmt            - Format code"
	@echo "  fmt-check      - Check formatting"
	@echo "  clippy         - Run clippy lints"
	@echo "  test           - Run Rust tests"
	@echo "  test-python    - Run Python tests"
	@echo "  check          - Run fmt-check + clippy"
	@echo "  build          - Build debug"
	@echo "  build-release  - Build release"
	@echo "  prepush        - Pre-push checks (fmt + clippy + test)"
	@echo "  clean          - Clean build artifacts"
	@echo "  header         - Generate C/C++ header"
	@echo "  install-python - Install Python package (dev mode)"
	@echo "  help           - Show this help"
