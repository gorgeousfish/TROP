# Makefile for TROP Stata Plugin
# 
# This Makefile builds the C bridge layer that connects Stata to the Rust core library.
# Supports macOS (ARM64/x64) and Linux (x64).
#
# Usage:
#   make          - Build the plugin
#   make clean    - Remove build artifacts
#   make test     - Run basic tests (requires Stata)
#   make rust     - Build the Rust library first
#   make all      - Build Rust library and plugin

# ============================================================================
# Platform Detection
# ============================================================================

UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

# ============================================================================
# Rust Library Configuration
# ============================================================================

RUST_DIR := ../rust
RUST_TARGET_DIR := $(RUST_DIR)/target/release
RUST_LIB := $(RUST_TARGET_DIR)/libtrop_core.a

# ============================================================================
# Platform-Specific Configuration
# ============================================================================

ifeq ($(UNAME_S),Darwin)
    # macOS
    CC := clang
    # Define SYSTEM=APPLEMAC (3) for macOS
    SYSTEM_DEF := -DSYSTEM=3
    
    ifeq ($(UNAME_M),arm64)
        # Apple Silicon (M1/M2/M3)
        TARGET := trop_macos_arm64.plugin
        ARCH_FLAGS := -arch arm64
        RUST_TARGET := aarch64-apple-darwin
    else
        # Intel Mac
        TARGET := trop_macos_x64.plugin
        ARCH_FLAGS := -arch x86_64
        RUST_TARGET := x86_64-apple-darwin
    endif
    
    # macOS linker flags for plugin bundle
    LDFLAGS := -bundle -undefined dynamic_lookup $(ARCH_FLAGS)
    
    # macOS needs to link against system frameworks for Rust
    LIBS := -lSystem -lc -lm
    
else ifeq ($(UNAME_S),Linux)
    # Linux
    CC := gcc
    # Define SYSTEM=OPUNIX (2) for Linux
    SYSTEM_DEF := -DSYSTEM=2
    TARGET := trop_linux_x64.plugin
    ARCH_FLAGS :=
    RUST_TARGET := x86_64-unknown-linux-gnu
    
    # Linux linker flags for shared library
    LDFLAGS := -shared -fPIC
    LIBS := -lm -lpthread -ldl
    
else
    $(error Unsupported platform: $(UNAME_S))
endif

# ============================================================================
# Compiler Flags
# ============================================================================

# Common flags
CFLAGS := -O2 -fPIC -Wall -Wextra $(ARCH_FLAGS) $(SYSTEM_DEF)

# Debug build (use: make DEBUG=1)
ifdef DEBUG
    CFLAGS := -g -O0 -fPIC -Wall -Wextra -DDEBUG $(ARCH_FLAGS) $(SYSTEM_DEF)
endif

# Include paths
INCLUDES := -I.

# ============================================================================
# Source Files
# ============================================================================

SRCS := stata_bridge.c stplugin.c
OBJS := $(SRCS:.c=.o)
HEADERS := stata_bridge.h stplugin.h

# ============================================================================
# Build Rules
# ============================================================================

.PHONY: all clean rust test info

# Default target: build plugin
all: rust $(TARGET)

# Build plugin only (assumes Rust library exists)
$(TARGET): $(OBJS) $(RUST_LIB)
	@echo "Linking $(TARGET)..."
	$(CC) $(LDFLAGS) -o $@ $(OBJS) $(RUST_LIB) $(LIBS)
	@echo "Build complete: $(TARGET)"

# Compile C source files
%.o: %.c $(HEADERS)
	@echo "Compiling $<..."
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Build Rust library
rust:
	@echo "Building Rust library..."
	cd $(RUST_DIR) && cargo build --release
	@echo "Rust library built: $(RUST_LIB)"

# Clean build artifacts
clean:
	@echo "Cleaning..."
	rm -f $(OBJS) $(TARGET)
	rm -f *.plugin
	@echo "Clean complete"

# Deep clean (including Rust)
distclean: clean
	cd $(RUST_DIR) && cargo clean

# Show build information
info:
	@echo "Platform: $(UNAME_S) $(UNAME_M)"
	@echo "Target: $(TARGET)"
	@echo "Compiler: $(CC)"
	@echo "CFLAGS: $(CFLAGS)"
	@echo "LDFLAGS: $(LDFLAGS)"
	@echo "Rust library: $(RUST_LIB)"

# Check if Rust library exists
$(RUST_LIB):
	@echo "Rust library not found. Building..."
	$(MAKE) rust

# ============================================================================
# Installation
# ============================================================================

# Install to Stata PERSONAL directory
install: $(TARGET)
	@echo "Installing $(TARGET) to Stata PERSONAL directory..."
	@if [ -z "$(STATA_PERSONAL)" ]; then \
		echo "Error: STATA_PERSONAL environment variable not set"; \
		echo "Set it to your Stata PERSONAL directory, e.g.:"; \
		echo "  export STATA_PERSONAL=~/Library/Application Support/Stata/ado/personal"; \
		exit 1; \
	fi
	mkdir -p "$(STATA_PERSONAL)/t"
	cp $(TARGET) "$(STATA_PERSONAL)/t/trop.plugin"
	@echo "Installed to $(STATA_PERSONAL)/t/trop.plugin"

# ============================================================================
# Testing
# ============================================================================

# Basic test (requires Stata)
test: $(TARGET)
	@echo "Running basic plugin test..."
	@echo "Note: This requires Stata to be installed and in PATH"
	@if command -v stata-mp >/dev/null 2>&1; then \
		stata-mp -b do test_plugin.do; \
	elif command -v stata-se >/dev/null 2>&1; then \
		stata-se -b do test_plugin.do; \
	elif command -v stata >/dev/null 2>&1; then \
		stata -b do test_plugin.do; \
	else \
		echo "Stata not found in PATH. Skipping test."; \
	fi

# ============================================================================
# Development Helpers
# ============================================================================

# Watch for changes and rebuild (requires fswatch on macOS)
watch:
	@echo "Watching for changes..."
	fswatch -o $(SRCS) $(HEADERS) | xargs -n1 -I{} make $(TARGET)

# Format code (requires clang-format)
format:
	clang-format -i $(SRCS) $(HEADERS)

# Static analysis (requires clang)
analyze:
	clang --analyze $(CFLAGS) $(INCLUDES) $(SRCS)
