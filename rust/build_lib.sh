#!/bin/bash
# Build Rust library for TROP Stata plugin

cd "$(dirname "$0")"
cargo build --release

# Copy the library to plugin directory
cp target/release/libtrop_core.dylib ../plugin/

echo "Build complete!"
