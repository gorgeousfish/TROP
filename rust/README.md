# TROP Core - Rust Implementation

[![Rust Tests](https://github.com/gorgeousfish/TROP/actions/workflows/rust-tests.yml/badge.svg)](https://github.com/gorgeousfish/TROP/actions/workflows/rust-tests.yml)

Pure Rust implementation of the TROP (Triply Robust Panel) estimator core algorithms, designed for integration with Stata via C ABI.

## Overview

This library implements the algorithms from:

> Athey, Imbens, Qu & Viviano (2025) "Triply Robust Panel Estimators"

## Features

- **Distance computation**: Unit and time distance matrices
- **Weight computation**: Twostep and Joint weighting schemes
- **LOOCV search**: Cross-validation for hyperparameter selection
- **Model estimation**: Alternating minimization with SVD-based nuclear norm regularization
- **Bootstrap inference**: Unit-level block bootstrap with stratified sampling

## Building

```bash
# Debug build
cargo build

# Release build (optimized)
cargo build --release
```

## Testing

```bash
# Run all tests
cargo test --release -- --test-threads=1

# Run with output
cargo test -- --nocapture

# Run specific module tests
cargo test bootstrap::
cargo test estimation::
cargo test weights::
```

### Test Coverage

| Module | Tests | Tolerance |
|--------|-------|-----------|
| distance.rs | 16 | < 1e-12 |
| weights.rs | 13 | < 1e-10 |
| estimation.rs | 17 | < 1e-10 (SVD), < 1e-6 (ATT) |
| loocv.rs | 12 | < 1e-8 |
| bootstrap.rs | 20 | < 1e-4 (SE), < 1e-12 (variance) |

Total: 110+ tests covering unit tests, integration tests, and property-based tests.

## Dependencies

- `ndarray` - N-dimensional arrays
- `faer` - Pure Rust linear algebra (SVD)
- `rayon` - Parallel iteration
- `rand` / `rand_xoshiro` - Random number generation
- `thiserror` - Error handling

## C ABI Exports

The library exports C-compatible functions for Stata integration:

- `stata_loocv_grid_search` - LOOCV hyperparameter search (Twostep)
- `stata_estimate_twostep` - Twostep model estimation
- `stata_bootstrap_trop_variance` - Bootstrap variance (Twostep)
- `stata_loocv_grid_search_joint` - LOOCV search (Joint)
- `stata_estimate_joint` - Joint model estimation
- `stata_bootstrap_trop_variance_joint` - Bootstrap variance (Joint)
- `stata_compute_unit_distance_matrix` - Distance matrix computation

## Module Structure

```
src/
├── lib.rs          # C ABI exports and panic protection
├── error.rs        # Error types (TropError enum)
├── distance.rs     # Distance computation
├── weights.rs      # Weight computation
├── estimation.rs   # Model estimation with SVD
├── loocv.rs        # LOOCV search
└── bootstrap.rs    # Bootstrap variance estimation
```

## License

MIT
