# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.0] - 2026-04-12

### Added
- Initial release of TROP Stata package
- Twostep estimation method with triply robust properties
- Joint estimation method for simultaneous optimization
- Precompiled plugin for macOS ARM64; other platforms buildable from Rust source
- Bootstrap inference with parallel computation
- LOOCV cross-validation for regularization selection
- Comprehensive predict functionality (y0, y1, te, residuals, mu, alpha, beta)
- estat subcommands (summarize, weights, loocv, factors, bootstrap, sensitivity, vce)
- Full help documentation

### Fixed
- Numerical consistency verified against reference implementation

### Notes
- Requires Stata 17.0 or higher
- Rust-based computational backend for performance
