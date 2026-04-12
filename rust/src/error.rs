//! Error types and FFI-compatible status codes.
//!
//! Each variant carries a fixed `i32` discriminant that is returned across the
//! C ABI boundary to the calling plugin layer.

use thiserror::Error;

/// Runtime error variants for the computational backend.
///
/// Represented as `#[repr(i32)]` so that the discriminant can be passed
/// directly through the C foreign-function interface.
#[derive(Error, Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum TropError {
    /// No error (code 0).
    #[error("Success")]
    Success = 0,

    /// A null pointer was passed across the FFI boundary (code 1).
    #[error("Null pointer")]
    NullPointer = 1,

    /// Matrix or vector dimensions are incompatible (code 2).
    #[error("Invalid dimension")]
    InvalidDimension = 2,

    /// The dataset contains no control (untreated) units (code 3).
    #[error("No control units")]
    NoControl = 3,

    /// The dataset contains no treated units (code 4).
    #[error("No treated units")]
    NoTreated = 4,

    /// An iterative solver did not converge within the allowed iterations (code 5).
    #[error("Convergence failure")]
    Convergence = 5,

    /// A matrix required for inversion or decomposition is singular (code 6).
    #[error("Singular matrix")]
    Singular = 6,

    /// Heap allocation failed (code 7).
    #[error("Memory allocation failure")]
    Memory = 7,

    /// An unrecoverable panic was caught at the FFI boundary (code 8).
    #[error("Rust panic")]
    RustPanic = 8,

    /// Leave-one-out cross-validation could not complete (code 9).
    #[error("LOOCV failure")]
    LoocvFail = 9,

    /// Bootstrap variance estimation failed (code 10).
    #[error("Bootstrap failure")]
    BootstrapFail = 10,

    /// Unclassified numerical or logical error (code 11).
    #[error("Computation failure")]
    Computation = 11,
}

impl TropError {
    /// Returns the `i32` discriminant for this variant.
    #[inline]
    pub fn code(&self) -> i32 {
        *self as i32
    }

    /// Maps an `i32` status code back to a variant.
    ///
    /// Returns `None` when `code` does not correspond to any defined variant.
    pub fn from_code(code: i32) -> Option<Self> {
        match code {
            0 => Some(TropError::Success),
            1 => Some(TropError::NullPointer),
            2 => Some(TropError::InvalidDimension),
            3 => Some(TropError::NoControl),
            4 => Some(TropError::NoTreated),
            5 => Some(TropError::Convergence),
            6 => Some(TropError::Singular),
            7 => Some(TropError::Memory),
            8 => Some(TropError::RustPanic),
            9 => Some(TropError::LoocvFail),
            10 => Some(TropError::BootstrapFail),
            11 => Some(TropError::Computation),
            _ => None,
        }
    }

    /// Returns `true` when the variant is [`Success`](TropError::Success).
    #[inline]
    pub fn is_success(&self) -> bool {
        matches!(self, TropError::Success)
    }
}

impl From<TropError> for i32 {
    #[inline]
    fn from(err: TropError) -> i32 {
        err.code()
    }
}

/// Convenience alias: `Result<T, TropError>`.
pub type TropResult<T> = Result<T, TropError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_codes() {
        assert_eq!(TropError::Success.code(), 0);
        assert_eq!(TropError::NullPointer.code(), 1);
        assert_eq!(TropError::InvalidDimension.code(), 2);
        assert_eq!(TropError::NoControl.code(), 3);
        assert_eq!(TropError::NoTreated.code(), 4);
        assert_eq!(TropError::Convergence.code(), 5);
        assert_eq!(TropError::Singular.code(), 6);
        assert_eq!(TropError::Memory.code(), 7);
        assert_eq!(TropError::RustPanic.code(), 8);
        assert_eq!(TropError::LoocvFail.code(), 9);
        assert_eq!(TropError::BootstrapFail.code(), 10);
        assert_eq!(TropError::Computation.code(), 11);
    }

    #[test]
    fn test_from_code() {
        assert_eq!(TropError::from_code(0), Some(TropError::Success));
        assert_eq!(TropError::from_code(5), Some(TropError::Convergence));
        assert_eq!(TropError::from_code(99), None);
    }

    #[test]
    fn test_into_i32() {
        let err: i32 = TropError::Convergence.into();
        assert_eq!(err, 5);
    }

    #[test]
    fn test_is_success() {
        assert!(TropError::Success.is_success());
        assert!(!TropError::Convergence.is_success());
    }
}
