//! Optional authoritative native backends. Each instance owns its native handles;
//! no thread-sharing promises are made. See docs/native_backends.md for requirements.
#[cfg(feature = "gdal-backend")]
pub mod gdal;
#[cfg(feature = "proj-backend")]
pub mod proj;

use crate::{Error, Result};
use std::ffi::{CStr, CString, c_char, c_void};
use std::ptr::NonNull;

pub(crate) fn backend_error(backend: &str, message: impl std::fmt::Display) -> Error {
    Error::Backend {
        backend: backend.into(),
        message: message.to_string(),
    }
}
pub(crate) fn cstring(value: &str) -> Result<CString> {
    CString::new(value).map_err(|_| Error::InvalidInput("native argument contains NUL".into()))
}
/// Copy a nullable native string before its owning object is destroyed.
pub(crate) unsafe fn copy_string(pointer: *const c_char) -> String {
    if pointer.is_null() {
        String::new()
    } else {
        // SAFETY: caller guarantees a live native NUL-terminated string.
        unsafe { CStr::from_ptr(pointer) }
            .to_string_lossy()
            .into_owned()
    }
}
pub(crate) struct Handle {
    pub pointer: NonNull<c_void>,
    destroy: unsafe extern "C" fn(*mut c_void),
}
impl Handle {
    /// Caller transfers sole ownership and supplies the matching destructor.
    pub unsafe fn new(
        pointer: *mut c_void,
        destroy: unsafe extern "C" fn(*mut c_void),
        backend: &str,
        message: &str,
    ) -> Result<Self> {
        Ok(Self {
            pointer: NonNull::new(pointer).ok_or_else(|| backend_error(backend, message))?,
            destroy,
        })
    }
    pub fn ptr(&self) -> *mut c_void {
        self.pointer.as_ptr()
    }
}
impl Drop for Handle {
    fn drop(&mut self) {
        // SAFETY: Handle owns the live allocation; Drop runs exactly once.
        unsafe { (self.destroy)(self.ptr()) };
    }
}
