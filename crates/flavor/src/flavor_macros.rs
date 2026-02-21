//! Port of `flavor_macros.hpp` — macros for defining flavor entity structs.
//!
//! The C++ `DEFINE_FLAVOR_MEMBERS` creates named fields with `get_all()` and `get_labels()`.
//! In Rust, we use declarative macros to generate the same pattern.

/// Define a flavor entity struct with named fields and iteration support.
///
/// Generates:
/// - A struct with the given name and named fields of type `$data_type`
/// - `get_all(&self) -> [&$data_type; N]` for read access
/// - `get_all_mut(&mut self) -> [&mut $data_type; N]` for write access
/// - `get_labels() -> &'static [&'static str]` for field names
/// - `size() -> usize` for the number of fields
///
/// # Example
/// ```ignore
/// define_flavor_members! {
///     pub struct PrecomputedEntities<T> {
///         q_m: T,
///         q_l: T,
///         q_r: T,
///     }
/// }
/// ```
#[macro_export]
macro_rules! define_flavor_members {
    (
        $(#[$meta:meta])*
        $vis:vis struct $name:ident<$T:ident> {
            $(
                $(#[$field_meta:meta])*
                $field_vis:vis $field:ident : $field_ty:ty
            ),* $(,)?
        }
    ) => {
        $(#[$meta])*
        $vis struct $name<$T> {
            $(
                $(#[$field_meta])*
                $field_vis $field : $field_ty,
            )*
        }

        impl<$T> $name<$T> {
            /// Number of fields in this entity struct.
            pub const fn size() -> usize {
                $crate::flavor_macros::count!($($field),*)
            }

            /// Get field labels (matching C++ `get_labels()`).
            pub fn get_labels() -> &'static [&'static str] {
                &[$(stringify!($field)),*]
            }
        }

        impl<$T: Copy> $name<$T> {
            /// Get a read-only array of all fields.
            pub fn get_all(&self) -> [&$T; $crate::flavor_macros::count!($($field),*)] {
                [$(&self.$field),*]
            }

            /// Get a mutable array of all fields.
            pub fn get_all_mut(&mut self) -> [&mut $T; $crate::flavor_macros::count!($($field),*)] {
                [$(&mut self.$field),*]
            }
        }
    };
}

/// Count the number of identifiers in a list.
#[macro_export]
macro_rules! _count {
    () => { 0usize };
    ($head:ident $(, $tail:ident)*) => { 1usize + $crate::flavor_macros::count!($($tail),*) };
}

// Re-export the macros at module level
pub use _count as count;
pub use define_flavor_members;
