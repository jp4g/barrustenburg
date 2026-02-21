use crate::ecc::fields::field::Field;
use crate::ecc::fields::field_params::FieldParams;

/// A compact array with three logical zones:
///   [0, start)           -> virtual zeros
///   [start, end)         -> real data backed by Vec
///   [end, virtual_size)  -> virtual zeros
///
/// Rust port of Barretenberg's `SharedShiftedVirtualZeroesArray`.
pub struct SharedShiftedVirtualZeroesArray<P: FieldParams> {
    data: Vec<Field<P>>,
    start: usize,
    end: usize,
    virtual_size: usize,
}

// Manual Clone because derive requires P: Clone unnecessarily (Field<P> is Copy).
impl<P: FieldParams> Clone for SharedShiftedVirtualZeroesArray<P> {
    fn clone(&self) -> Self {
        Self {
            data: self.data.clone(),
            start: self.start,
            end: self.end,
            virtual_size: self.virtual_size,
        }
    }
}

impl<P: FieldParams> SharedShiftedVirtualZeroesArray<P> {
    /// Create an array of `size` zeros starting at `start_index`.
    pub fn new(size: usize, virtual_size: usize, start_index: usize) -> Self {
        assert!(
            start_index + size <= virtual_size,
            "start_index ({}) + size ({}) exceeds virtual_size ({})",
            start_index,
            size,
            virtual_size,
        );
        Self {
            data: vec![Field::zero(); size],
            start: start_index,
            end: start_index + size,
            virtual_size,
        }
    }

    /// Wrap existing data starting at `start_index`.
    pub fn from_vec(data: Vec<Field<P>>, virtual_size: usize, start_index: usize) -> Self {
        let end = start_index + data.len();
        assert!(
            end <= virtual_size,
            "start_index ({}) + data.len() ({}) exceeds virtual_size ({})",
            start_index,
            data.len(),
            virtual_size,
        );
        Self {
            data,
            start: start_index,
            end,
            virtual_size,
        }
    }

    /// Read element at logical `index`. Returns zero for virtual zones.
    #[inline]
    pub fn get(&self, index: usize) -> Field<P> {
        if index >= self.start && index < self.end {
            self.data[index - self.start]
        } else {
            Field::zero()
        }
    }

    /// Write element at logical `index`. Panics if outside real data range.
    #[inline]
    pub fn set(&mut self, index: usize, value: Field<P>) {
        assert!(
            index >= self.start && index < self.end,
            "index {} out of real range [{}..{})",
            index,
            self.start,
            self.end,
        );
        self.data[index - self.start] = value;
    }

    /// Backing data slice.
    #[inline]
    pub fn data(&self) -> &[Field<P>] {
        &self.data
    }

    /// Mutable backing data slice.
    #[inline]
    pub fn data_mut(&mut self) -> &mut [Field<P>] {
        &mut self.data
    }

    /// Length of real data (end - start).
    #[inline]
    pub fn size(&self) -> usize {
        self.end - self.start
    }

    /// Total logical size including virtual zeros.
    #[inline]
    pub fn virtual_size(&self) -> usize {
        self.virtual_size
    }

    /// First real index.
    #[inline]
    pub fn start_index(&self) -> usize {
        self.start
    }

    /// One past the last real index.
    #[inline]
    pub fn end_index(&self) -> usize {
        self.end
    }
}
