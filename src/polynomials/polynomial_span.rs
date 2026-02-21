use crate::ecc::fields::field::Field;
use crate::ecc::fields::field_params::FieldParams;

pub struct PolynomialSpan<'a, P: FieldParams> {
    pub start_index: usize,
    pub span: &'a [Field<P>],
}

impl<'a, P: FieldParams> Clone for PolynomialSpan<'a, P> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, P: FieldParams> Copy for PolynomialSpan<'a, P> {}

impl<'a, P: FieldParams> PolynomialSpan<'a, P> {
    pub fn new(span: &'a [Field<P>], start_index: usize) -> Self {
        Self { start_index, span }
    }

    pub fn end_index(&self) -> usize {
        self.start_index + self.span.len()
    }

    pub fn size(&self) -> usize {
        self.span.len()
    }

    pub fn get(&self, index: usize) -> Field<P> {
        if index >= self.start_index && index < self.end_index() {
            self.span[index - self.start_index]
        } else {
            Field::zero()
        }
    }

    pub fn subspan(&self, offset: usize, length: usize) -> Self {
        Self {
            start_index: self.start_index + offset,
            span: &self.span[offset..offset + length],
        }
    }
}
