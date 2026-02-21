use bbrs_ecc::curves::bn254::Fr;

use super::params::{INTERNAL_DIAGONAL, ROUND_CONSTANTS};

const T: usize = 4;
const ROUNDS_F: usize = 8;
const ROUNDS_P: usize = 56;

/// S-box: x -> x^5
#[inline]
fn apply_single_sbox(x: &mut Fr) {
    let xx = x.sqr();
    let xxxx = xx.sqr();
    *x = *x * xxxx;
}

/// Apply S-box to all elements of state.
#[inline]
fn apply_sbox(state: &mut [Fr; T]) {
    for x in state.iter_mut() {
        apply_single_sbox(x);
    }
}

/// Add round constants to state.
#[inline]
fn add_round_constants(state: &mut [Fr; T], rc: &[Fr; T]) {
    for i in 0..T {
        state[i] = state[i] + rc[i];
    }
}

/// External (full) MDS matrix multiplication using the hardcoded 4x4 matrix:
///   | 5 7 1 3 |
///   | 4 6 1 1 |
///   | 1 3 5 7 |
///   | 1 1 4 6 |
///
/// Algorithm taken directly from the Poseidon2 paper.
#[inline]
fn matrix_multiplication_external(state: &mut [Fr; T]) {
    let t0 = state[0] + state[1]; // A + B
    let t1 = state[2] + state[3]; // C + D
    let mut t2 = state[1] + state[1]; // 2B
    t2 = t2 + t1; // 2B + C + D
    let mut t3 = state[3] + state[3]; // 2D
    t3 = t3 + t0; // 2D + A + B
    let mut t4 = t1 + t1;
    t4 = t4 + t4;
    t4 = t4 + t3; // A + B + 4C + 6D
    let mut t5 = t0 + t0;
    t5 = t5 + t5;
    t5 = t5 + t2; // 4A + 6B + C + D
    let t6 = t3 + t5; // 5A + 7B + C + 3D
    let t7 = t2 + t4; // A + 3B + 5C + 7D

    state[0] = t6;
    state[1] = t5;
    state[2] = t7;
    state[3] = t4;
}

/// Internal matrix multiplication: M_I * state where M_I has
/// diagonal = internal_matrix_diagonal and all off-diagonal = 1.
/// Equivalent to: state[i] = diagonal[i] * state[i] + sum(all states)
#[inline]
fn matrix_multiplication_internal(state: &mut [Fr; T]) {
    let diag = &*INTERNAL_DIAGONAL;
    let mut sum = state[0];
    for i in 1..T {
        sum = sum + state[i];
    }
    for i in 0..T {
        state[i] = state[i] * diag[i];
        state[i] = state[i] + sum;
    }
}

/// Poseidon2 permutation over BN254 Fr with t=4.
///
/// Structure: initial linear layer, then 4 full external rounds,
/// 56 partial internal rounds, 4 full external rounds.
pub fn permutation(input: &[Fr; T]) -> [Fr; T] {
    let rc = &*ROUND_CONSTANTS;
    let mut state = *input;

    // Initial linear layer
    matrix_multiplication_external(&mut state);

    // First set of external (full) rounds
    let rounds_f_beginning = ROUNDS_F / 2; // = 4
    for i in 0..rounds_f_beginning {
        add_round_constants(&mut state, &rc[i]);
        apply_sbox(&mut state);
        matrix_multiplication_external(&mut state);
    }

    // Internal (partial) rounds
    let p_end = rounds_f_beginning + ROUNDS_P; // = 60
    for i in rounds_f_beginning..p_end {
        state[0] = state[0] + rc[i][0];
        apply_single_sbox(&mut state[0]);
        matrix_multiplication_internal(&mut state);
    }

    // Remaining external (full) rounds
    let num_rounds = ROUNDS_F + ROUNDS_P; // = 64
    for i in p_end..num_rounds {
        add_round_constants(&mut state, &rc[i]);
        apply_sbox(&mut state);
        matrix_multiplication_external(&mut state);
    }

    state
}
