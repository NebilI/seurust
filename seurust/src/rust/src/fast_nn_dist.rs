use crate::utils::use_parallel;
use extendr_api::prelude::*;
use extendr_ffi::{
    R_NamesSymbol, Rboolean, Rf_allocVector, Rf_getAttrib, Rf_isNull, Rf_protect, Rf_setAttrib,
    Rf_unprotect, INTEGER, LENGTH, REAL, SET_VECTOR_ELT, SEXPTYPE, TYPEOF, VECTOR_ELT,
};
use rayon::prelude::*;

fn neighbor_id(raw: f64, cell: usize, nrows_y: usize) -> extendr_api::Result<usize> {
    if !raw.is_finite() || raw < 1.0 || raw > nrows_y as f64 {
        return Err(extendr_api::Error::Other(format!(
            "fast_dist: neighbor index {raw} for cell {} is outside 1..{nrows_y}.",
            cell + 1
        )));
    }
    Ok(raw as usize - 1)
}

fn parse_neighbor_ids(
    elt: extendr_ffi::SEXP,
    cell: usize,
    nrows_y: usize,
) -> extendr_api::Result<Vec<usize>> {
    unsafe {
        match TYPEOF(elt) {
            SEXPTYPE::INTSXP => {
                let n = LENGTH(elt) as usize;
                let ptr = INTEGER(elt);
                let mut ids = Vec::with_capacity(n);
                for j in 0..n {
                    let v = *ptr.add(j);
                    let raw = if v == i32::MIN { f64::NAN } else { v as f64 };
                    ids.push(neighbor_id(raw, cell, nrows_y)?);
                }
                Ok(ids)
            }
            SEXPTYPE::REALSXP => {
                let n = LENGTH(elt) as usize;
                let ptr = REAL(elt);
                let mut ids = Vec::with_capacity(n);
                for j in 0..n {
                    ids.push(neighbor_id(*ptr.add(j), cell, nrows_y)?);
                }
                Ok(ids)
            }
            _ => Err(extendr_api::Error::Other(format!(
                "fast_dist: neighbors for cell {} must be an integer or numeric vector.",
                cell + 1
            ))),
        }
    }
}

/// Euclidean distances from one row of `x` to the listed rows of `y`.
/// Four neighbors share a pass over the dimensions so the CPU has independent
/// add chains; each distance is still summed in column order.
pub(crate) fn fill_distances(
    out: &mut [f64],
    neighbors: &[usize],
    x: &[f64],
    y: &[f64],
    x_row: usize,
    ncols: usize,
    nrows_x: usize,
    nrows_y: usize,
) {
    debug_assert_eq!(out.len(), neighbors.len());
    let mut j = 0;
    while j + 4 <= neighbors.len() {
        let (a, b, c, d) = (
            neighbors[j],
            neighbors[j + 1],
            neighbors[j + 2],
            neighbors[j + 3],
        );
        let (mut s0, mut s1, mut s2, mut s3) = (0.0, 0.0, 0.0, 0.0);
        for col in 0..ncols {
            let xv = x[x_row + col * nrows_x];
            let base = col * nrows_y;
            let t0 = xv - y[a + base];
            let t1 = xv - y[b + base];
            let t2 = xv - y[c + base];
            let t3 = xv - y[d + base];
            s0 += t0 * t0;
            s1 += t1 * t1;
            s2 += t2 * t2;
            s3 += t3 * t3;
        }
        out[j] = s0.sqrt();
        out[j + 1] = s1.sqrt();
        out[j + 2] = s2.sqrt();
        out[j + 3] = s3.sqrt();
        j += 4;
    }
    for j in j..neighbors.len() {
        let mut sum = 0.0;
        let id = neighbors[j];
        for col in 0..ncols {
            let t = x[x_row + col * nrows_x] - y[id + col * nrows_y];
            sum += t * t;
        }
        out[j] = sum.sqrt();
    }
}

pub fn fast_dist_impl(x: &RMatrix<f64>, y: &RMatrix<f64>, n: &List) -> extendr_api::Result<Robj> {
    let ngraph_size = n.len();
    if x.nrows() != ngraph_size {
        return Ok(Robj::from(List::new(0)));
    }
    if x.ncols() != y.ncols() {
        return Err(extendr_api::Error::Other(format!(
            "fast_dist: x has {} columns but y has {}.",
            x.ncols(),
            y.ncols()
        )));
    }

    let ncols = x.ncols();
    let nrows_x = x.nrows();
    let nrows_y = y.nrows();
    let x_data = x.data();
    let y_data = y.data();
    let n_sexp = unsafe { n.get() };

    let mut neighbor_ids = Vec::with_capacity(ngraph_size);
    for i in 0..ngraph_size {
        let elt = unsafe { VECTOR_ELT(n_sexp, i as extendr_ffi::R_xlen_t) };
        neighbor_ids.push(parse_neighbor_ids(elt, i, nrows_y)?);
    }

    let work: usize = neighbor_ids
        .iter()
        .map(|ids| ids.len().saturating_mul(ncols))
        .sum();
    let parallel = use_parallel(ngraph_size, work);
    let precomputed = if parallel {
        let rows: Vec<Vec<f64>> = neighbor_ids
            .par_iter()
            .enumerate()
            .map(|(i, ids)| {
                let mut out = vec![0.0; ids.len()];
                fill_distances(&mut out, ids, x_data, y_data, i, ncols, nrows_x, nrows_y);
                out
            })
            .collect();
        Some(rows)
    } else {
        None
    };

    unsafe {
        let out = Rf_protect(Rf_allocVector(
            SEXPTYPE::VECSXP,
            ngraph_size as extendr_ffi::R_xlen_t,
        ));
        for i in 0..ngraph_size {
            let ids = &neighbor_ids[i];
            let elt = Rf_protect(Rf_allocVector(
                SEXPTYPE::REALSXP,
                ids.len() as extendr_ffi::R_xlen_t,
            ));
            let dest = if ids.is_empty() {
                &mut []
            } else {
                std::slice::from_raw_parts_mut(REAL(elt), ids.len())
            };
            if let Some(rows) = precomputed.as_ref() {
                dest.copy_from_slice(&rows[i]);
            } else {
                fill_distances(dest, ids, x_data, y_data, i, ncols, nrows_x, nrows_y);
            }
            SET_VECTOR_ELT(out, i as extendr_ffi::R_xlen_t, elt);
            Rf_unprotect(1);
        }
        let names = Rf_getAttrib(n_sexp, R_NamesSymbol);
        if Rf_isNull(names) != Rboolean::TRUE {
            Rf_setAttrib(out, R_NamesSymbol, names);
        }
        let robj = Robj::from_sexp(out);
        Rf_unprotect(1);
        Ok(robj)
    }
}

#[cfg(test)]
mod tests {
    use super::fill_distances;

    fn naive(
        neighbors: &[usize],
        x: &[f64],
        y: &[f64],
        x_row: usize,
        ncols: usize,
        nrows_x: usize,
        nrows_y: usize,
    ) -> Vec<f64> {
        neighbors
            .iter()
            .map(|&id| {
                let mut sum = 0.0;
                for col in 0..ncols {
                    let t = x[x_row + col * nrows_x] - y[id + col * nrows_y];
                    sum += t * t;
                }
                sum.sqrt()
            })
            .collect()
    }

    #[test]
    fn four_neighbor_pass_is_bit_identical() {
        let nrows = 6;
        let ncols = 5;
        let x: Vec<f64> = (0..nrows * ncols).map(|i| i as f64 * 0.1).collect();
        let y: Vec<f64> = (0..nrows * ncols)
            .map(|i| (i as f64) * 0.07 - 0.2)
            .collect();
        let neighbors = vec![0, 2, 4, 1, 3];
        let expected = naive(&neighbors, &x, &y, 2, ncols, nrows, nrows);
        let mut got = vec![0.0; neighbors.len()];
        fill_distances(&mut got, &neighbors, &x, &y, 2, ncols, nrows, nrows);
        for (a, b) in got.iter().zip(&expected) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }
}
