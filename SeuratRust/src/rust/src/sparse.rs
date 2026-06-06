//! dgCMatrix / dgRMatrix slot helpers.
use extendr_api::prelude::*;
use sprs::{CsMat, TriMat};

pub fn vec_from_doubles(x: &Doubles) -> Vec<f64> {
    x.iter().map(|v| v.0).collect()
}

pub fn vec_from_integers(i: &Integers) -> Vec<i32> {
    i.iter().map(|v| v.0).collect()
}

fn vec_to_doubles(values: Vec<f64>) -> Doubles {
    let len = values.len();
    let mut out = Doubles::new(len);
    if len > 0 {
        out.as_robj_mut()
            .as_real_slice_mut()
            .expect("numeric output")
            .copy_from_slice(&values);
    }
    out
}

fn vec_to_integers(values: Vec<i32>) -> Integers {
    let len = values.len();
    let mut out = Integers::new(len);
    if len > 0 {
        out.as_robj_mut()
            .as_integer_slice_mut()
            .expect("integer output")
            .copy_from_slice(&values);
    }
    out
}

/// Build a dgCMatrix from preallocated slot vectors (no extra copy).
pub fn dgcmatrix_from_buffers(
    x: Doubles,
    i: Integers,
    p: Integers,
    dim: Integers,
) -> extendr_api::Result<Robj> {
    call!(
        "methods::new",
        "dgCMatrix",
        x = x,
        i = i,
        p = p,
        Dim = dim
    )
}

/// Sort/merge coordinate triplets and write CSC slots directly into R memory.
pub fn dgcmatrix_from_triplets(
    nrows: i32,
    ncols: i32,
    mut triplets: Vec<(usize, usize, f64)>,
) -> extendr_api::Result<Robj> {
    let ncols_usize = ncols as usize;
    let dim = Integers::from_values(vec![nrows, ncols]);

    if triplets.is_empty() {
        let mut p = Integers::new(ncols_usize + 1);
        p.as_robj_mut()
            .as_integer_slice_mut()
            .expect("p")
            .fill(0);
        return dgcmatrix_from_buffers(Doubles::new(0), Integers::new(0), p, dim);
    }

    triplets.sort_unstable_by_key(|&(r, c, _)| (c, r));

    let mut merged: Vec<(usize, usize, f64)> = Vec::with_capacity(triplets.len());
    for (r, c, v) in triplets {
        if let Some(last) = merged.last_mut() {
            if last.0 == r && last.1 == c {
                last.2 += v;
                continue;
            }
        }
        merged.push((r, c, v));
    }

    let nnz = merged.len();
    let mut x_out = Doubles::new(nnz);
    let mut i_out = Integers::new(nnz);
    let mut p_out = Integers::new(ncols_usize + 1);

    let x = x_out
        .as_robj_mut()
        .as_real_slice_mut()
        .expect("numeric x");
    let i = i_out
        .as_robj_mut()
        .as_integer_slice_mut()
        .expect("integer i");
    let p = p_out
        .as_robj_mut()
        .as_integer_slice_mut()
        .expect("integer p");

    let mut nz = 0usize;
    for col in 0..ncols_usize {
        p[col] = nz as i32;
        while nz < nnz && merged[nz].1 == col {
            i[nz] = merged[nz].0 as i32;
            x[nz] = merged[nz].2;
            nz += 1;
        }
    }
    p[ncols_usize] = nnz as i32;

    dgcmatrix_from_buffers(x_out, i_out, p_out, dim)
}

/// Borrowed dgCMatrix CSC slots backed by R memory (zero-copy input).
#[derive(Clone, Copy, Debug)]
pub struct CscView<'a> {
    pub x: &'a [f64],
    pub i: &'a [i32],
    pub p: &'a [i32],
    pub nrows: i32,
    pub ncols: i32,
}

impl<'a> CscView<'a> {
    pub fn from_slots(
        x: &'a Doubles,
        i: &'a Integers,
        p: &'a Integers,
        nrows: i32,
        ncols: i32,
    ) -> Self {
        Self {
            x: x.as_robj().as_real_slice().expect("numeric x"),
            i: i.as_robj().as_integer_slice().expect("integer i"),
            p: p.as_robj().as_integer_slice().expect("integer p"),
            nrows,
            ncols,
        }
    }

    pub fn col_sums(&self) -> Vec<f64> {
        let ncols = self.ncols as usize;
        let mut sums = vec![0.0; ncols];
        for col in 0..ncols {
            for idx in self.p[col] as usize..self.p[col + 1] as usize {
                sums[col] += self.x[idx];
            }
        }
        sums
    }
}

/// Column-compressed sparse matrix slots (dgCMatrix).
#[derive(Clone, Debug)]
pub struct CscSlots {
    pub x: Vec<f64>,
    pub i: Vec<i32>,
    pub p: Vec<i32>,
    pub nrows: i32,
    pub ncols: i32,
}

impl CscSlots {
    pub fn from_r(x: Doubles, i: Integers, p: Integers, nrows: i32, ncols: i32) -> Self {
        Self {
            x: vec_from_doubles(&x),
            i: vec_from_integers(&i),
            p: vec_from_integers(&p),
            nrows,
            ncols,
        }
    }

    pub fn to_r_list(&self) -> List {
        list!(
            x = vec_to_doubles(self.x.clone()),
            i = vec_to_integers(self.i.clone()),
            p = vec_to_integers(self.p.clone()),
            Dim = Integers::from_values(vec![self.nrows, self.ncols])
        )
    }

    pub fn into_r_list(self) -> List {
        list!(
            x = vec_to_doubles(self.x),
            i = vec_to_integers(self.i),
            p = vec_to_integers(self.p),
            Dim = Integers::from_values(vec![self.nrows, self.ncols])
        )
    }

    /// Build a dgCMatrix SEXP directly (avoids R-side `Matrix::sparseMatrix`).
    pub fn into_r_dgcmatrix(self) -> extendr_api::Result<Robj> {
        dgcmatrix_from_buffers(
            vec_to_doubles(self.x),
            vec_to_integers(self.i),
            vec_to_integers(self.p),
            Integers::from_values(vec![self.nrows, self.ncols]),
        )
    }

    pub fn col_sums(&self) -> Vec<f64> {
        let ncols = self.ncols as usize;
        let mut sums = vec![0.0; ncols];
        for col in 0..ncols {
            for idx in self.p[col] as usize..self.p[col + 1] as usize {
                sums[col] += self.x[idx];
            }
        }
        sums
    }

    pub fn get(&self, row: usize, col: usize) -> f64 {
        for idx in self.p[col] as usize..self.p[col + 1] as usize {
            if self.i[idx] as usize == row {
                return self.x[idx];
            }
        }
        0.0
    }

    pub fn to_cs_mat(&self) -> CsMat<f64> {
        let shape = (self.nrows as usize, self.ncols as usize);
        let indptr: Vec<usize> = self.p.iter().map(|&v| v as usize).collect();
        let indices: Vec<usize> = self.i.iter().map(|&v| v as usize).collect();
        CsMat::new_csc(shape, indptr, indices, self.x.clone())
    }

    pub fn from_cs_mat(mat: &CsMat<f64>) -> Self {
        let (nrows, ncols) = mat.shape();
        let indptr = mat.indptr();
        let ip = indptr.raw_storage();
        let p: Vec<i32> = (0..=ncols).map(|col| ip[col] as i32).collect();
        let i: Vec<i32> = mat.indices().iter().map(|&v| v as i32).collect();
        let x = mat.data().to_vec();
        Self {
            x,
            i,
            p,
            nrows: nrows as i32,
            ncols: ncols as i32,
        }
    }
}

/// Row-compressed sparse matrix slots (dgRMatrix).
#[derive(Clone, Debug)]
pub struct CsrSlots {
    pub x: Vec<f64>,
    pub j: Vec<i32>,
    pub p: Vec<i32>,
    pub nrows: i32,
    pub ncols: i32,
}

impl CsrSlots {
    pub fn from_r(x: Doubles, j: Integers, p: Integers, nrows: i32, ncols: i32) -> Self {
        Self {
            x: vec_from_doubles(&x),
            j: vec_from_integers(&j),
            p: vec_from_integers(&p),
            nrows,
            ncols,
        }
    }

    pub fn to_cs_mat(&self) -> CsMat<f64> {
        let mut tri = TriMat::new((self.nrows as usize, self.ncols as usize));
        for row in 0..self.nrows as usize {
            for idx in self.p[row] as usize..self.p[row + 1] as usize {
                tri.add_triplet(row, self.j[idx] as usize, self.x[idx]);
            }
        }
        tri.to_csr()
    }
}

pub fn csc_from_triplets(
    nrows: usize,
    ncols: usize,
    triplets: &[(usize, usize, f64)],
) -> CscSlots {
    let mut tri = TriMat::new((nrows, ncols));
    for &(r, c, v) in triplets {
        tri.add_triplet(r, c, v);
    }
    CscSlots::from_cs_mat(&tri.to_csc())
}

/// Build dgCMatrix slots from coordinate triplets without sprs conversion.
pub fn csc_slots_from_triplets(
    nrows: i32,
    ncols: i32,
    mut triplets: Vec<(usize, usize, f64)>,
) -> CscSlots {
    let ncols_usize = ncols as usize;
    if triplets.is_empty() {
        return CscSlots {
            x: Vec::new(),
            i: Vec::new(),
            p: vec![0; ncols_usize + 1],
            nrows,
            ncols,
        };
    }

    triplets.sort_unstable_by_key(|&(r, c, _)| (c, r));

    // Sum duplicate (row, col) entries — matches Eigen::setFromTriplets behavior.
    let mut merged: Vec<(usize, usize, f64)> = Vec::with_capacity(triplets.len());
    for (r, c, v) in triplets {
        if let Some(last) = merged.last_mut() {
            if last.0 == r && last.1 == c {
                last.2 += v;
                continue;
            }
        }
        merged.push((r, c, v));
    }
    let triplets = merged;

    let nnz = triplets.len();
    let mut p = vec![0i32; ncols_usize + 1];
    let mut i = Vec::with_capacity(nnz);
    let mut x = Vec::with_capacity(nnz);

    let mut nz = 0;
    for col in 0..ncols_usize {
        p[col] = nz as i32;
        while nz < nnz && triplets[nz].1 == col {
            i.push(triplets[nz].0 as i32);
            x.push(triplets[nz].2);
            nz += 1;
        }
    }
    p[ncols_usize] = nnz as i32;

    CscSlots { x, i, p, nrows, ncols }
}

pub fn rmatrix_from_ndarray(values: ndarray::ArrayView2<f64>) -> RMatrix<f64> {
    let (nrows, ncols) = values.dim();
    RMatrix::new_matrix(nrows, ncols, |r, c| values[[r, c]])
}

pub fn ndarray_from_rmatrix(mat: &RMatrix<f64>) -> ndarray::Array2<f64> {
    let mut values = ndarray::Array2::zeros((mat.nrows(), mat.ncols()));
    for r in 0..mat.nrows() {
        for c in 0..mat.ncols() {
            values[[r, c]] = mat[[r, c]];
        }
    }
    values
}

pub fn strings_to_str_vec(names: Strings) -> Vec<String> {
    names.into_iter().map(|s| s.to_string()).collect()
}
