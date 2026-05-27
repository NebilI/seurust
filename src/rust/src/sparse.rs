//! dgCMatrix / dgRMatrix slot helpers.
use extendr_api::prelude::*;
use sprs::{CsMat, TriMat};

pub fn vec_from_doubles(x: &Doubles) -> Vec<f64> {
    x.iter().map(|v| v.0).collect()
}

pub fn vec_from_integers(i: &Integers) -> Vec<i32> {
    i.iter().map(|v| v.0).collect()
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
            x = Doubles::from_values(self.x.clone()),
            i = Integers::from_values(self.i.clone()),
            p = Integers::from_values(self.p.clone()),
            Dim = Integers::from_values(vec![self.nrows, self.ncols])
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
