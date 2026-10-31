use rayon::prelude::*;
use std::time::Instant;

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn unif(&mut self) -> f64 {
        (self.next() as f64) / ((1u64 << 31) as f64)
    }
}

struct Csc {
    x: Vec<f64>,
    i: Vec<i32>,
    p: Vec<i32>,
    nrows: usize,
    ncols: usize,
}

fn sparse(nrows: usize, ncols: usize, density: f64, seed: u64) -> Csc {
    let mut r = Lcg(seed);
    let mut x = Vec::new();
    let mut i = Vec::new();
    let mut p = vec![0i32];
    for _ in 0..ncols {
        for row in 0..nrows {
            if r.unif() < density {
                i.push(row as i32);
                x.push(1.0 + (r.next() % 9) as f64);
            }
        }
        p.push(x.len() as i32);
    }
    Csc { x, i, p, nrows, ncols }
}

/// Fresh, untouched allocation, like R's allocVector for large vectors.
fn alloc_uninit(n: usize) -> Vec<f64> {
    let mut v = Vec::with_capacity(n);
    unsafe { v.set_len(n) };
    v
}

fn median_ms(reps: usize, mut f: impl FnMut() -> f64) -> f64 {
    let mut sink = 0.0;
    sink += f();
    let mut t: Vec<f64> = (0..reps)
        .map(|_| {
            let s = Instant::now();
            sink += f();
            s.elapsed().as_secs_f64() * 1e3
        })
        .collect();
    t.sort_by(|a, b| a.partial_cmp(b).unwrap());
    std::hint::black_box(sink);
    t[reps / 2]
}

#[derive(Clone, Copy)]
struct SendPtr(*mut f64);
unsafe impl Send for SendPtr {}
unsafe impl Sync for SendPtr {}

// ---------------- LogNorm ----------------
fn lognorm_current(m: &Csc) -> f64 {
    let mut sums = vec![0.0; m.ncols];
    for c in 0..m.ncols {
        for k in m.p[c] as usize..m.p[c + 1] as usize {
            sums[c] += m.x[k];
        }
    }
    let mut out = alloc_uninit(m.x.len());
    for v in out.iter_mut() {
        *v = 0.0;
    }
    out.copy_from_slice(&m.x);
    for c in 0..m.ncols {
        for k in m.p[c] as usize..m.p[c + 1] as usize {
            out[k] = (out[k] / sums[c] * 1e4).ln_1p();
        }
    }
    out[out.len() / 2]
}

fn lognorm_fused_par(m: &Csc) -> f64 {
    let mut out = alloc_uninit(m.x.len());
    let ptr = SendPtr(out.as_mut_ptr());
    (0..m.ncols).into_par_iter().with_min_len(64).for_each(|c| {
        let ptr = ptr;
        let s = m.p[c] as usize;
        let e = m.p[c + 1] as usize;
        let col = &m.x[s..e];
        let sum: f64 = col.iter().sum();
        for (k, &v) in col.iter().enumerate() {
            unsafe { *ptr.0.add(s + k) = (v / sum * 1e4).ln_1p() };
        }
    });
    out[out.len() / 2]
}

// ---------------- FastSparseRowScale ----------------
fn row_stats(m: &Csc) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = m.nrows;
    let mut s = vec![0.0; n];
    let mut sq = vec![0.0; n];
    for k in 0..m.x.len() {
        let r = m.i[k] as usize;
        s[r] += m.x[k];
        sq[r] += m.x[k] * m.x[k];
    }
    let cols = m.ncols as f64;
    let mut mu = vec![0.0; n];
    let mut inv = vec![0.0; n];
    let mut zero = vec![0.0; n];
    for r in 0..n {
        mu[r] = s[r] / cols;
        let sd = ((sq[r] - s[r] * s[r] / cols) / (cols - 1.0)).sqrt();
        inv[r] = 1.0 / sd;
        zero[r] = ((0.0 - mu[r]) * inv[r]).min(10.0);
    }
    (mu, inv, zero)
}

fn fill_col(m: &Csc, col: usize, out_col: &mut [f64], mu: &[f64], inv: &[f64], zero: &[f64]) {
    out_col.copy_from_slice(zero);
    for k in m.p[col] as usize..m.p[col + 1] as usize {
        let r = m.i[k] as usize;
        out_col[r] = ((m.x[k] - mu[r]) * inv[r]).min(10.0);
    }
}

fn rowscale_current(m: &Csc) -> f64 {
    let (mu, inv, zero) = row_stats(m);
    let mut data = vec![0.0; m.nrows * m.ncols];
    data.par_chunks_mut(m.nrows)
        .enumerate()
        .for_each(|(c, oc)| fill_col(m, c, oc, &mu, &inv, &zero));
    let mut out = alloc_uninit(data.len());
    out.copy_from_slice(&data);
    out[out.len() / 3]
}

fn rowscale_direct(m: &Csc) -> f64 {
    let (mu, inv, zero) = row_stats(m);
    let mut out = alloc_uninit(m.nrows * m.ncols);
    out.par_chunks_mut(m.nrows)
        .enumerate()
        .for_each(|(c, oc)| fill_col(m, c, oc, &mu, &inv, &zero));
    out[out.len() / 3]
}

// Known stats: current = RowIndex + gene-parallel strided writes + memcpy.
fn knownstats_current(m: &Csc, mu: &[f64], sd: &[f64]) -> f64 {
    let n = m.nrows;
    let mut cnt = vec![0usize; n];
    for &r in &m.i {
        cnt[r as usize] += 1;
    }
    let mut ptr = vec![0usize; n + 1];
    for r in 0..n {
        ptr[r + 1] = ptr[r] + cnt[r];
    }
    let mut cols = vec![0usize; ptr[n]];
    let mut xi = vec![0usize; ptr[n]];
    let mut cur = ptr.clone();
    for c in 0..m.ncols {
        for k in m.p[c] as usize..m.p[c + 1] as usize {
            let r = m.i[k] as usize;
            cols[cur[r]] = c;
            xi[cur[r]] = k;
            cur[r] += 1;
        }
    }
    let mut data = vec![0.0; n * m.ncols];
    let out = SendPtr(data.as_mut_ptr());
    (0..n).into_par_iter().for_each(|g| {
        let out = out;
        let inv = 1.0 / sd[g];
        unsafe {
            for c in 0..m.ncols {
                *out.0.add(g + c * n) = ((0.0 - mu[g]) * inv).min(10.0);
            }
            for pos in ptr[g]..ptr[g + 1] {
                *out.0.add(g + cols[pos] * n) = ((m.x[xi[pos]] - mu[g]) * inv).min(10.0);
            }
        }
    });
    let mut r = alloc_uninit(data.len());
    r.copy_from_slice(&data);
    r[r.len() / 3]
}

fn knownstats_colpar(m: &Csc, mu: &[f64], sd: &[f64]) -> f64 {
    let inv: Vec<f64> = sd.iter().map(|s| 1.0 / s).collect();
    let zero: Vec<f64> = (0..m.nrows).map(|g| ((0.0 - mu[g]) * inv[g]).min(10.0)).collect();
    let mut out = alloc_uninit(m.nrows * m.ncols);
    out.par_chunks_mut(m.nrows)
        .enumerate()
        .for_each(|(c, oc)| fill_col(m, c, oc, mu, &inv, &zero));
    out[out.len() / 3]
}

// ---------------- ComputeSNN ----------------
fn nn_index(n: usize, k: usize, seed: u64) -> Vec<f64> {
    let mut r = Lcg(seed);
    (0..n * k).map(|_| (1 + r.next() as usize % n) as f64).collect()
}

fn snn_build_inverted(data: &[f64], n: usize, k: usize) -> (Vec<Vec<(usize, u32)>>, Vec<usize>, Vec<usize>, Vec<u32>) {
    let mut neighbors = vec![Vec::with_capacity(k); n];
    for rank in 0..k {
        for i in 0..n {
            let m = data[rank * n + i] as usize - 1;
            let row: &mut Vec<(usize, u32)> = &mut neighbors[i];
            if let Some((_, c)) = row.iter_mut().find(|(idx, _)| *idx == m) {
                *c += 1;
            } else {
                row.push((m, 1));
            }
        }
    }
    let mut cnt = vec![0usize; n];
    for nb in &neighbors {
        for &(m, _) in nb {
            cnt[m] += 1;
        }
    }
    let mut ptr = vec![0usize; n + 1];
    for m in 0..n {
        ptr[m + 1] = ptr[m] + cnt[m];
    }
    let mut cells = vec![0usize; ptr[n]];
    let mut w = vec![0u32; ptr[n]];
    let mut nx = ptr.clone();
    for (cell, nb) in neighbors.iter().enumerate() {
        for &(m, wt) in nb {
            cells[nx[m]] = cell;
            w[nx[m]] = wt;
            nx[m] += 1;
        }
    }
    (neighbors, ptr, cells, w)
}

fn snn(data: &[f64], n: usize, k: usize, prune: f64, prune_first: bool) -> f64 {
    let (neighbors, ptr, cells, w) = snn_build_inverted(data, n, k);
    let kf = k as f64;
    let cols: Vec<(Vec<i32>, Vec<f64>)> = (0..n)
        .into_par_iter()
        .map_init(
            || (vec![0u32; n], vec![0u32; n], Vec::<usize>::new(), 1u32),
            |(acc, marks, touched, stamp), col| {
                *stamp += 1;
                for &(m, cw) in &neighbors[col] {
                    for pos in ptr[m]..ptr[m + 1] {
                        let row = cells[pos];
                        let wt = w[pos] * cw;
                        if marks[row] == *stamp {
                            acc[row] += wt;
                        } else {
                            marks[row] = *stamp;
                            acc[row] = wt;
                            touched.push(row);
                        }
                    }
                }
                let mut ic = Vec::new();
                let mut xc = Vec::new();
                if prune_first {
                    touched.retain(|&r| {
                        let v = acc[r] as f64;
                        v / (kf + (kf - v)) >= prune
                    });
                }
                touched.sort_unstable();
                for r in touched.drain(..) {
                    let v = acc[r] as f64;
                    let s = v / (kf + (kf - v));
                    if s >= prune {
                        ic.push(r as i32);
                        xc.push(s);
                    }
                }
                (ic, xc)
            },
        )
        .collect();
    let nnz: usize = cols.iter().map(|c| c.1.len()).sum();
    let mut x = Vec::with_capacity(nnz);
    let mut i = Vec::with_capacity(nnz);
    for (ic, xc) in cols {
        i.extend(ic);
        x.extend(xc);
    }
    // into_r_dgcmatrix: Doubles::new zero-fill + copy
    let mut rx = alloc_uninit(nnz);
    rx.iter_mut().for_each(|v| *v = 0.0);
    rx.copy_from_slice(&x);
    rx.len() as f64 + i.len() as f64
}

// ---------------- Standardize ----------------
fn standardize(values: &[f64], nrows: usize, ncols: usize, direct: bool) -> f64 {
    let body = |c: usize, oc: &mut [f64]| {
        let col = &values[c * nrows..(c + 1) * nrows];
        let mean = col.iter().sum::<f64>() / nrows as f64;
        let sd = (col.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / (nrows as f64 - 1.0)).sqrt();
        for r in 0..nrows {
            oc[r] = (col[r] - mean) / sd;
        }
    };
    if direct {
        let mut out = alloc_uninit(nrows * ncols);
        out.par_chunks_mut(nrows).enumerate().for_each(|(c, oc)| body(c, oc));
        out[7]
    } else {
        let mut tmp = vec![0.0; nrows * ncols];
        tmp.par_chunks_mut(nrows).enumerate().for_each(|(c, oc)| body(c, oc));
        let mut out = alloc_uninit(tmp.len());
        out.copy_from_slice(&tmp);
        out[7]
    }
}

// ---------------- CcaCrossprodMultiply ----------------
fn cca(left: &[f64], right: &[f64], x: &[f64], nf: usize, nl: usize, par: bool) -> f64 {
    let mut fs = vec![0.0; nf];
    for (c, &xv) in x.iter().enumerate() {
        for r in 0..nf {
            fs[r] += right[c * nf + r] * xv;
        }
    }
    let f = |c: usize| (0..nf).map(|r| left[c * nf + r] * fs[r]).sum::<f64>();
    let out: Vec<f64> = if par { (0..nl).into_par_iter().map(f).collect() } else { (0..nl).map(f).collect() };
    out[3]
}

fn main() {
    // Warm the rayon pool.
    (0..1000).into_par_iter().for_each(|_| {});
    println!("threads: {}", rayon::current_num_threads());

    let m = sparse(5000, 2000, 0.05, 1);
    println!("nnz {}", m.x.len());
    let a = lognorm_current(&m);
    let b = lognorm_fused_par(&m);
    assert!((a - b).abs() < 1e-12);
    println!("LogNorm  current {:.2} ms | fused+parallel+uninit {:.2} ms",
        median_ms(15, || lognorm_current(&m)), median_ms(15, || lognorm_fused_par(&m)));

    println!("FastSparseRowScale  current {:.2} ms | direct-into-output {:.2} ms",
        median_ms(9, || rowscale_current(&m)), median_ms(9, || rowscale_direct(&m)));

    let (mu, inv, _) = row_stats(&m);
    let sd: Vec<f64> = inv.iter().map(|v| 1.0 / v).collect();
    println!("KnownStats  current {:.2} ms | column-parallel direct {:.2} ms",
        median_ms(9, || knownstats_current(&m, &mu, &sd)), median_ms(9, || knownstats_colpar(&m, &mu, &sd)));

    let nn = nn_index(1500, 20, 5);
    let p = 1.0 / 15.0;
    assert_eq!(snn(&nn, 1500, 20, p, false), snn(&nn, 1500, 20, p, true));
    println!("  SNN serial setup (neighbor lists + inverted index) {:.3} ms",
        median_ms(21, || snn_build_inverted(&nn, 1500, 20).1[7] as f64));
    println!("ComputeSNN 1500/k20  current {:.3} ms | prune-before-sort {:.3} ms",
        median_ms(21, || snn(&nn, 1500, 20, p, false)), median_ms(21, || snn(&nn, 1500, 20, p, true)));

    let mut r = Lcg(9);
    let dense: Vec<f64> = (0..2000 * 120).map(|_| r.unif()).collect();
    println!("Standardize  current {:.3} ms | direct {:.3} ms",
        median_ms(31, || standardize(&dense, 2000, 120, false)), median_ms(31, || standardize(&dense, 2000, 120, true)));

    let left: Vec<f64> = (0..80 * 200).map(|_| r.unif()).collect();
    let right: Vec<f64> = (0..80 * 180).map(|_| r.unif()).collect();
    let xv: Vec<f64> = (0..180).map(|_| r.unif()).collect();
    println!("Cca  rayon {:.4} ms | serial {:.4} ms",
        median_ms(101, || cca(&left, &right, &xv, 80, 200, true)), median_ms(101, || cca(&left, &right, &xv, 80, 200, false)));
}
