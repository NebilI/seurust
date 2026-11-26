//! Prototypes for the four kernels still slower than Seurat 5.6 C++ after
//! roadmap items 1–7 (see `SLOW_KERNELS_PLAN.md`). Run with
//! `cargo run --release -- round2`.

use crate::{median_ms, Csc, Lcg};
use nalgebra::{DMatrix, SymmetricEigen, SymmetricTridiagonal};
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::Mutex;

fn rnorm(r: &mut Lcg, n: usize) -> Vec<f64> {
    (0..n)
        .map(|_| {
            let u1 = r.unif().max(1e-300);
            let u2 = r.unif();
            (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
        })
        .collect()
}

// ---------------- EigenGramPCA ----------------
fn gram_dgemm(x: &[f64], nf: usize, nc: usize) -> Vec<f64> {
    let mut g = vec![0.0; nf * nf];
    unsafe {
        matrixmultiply::dgemm(
            nf, nc, nf, 1.0,
            x.as_ptr(), 1, nf as isize,
            x.as_ptr(), nf as isize, 1,
            0.0, g.as_mut_ptr(), 1, nf as isize,
        );
    }
    g
}

/// Same dgemm, split over output column blocks on the rayon pool.
fn gram_dgemm_par(x: &[f64], nf: usize, nc: usize, blocks: usize) -> Vec<f64> {
    let mut g = vec![0.0; nf * nf];
    let bw = nf.div_ceil(blocks);
    g.par_chunks_mut(bw * nf).enumerate().for_each(|(b, out)| {
        let j0 = b * bw;
        let w = out.len() / nf;
        unsafe {
            matrixmultiply::dgemm(
                nf, nc, w, 1.0,
                x.as_ptr(), 1, nf as isize,
                x.as_ptr().add(j0), nf as isize, 1,
                0.0, out.as_mut_ptr(), 1, nf as isize,
            );
        }
    });
    g
}

/// Symmetric product: lower block triangle only, mirrored (≈ dsyrk flop count).
fn gram_syrk_blocked(x: &[f64], nf: usize, nc: usize, bs: usize) -> Vec<f64> {
    let mut g = vec![0.0; nf * nf];
    let nb = nf.div_ceil(bs);
    let pairs: Vec<(usize, usize)> = (0..nb).flat_map(|bj| (bj..nb).map(move |bi| (bi, bj))).collect();
    let gp = crate::SendPtr(g.as_mut_ptr());
    pairs.par_iter().for_each(|&(bi, bj)| {
        let gp = gp;
        let i0 = bi * bs;
        let j0 = bj * bs;
        let h = bs.min(nf - i0);
        let w = bs.min(nf - j0);
        unsafe {
            matrixmultiply::dgemm(
                h, nc, w, 1.0,
                x.as_ptr().add(i0), 1, nf as isize,
                x.as_ptr().add(j0), nf as isize, 1,
                0.0, gp.0.add(i0 + j0 * nf), 1, nf as isize,
            );
        }
    });
    for j in 0..nf {
        for i in 0..j {
            g[i + j * nf] = g[j + i * nf];
        }
    }
    g
}

fn eig_full_sorted(g: &[f64], nf: usize, npcs: usize) -> (Vec<f64>, Vec<f64>) {
    let e = SymmetricEigen::new(DMatrix::from_column_slice(nf, nf, g));
    let mut order: Vec<usize> = (0..nf).collect();
    order.sort_by(|&a, &b| e.eigenvalues[b].partial_cmp(&e.eigenvalues[a]).unwrap());
    let vals = order[..npcs].iter().map(|&j| e.eigenvalues[j]).collect();
    let mut vecs = vec![0.0; nf * npcs];
    for (k, &j) in order[..npcs].iter().enumerate() {
        for i in 0..nf {
            vecs[i + k * nf] = e.eigenvectors[(i, j)];
        }
    }
    (vals, vecs)
}

/// Householder tridiagonalization only: the O(n^3) part a LAPACK `dsyevr(RANGE='I')`
/// call still has to do. A lower bound for that approach in pure Rust.
fn tridiag_only(g: &[f64], nf: usize) -> f64 {
    let t = SymmetricTridiagonal::new(DMatrix::from_column_slice(nf, nf, g));
    let (d, _) = t.unpack_tridiagonal();
    d[nf / 2]
}

fn symv(g: &[f64], nf: usize, v: &[f64], out: &mut [f64]) {
    out.iter_mut().for_each(|o| *o = 0.0);
    for j in 0..nf {
        let vj = v[j];
        let col = &g[j * nf..(j + 1) * nf];
        for i in 0..nf {
            out[i] += col[i] * vj;
        }
    }
}

/// Thick-restart-free Lanczos with full reorthogonalization: `m` steps, then
/// the top `npcs` Ritz pairs. Returns (values, vectors, max relative residual).
fn lanczos(g: &[f64], nf: usize, npcs: usize, m: usize) -> (Vec<f64>, Vec<f64>, f64) {
    let m = m.min(nf);
    let mut v = vec![0.0; nf * m];
    let mut alpha = vec![0.0; m];
    let mut beta = vec![0.0; m];
    let inv = 1.0 / (nf as f64).sqrt();
    for i in 0..nf {
        v[i] = inv * (1.0 + 0.01 * (i % 7) as f64);
    }
    let n0: f64 = v[..nf].iter().map(|a| a * a).sum::<f64>().sqrt();
    v[..nf].iter_mut().for_each(|a| *a /= n0);
    let mut w = vec![0.0; nf];
    let mut steps = m;
    for j in 0..m {
        symv(g, nf, &v[j * nf..(j + 1) * nf], &mut w);
        alpha[j] = (0..nf).map(|i| w[i] * v[i + j * nf]).sum();
        for _pass in 0..2 {
            for k in 0..=j {
                let vk = &v[k * nf..(k + 1) * nf];
                let h: f64 = (0..nf).map(|i| w[i] * vk[i]).sum();
                for i in 0..nf {
                    w[i] -= h * vk[i];
                }
            }
        }
        if j + 1 == m {
            break;
        }
        let b: f64 = w.iter().map(|a| a * a).sum::<f64>().sqrt();
        if b < 1e-12 {
            steps = j + 1;
            break;
        }
        beta[j] = b;
        let (_, rest) = v.split_at_mut((j + 1) * nf);
        for i in 0..nf {
            rest[i] = w[i] / b;
        }
    }
    let mut t = DMatrix::<f64>::zeros(steps, steps);
    for j in 0..steps {
        t[(j, j)] = alpha[j];
        if j + 1 < steps {
            t[(j, j + 1)] = beta[j];
            t[(j + 1, j)] = beta[j];
        }
    }
    let e = SymmetricEigen::new(t);
    let mut order: Vec<usize> = (0..steps).collect();
    order.sort_by(|&a, &b| e.eigenvalues[b].partial_cmp(&e.eigenvalues[a]).unwrap());
    let mut vals = vec![0.0; npcs];
    let mut vecs = vec![0.0; nf * npcs];
    for (k, &j) in order[..npcs].iter().enumerate() {
        vals[k] = e.eigenvalues[j];
        for s in 0..steps {
            let y = e.eigenvectors[(s, j)];
            for i in 0..nf {
                vecs[i + k * nf] += v[i + s * nf] * y;
            }
        }
    }
    let mut worst: f64 = 0.0;
    let mut gv = vec![0.0; nf];
    for k in 0..npcs {
        let u = &vecs[k * nf..(k + 1) * nf];
        symv(g, nf, u, &mut gv);
        let r: f64 = (0..nf).map(|i| (gv[i] - vals[k] * u[i]).powi(2)).sum::<f64>().sqrt();
        worst = worst.max(r / vals[k].abs());
    }
    (vals, vecs, worst)
}

fn embeddings(x: &[f64], u: &[f64], nf: usize, nc: usize, npcs: usize) -> Vec<f64> {
    let mut e = vec![0.0; nc * npcs];
    unsafe {
        matrixmultiply::dgemm(
            nc, nf, npcs, 1.0,
            x.as_ptr(), nf as isize, 1,
            u.as_ptr(), 1, nf as isize,
            0.0, e.as_mut_ptr(), 1, nc as isize,
        );
    }
    e
}

fn bench_eigen_gram(nf: usize, nc: usize, npcs: usize, reps: usize) {
    let mut r = Lcg(25);
    let x = rnorm(&mut r, nf * nc);
    let g = gram_dgemm(&x, nf, nc);
    let gp = gram_dgemm_par(&x, nf, nc, 4);
    let gs = gram_syrk_blocked(&x, nf, nc, nf.div_ceil(4).max(32));
    let maxdiff = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(p, q)| (p - q).abs()).fold(0.0, f64::max);
    println!("EigenGramPCA {nf} x {nc}, {npcs} PCs");
    println!("  gram: dgemm serial {:.3} ms | dgemm 4 col-blocks {:.3} ms (max diff {:.1e}) | sym blocks {:.3} ms (max diff {:.1e})",
        median_ms(reps, || gram_dgemm(&x, nf, nc)[7]),
        median_ms(reps, || gram_dgemm_par(&x, nf, nc, 4)[7]), maxdiff(&g, &gp),
        median_ms(reps, || gram_syrk_blocked(&x, nf, nc, nf.div_ceil(4).max(32))[7]), maxdiff(&g, &gs));
    let (full_vals, _) = eig_full_sorted(&g, nf, npcs);
    println!("  eigen: nalgebra full SymmetricEigen + sort {:.3} ms | tridiagonalization only {:.3} ms",
        median_ms(reps, || eig_full_sorted(&g, nf, npcs).0[0]),
        median_ms(reps, || tridiag_only(&g, nf)));
    for m in [2 * npcs + 1, 40, 60] {
        if m > nf {
            continue;
        }
        let (vals, _, res) = lanczos(&g, nf, npcs, m);
        let rel = vals.iter().zip(&full_vals).map(|(a, b)| ((a - b) / b).abs()).fold(0.0, f64::max);
        println!("  eigen: Lanczos m={m:>3} {:.3} ms | eigenvalue rel err {:.1e}, max residual {:.1e}",
            median_ms(reps, || lanczos(&g, nf, npcs, m).0[0]), rel, res);
    }
    let (_, u) = eig_full_sorted(&g, nf, npcs);
    println!("  embeddings X^T U {:.3} ms",
        median_ms(reps, || embeddings(&x, &u, nf, nc, npcs)[3]));
    println!("  pipeline: current (dgemm + full eigen + emb) {:.3} ms",
        median_ms(reps, || {
            let g = gram_dgemm(&x, nf, nc);
            let (_, u) = eig_full_sorted(&g, nf, npcs);
            embeddings(&x, &u, nf, nc, npcs)[3]
        }));
}

// ---------------- CcaCrossprodMultiply ----------------
fn cca_scores(right: &[f64], x: &[f64], nf: usize) -> Vec<f64> {
    let mut fs = vec![0.0; nf];
    for (c, &xv) in x.iter().enumerate() {
        for r in 0..nf {
            fs[r] += right[c * nf + r] * xv;
        }
    }
    fs
}

fn cca_current(left: &[f64], right: &[f64], x: &[f64], nf: usize, nl: usize) -> Vec<f64> {
    let fs = cca_scores(right, x, nf);
    (0..nl)
        .map(|c| {
            let mut s = 0.0;
            for r in 0..nf {
                s += left[c * nf + r] * fs[r];
            }
            s
        })
        .collect()
}

/// Four output columns per pass: four independent add chains, each still summed
/// in row order, so every output is bit-identical to `cca_current`.
fn cca_blocked4(left: &[f64], right: &[f64], x: &[f64], nf: usize, nl: usize) -> Vec<f64> {
    let fs = cca_scores(right, x, nf);
    let mut out = vec![0.0; nl];
    let mut c = 0;
    while c + 4 <= nl {
        let (l0, l1, l2, l3) = (
            &left[c * nf..(c + 1) * nf],
            &left[(c + 1) * nf..(c + 2) * nf],
            &left[(c + 2) * nf..(c + 3) * nf],
            &left[(c + 3) * nf..(c + 4) * nf],
        );
        let (mut s0, mut s1, mut s2, mut s3) = (0.0, 0.0, 0.0, 0.0);
        for r in 0..nf {
            let f = fs[r];
            s0 += l0[r] * f;
            s1 += l1[r] * f;
            s2 += l2[r] * f;
            s3 += l3[r] * f;
        }
        out[c] = s0;
        out[c + 1] = s1;
        out[c + 2] = s2;
        out[c + 3] = s3;
        c += 4;
    }
    for c in c..nl {
        let mut s = 0.0;
        for r in 0..nf {
            s += left[c * nf + r] * fs[r];
        }
        out[c] = s;
    }
    out
}

fn bench_cca() {
    let mut r = Lcg(26);
    let (nf, nl, nr) = (80, 200, 180);
    let left = rnorm(&mut r, nf * nl);
    let right = rnorm(&mut r, nf * nr);
    let x = rnorm(&mut r, nr);
    let a = cca_current(&left, &right, &x, nf, nl);
    let b = cca_blocked4(&left, &right, &x, nf, nl);
    assert!(a.iter().zip(&b).all(|(p, q)| p.to_bits() == q.to_bits()));
    println!("CcaCrossprodMultiply 80 x 200 / 80 x 180 (compute only; bit-identical)");
    println!("  current serial {:.4} ms | 4-column blocked {:.4} ms | feature_scores pass alone {:.4} ms",
        median_ms(2001, || cca_current(&left, &right, &x, nf, nl)[3]),
        median_ms(2001, || cca_blocked4(&left, &right, &x, nf, nl)[3]),
        median_ms(2001, || cca_scores(&right, &x, nf)[3]));
}

// ---------------- SparseRowMeanVar ----------------
fn finish(sums: &[f64], sumsq: &[f64], n_cells: usize) -> (Vec<f64>, Vec<f64>) {
    let cols = n_cells as f64;
    let denom = (n_cells.max(2) - 1) as f64;
    let means = sums.iter().map(|s| s / cols).collect();
    let vars = sums.iter().zip(sumsq).map(|(s, q)| (q - s * s / cols) / denom).collect();
    (means, vars)
}

fn rmv_current(m: &Csc) -> (Vec<f64>, Vec<f64>) {
    let n = m.nrows;
    let (sums, sumsq) = (0..m.ncols)
        .into_par_iter()
        .fold(
            || (vec![0.0; n], vec![0.0; n]),
            |mut acc, col| {
                for k in m.p[col] as usize..m.p[col + 1] as usize {
                    let r = m.i[k] as usize;
                    let v = m.x[k];
                    acc.0[r] += v;
                    acc.1[r] += v * v;
                }
                acc
            },
        )
        .reduce(
            || (vec![0.0; n], vec![0.0; n]),
            |mut l, r| {
                for g in 0..n {
                    l.0[g] += r.0[g];
                    l.1[g] += r.1[g];
                }
                l
            },
        );
    finish(&sums, &sumsq, m.ncols)
}

/// Columns split into `t` contiguous nnz-balanced chunks, one interleaved
/// `[sum, sumsq]` accumulator per chunk, reduced by gene range in parallel.
fn rmv_chunked(m: &Csc, t: usize) -> (Vec<f64>, Vec<f64>) {
    let n = m.nrows;
    let nnz = m.x.len();
    let mut bounds = vec![0usize; t + 1];
    let mut c = 0;
    for k in 1..t {
        let target = (nnz * k / t) as i32;
        while c < m.ncols && m.p[c] < target {
            c += 1;
        }
        bounds[k] = c;
    }
    bounds[t] = m.ncols;
    let accs: Vec<Vec<[f64; 2]>> = (0..t)
        .into_par_iter()
        .map(|k| {
            let mut acc = vec![[0.0f64; 2]; n];
            let s = m.p[bounds[k]] as usize;
            let e = m.p[bounds[k + 1]] as usize;
            for (&r, &v) in m.i[s..e].iter().zip(&m.x[s..e]) {
                let a = &mut acc[r as usize];
                a[0] += v;
                a[1] += v * v;
            }
            acc
        })
        .collect();
    let cols = m.ncols as f64;
    let denom = (m.ncols.max(2) - 1) as f64;
    let mut means = vec![0.0; n];
    let mut vars = vec![0.0; n];
    let chunk = n.div_ceil(t);
    means
        .par_chunks_mut(chunk)
        .zip(vars.par_chunks_mut(chunk))
        .enumerate()
        .for_each(|(b, (mo, vo))| {
            let g0 = b * chunk;
            for (j, (mg, vg)) in mo.iter_mut().zip(vo.iter_mut()).enumerate() {
                let (mut s, mut q) = (0.0, 0.0);
                for acc in &accs {
                    s += acc[g0 + j][0];
                    q += acc[g0 + j][1];
                }
                *mg = s / cols;
                *vg = (q - s * s / cols) / denom;
            }
        });
    (means, vars)
}

fn rmv_serial(m: &Csc) -> (Vec<f64>, Vec<f64>) {
    let mut acc = vec![[0.0f64; 2]; m.nrows];
    for (&r, &v) in m.i.iter().zip(&m.x) {
        let a = &mut acc[r as usize];
        a[0] += v;
        a[1] += v * v;
    }
    let sums: Vec<f64> = acc.iter().map(|a| a[0]).collect();
    let sq: Vec<f64> = acc.iter().map(|a| a[1]).collect();
    finish(&sums, &sq, m.ncols)
}

fn bench_row_mean_var(m: &Csc) {
    let a = rmv_current(m);
    for t in [3, 4] {
        let b = rmv_chunked(m, t);
        assert!(a.0.iter().zip(&b.0).all(|(p, q)| p.to_bits() == q.to_bits()));
        assert!(a.1.iter().zip(&b.1).all(|(p, q)| p.to_bits() == q.to_bits()));
    }
    println!("SparseRowMeanVar 5000 x 2000, 5% (integer counts, so sums are order-independent)");
    let pool3 = rayon::ThreadPoolBuilder::new().num_threads(3).build().unwrap();
    println!("  current fold/reduce {:.3} ms | chunked t=4 {:.3} ms | chunked t=3 (3-thread pool) {:.3} ms | serial {:.3} ms",
        median_ms(201, || rmv_current(m).1[7]),
        median_ms(201, || rmv_chunked(m, 4).1[7]),
        pool3.install(|| median_ms(201, || rmv_chunked(m, 3).1[7])),
        median_ms(201, || rmv_serial(m).1[7]));
}

// ---------------- fast_dist ----------------
fn fast_dist_compute(x: &[f64], y: &[f64], nb: &[Vec<usize>], n: usize, d: usize, out: &mut [Vec<f64>]) {
    for (cell, (ids, o)) in nb.iter().zip(out.iter_mut()).enumerate() {
        for (j, &id) in ids.iter().enumerate() {
            let mut s = 0.0;
            for c in 0..d {
                let t = x[cell + c * n] - y[id + c * n];
                s += t * t;
            }
            o[j] = s.sqrt();
        }
    }
}

/// Transpose `x` and `y` to row-major once so each distance reads two
/// contiguous `d`-length rows. Same per-element summation order.
fn fast_dist_transposed(x: &[f64], y: &[f64], nb: &[Vec<usize>], n: usize, d: usize, out: &mut [Vec<f64>]) {
    let tr = |m: &[f64]| {
        let mut t = vec![0.0; n * d];
        for c in 0..d {
            for r in 0..n {
                t[r * d + c] = m[r + c * n];
            }
        }
        t
    };
    let xt = tr(x);
    let yt = tr(y);
    for (cell, (ids, o)) in nb.iter().zip(out.iter_mut()).enumerate() {
        let xs = &xt[cell * d..(cell + 1) * d];
        for (j, &id) in ids.iter().enumerate() {
            let ys = &yt[id * d..(id + 1) * d];
            let mut s = 0.0;
            for c in 0..d {
                let t = xs[c] - ys[c];
                s += t * t;
            }
            o[j] = s.sqrt();
        }
    }
}

/// Four neighbors per pass over the dimensions: four independent add chains,
/// each in the same column order as `fast_dist_compute`.
fn fast_dist_interleaved(x: &[f64], y: &[f64], nb: &[Vec<usize>], n: usize, d: usize, out: &mut [Vec<f64>]) {
    for (cell, (ids, o)) in nb.iter().zip(out.iter_mut()).enumerate() {
        let mut j = 0;
        while j + 4 <= ids.len() {
            let (a, b, c4, e) = (ids[j], ids[j + 1], ids[j + 2], ids[j + 3]);
            let (mut s0, mut s1, mut s2, mut s3) = (0.0, 0.0, 0.0, 0.0);
            for c in 0..d {
                let xv = x[cell + c * n];
                let base = c * n;
                let t0 = xv - y[a + base];
                let t1 = xv - y[b + base];
                let t2 = xv - y[c4 + base];
                let t3 = xv - y[e + base];
                s0 += t0 * t0;
                s1 += t1 * t1;
                s2 += t2 * t2;
                s3 += t3 * t3;
            }
            o[j] = s0.sqrt();
            o[j + 1] = s1.sqrt();
            o[j + 2] = s2.sqrt();
            o[j + 3] = s3.sqrt();
            j += 4;
        }
        for j in j..ids.len() {
            let mut s = 0.0;
            for c in 0..d {
                let t = x[cell + c * n] - y[ids[j] + c * n];
                s += t * t;
            }
            o[j] = s.sqrt();
        }
    }
}

/// Stand-in for extendr 0.9 `ownership::protect`/`unprotect`: a global
/// `Mutex<HashMap<SEXP, Object>>` (SipHash) touched once per R object.
struct Ownership {
    objects: HashMap<usize, (usize, usize)>,
    preservation: Vec<usize>,
}

fn bench_fast_dist() {
    let (n, d, k) = (800, 20, 10);
    let mut r = Lcg(19);
    let x = rnorm(&mut r, n * d);
    let y = rnorm(&mut r, n * d);
    let nb: Vec<Vec<usize>> = (0..n).map(|_| (0..k).map(|_| r.next() as usize % n).collect()).collect();
    let mut out = vec![vec![0.0; k]; n];
    println!("fast_dist 800 x 20, 10 neighbors");
    let mut out2 = vec![vec![0.0; k]; n];
    fast_dist_compute(&x, &y, &nb, n, d, &mut out);
    fast_dist_transposed(&x, &y, &nb, n, d, &mut out2);
    assert!(out.iter().flatten().zip(out2.iter().flatten()).all(|(a, b)| a.to_bits() == b.to_bits()));
    let mut out3 = vec![vec![0.0; k]; n];
    fast_dist_interleaved(&x, &y, &nb, n, d, &mut out3);
    assert!(out.iter().flatten().zip(out3.iter().flatten()).all(|(a, b)| a.to_bits() == b.to_bits()));
    println!("  distance arithmetic: current strided {:.4} ms | transpose once {:.4} ms | 4 neighbors per pass {:.4} ms (all bit-identical)",
        median_ms(2001, || { fast_dist_compute(&x, &y, &nb, n, d, &mut out); out[3][2] }),
        median_ms(2001, || { fast_dist_transposed(&x, &y, &nb, n, d, &mut out2); out2[3][2] }),
        median_ms(2001, || { fast_dist_interleaved(&x, &y, &nb, n, d, &mut out3); out3[3][2] }));
    let own = Mutex::new(Ownership { objects: HashMap::with_capacity(100_000), preservation: vec![0; 100_000] });
    let mut next = 1usize << 20;
    // Per call: 800 list elements read as Robj + 800 Doubles allocated + the list,
    // each protected once and unprotected once.
    println!("  modelled extendr protect/unprotect for 1,601 objects {:.4} ms",
        median_ms(2001, || {
            let mut acc = 0usize;
            for _ in 0..1601 {
                next = next.wrapping_add(4096) % (1usize << 40);
                let sexp = next;
                {
                    let mut o = own.lock().unwrap();
                    let idx = o.objects.len() % 100_000;
                    o.preservation[idx] = sexp;
                    o.objects.entry(sexp).or_insert((1, idx)).0 += 1;
                }
                {
                    let mut o = own.lock().unwrap();
                    if let Some(e) = o.objects.get_mut(&sexp) {
                        e.0 -= 1;
                        acc += e.1;
                    }
                    if o.objects.len() > 90_000 {
                        o.objects.clear();
                    }
                }
            }
            acc as f64
        }));
}

pub fn run(m: &Csc) {
    bench_eigen_gram(120, 400, 10, 101);
    bench_eigen_gram(1000, 3000, 30, 3);
    bench_cca();
    bench_row_mean_var(m);
    bench_fast_dist();
}
