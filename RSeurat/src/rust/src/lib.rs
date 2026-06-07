mod data_manipulation;
mod fast_nn_dist;
mod integration;
mod modularity;
mod snn;
mod sparse;
mod stats;
mod utils;

use extendr_api::prelude::*;

use data_manipulation::{
    fast_cov_impl, fast_cov_mats_impl, fast_exp_mean_impl, fast_log_vmr_impl,
    fast_rbind_impl, fast_sparse_row_scale_impl, fast_sparse_row_scale_with_known_stats_impl,
    graph_to_neighbor_helper_impl, log_norm_impl, replace_cols_impl, row_merge_matrices_impl,
    row_var_impl, run_umi_sampling_impl, run_umi_sampling_per_cell_impl, sparse_row_var2_impl,
    sparse_row_var_impl, sparse_row_var_std_impl, standardize_impl,
};
use fast_nn_dist::fast_dist_impl;
use integration::{find_weights_impl, integrate_data_impl, score_helper_impl};
use modularity::run_modularity_clustering_impl;
use snn::{
    compute_snn_impl, compute_snn_to_r_impl, direct_snn_to_file_impl,
    snn_smallest_nonzero_dist_impl, write_edge_file_impl,
};
use sparse::{strings_to_str_vec, vec_from_doubles, vec_from_integers, CscSlots, CscView, CsrSlots};
use stats::{
    row_mean_dgcmatrix_impl, row_sum_dgcmatrix_impl, row_var_dgcmatrix_impl,
};

#[extendr]
fn row_sum_dgcmatrix(x: Doubles, i: Integers, rows: i32, _cols: i32) -> Doubles {
    row_sum_dgcmatrix_impl(&x, &i, rows)
}

#[extendr]
fn row_mean_dgcmatrix(x: Doubles, i: Integers, rows: i32, cols: i32) -> Doubles {
    row_mean_dgcmatrix_impl(&x, &i, rows, cols)
}

#[extendr]
fn row_var_dgcmatrix(x: Doubles, i: Integers, rows: i32, cols: i32) -> Doubles {
    row_var_dgcmatrix_impl(&x, &i, rows, cols)
}

#[extendr]
fn log_norm(
    x: Doubles,
    i: Integers,
    p: Integers,
    nrows: i32,
    ncols: i32,
    scale_factor: i32,
    display_progress: bool,
) -> Doubles {
    let col_sums = {
        let view = CscView::from_slots(&x, &i, &p, nrows, ncols);
        view.col_sums()
    };
    let p_slice = p.as_robj().as_integer_slice().expect("integer p");
    let x_in = x.as_robj().as_real_slice().expect("numeric x");
    let mut out = Doubles::new(x_in.len());
    let x_out = out
        .as_robj_mut()
        .as_real_slice_mut()
        .expect("numeric output");
    x_out.copy_from_slice(x_in);
    log_norm_impl(
        x_out,
        p_slice,
        &col_sums,
        ncols as usize,
        scale_factor,
        display_progress,
    );
    out
}

#[extendr]
fn standardize(mat: RMatrix<f64>, display_progress: bool) -> RMatrix<f64> {
    standardize_impl(&mat, display_progress)
}

#[extendr]
fn fast_cov(mat: RMatrix<f64>, center: bool) -> RMatrix<f64> {
    fast_cov_impl(&mat, center)
}

#[extendr]
fn fast_cov_mats(mat1: RMatrix<f64>, mat2: RMatrix<f64>, center: bool) -> RMatrix<f64> {
    fast_cov_mats_impl(&mat1, &mat2, center)
}

#[extendr]
fn fast_rbind(mat1: RMatrix<f64>, mat2: RMatrix<f64>) -> RMatrix<f64> {
    fast_rbind_impl(&mat1, &mat2)
}

#[extendr]
fn row_var(mat: RMatrix<f64>) -> Doubles {
    row_var_impl(&mat)
}

#[extendr]
fn fast_exp_mean(
    x: Doubles,
    i: Integers,
    p: Integers,
    nrows: i32,
    ncols: i32,
    display_progress: bool,
) -> Doubles {
    let view = CscView::from_slots(&x, &i, &p, nrows, ncols);
    fast_exp_mean_impl(view, display_progress)
}

#[extendr]
fn sparse_row_var(
    x: Doubles,
    i: Integers,
    p: Integers,
    nrows: i32,
    ncols: i32,
    display_progress: bool,
) -> Doubles {
    let view = CscView::from_slots(&x, &i, &p, nrows, ncols);
    sparse_row_var_impl(view, display_progress)
}

#[extendr]
fn sparse_row_var2(
    x: Doubles,
    i: Integers,
    p: Integers,
    nrows: i32,
    ncols: i32,
    mu: Doubles,
    display_progress: bool,
) -> Doubles {
    let view = CscView::from_slots(&x, &i, &p, nrows, ncols);
    let mu_vec = vec_from_doubles(&mu);
    sparse_row_var2_impl(view, &mu_vec, display_progress)
}

#[extendr]
fn sparse_row_var_std(
    x: Doubles,
    i: Integers,
    p: Integers,
    nrows: i32,
    ncols: i32,
    mu: Doubles,
    sd: Doubles,
    vmax: f64,
    display_progress: bool,
) -> Doubles {
    let view = CscView::from_slots(&x, &i, &p, nrows, ncols);
    let mu_vec = vec_from_doubles(&mu);
    let sd_vec = vec_from_doubles(&sd);
    sparse_row_var_std_impl(view, &mu_vec, &sd_vec, vmax, display_progress)
}

#[extendr]
fn fast_log_vmr(
    x: Doubles,
    i: Integers,
    p: Integers,
    nrows: i32,
    ncols: i32,
    display_progress: bool,
) -> Doubles {
    let view = CscView::from_slots(&x, &i, &p, nrows, ncols);
    fast_log_vmr_impl(view, display_progress)
}

#[extendr]
fn fast_sparse_row_scale(
    x: Doubles,
    i: Integers,
    p: Integers,
    nrows: i32,
    ncols: i32,
    scale: bool,
    center: bool,
    scale_max: f64,
    display_progress: bool,
) -> RMatrix<f64> {
    let view = CscView::from_slots(&x, &i, &p, nrows, ncols);
    fast_sparse_row_scale_impl(view, scale, center, scale_max, display_progress)
}

#[extendr]
fn fast_sparse_row_scale_with_known_stats(
    x: Doubles,
    i: Integers,
    p: Integers,
    nrows: i32,
    ncols: i32,
    mu: Doubles,
    sigma: Doubles,
    scale: bool,
    center: bool,
    scale_max: f64,
    display_progress: bool,
) -> RMatrix<f64> {
    let view = CscView::from_slots(&x, &i, &p, nrows, ncols);
    let mu_vec = vec_from_doubles(&mu);
    let sigma_vec = vec_from_doubles(&sigma);
    fast_sparse_row_scale_with_known_stats_impl(
        view,
        &mu_vec,
        &sigma_vec,
        scale,
        center,
        scale_max,
        display_progress,
    )
}

#[extendr(use_rng = true)]
fn run_umi_sampling(
    x: Doubles,
    i: Integers,
    p: Integers,
    nrows: i32,
    ncols: i32,
    sample_val: i32,
    upsample: bool,
    _display_progress: bool,
) -> List {
    let mat = CscSlots::from_r(x, i, p, nrows, ncols);
    run_umi_sampling_impl(mat, sample_val, upsample).to_r_list()
}

#[extendr(use_rng = true)]
fn run_umi_sampling_per_cell(
    x: Doubles,
    i: Integers,
    p: Integers,
    nrows: i32,
    ncols: i32,
    sample_val: Doubles,
    upsample: bool,
    _display_progress: bool,
) -> List {
    let mat = CscSlots::from_r(x, i, p, nrows, ncols);
    let sample_vec = vec_from_doubles(&sample_val);
    run_umi_sampling_per_cell_impl(mat, &sample_vec, upsample).to_r_list()
}

#[extendr]
fn row_merge_matrices(
    x1: Doubles,
    j1: Integers,
    p1: Integers,
    nrows1: i32,
    ncols1: i32,
    x2: Doubles,
    j2: Integers,
    p2: Integers,
    nrows2: i32,
    ncols2: i32,
    mat1_rownames: Strings,
    mat2_rownames: Strings,
    all_rownames: Strings,
) -> List {
    let mat1 = CsrSlots::from_r(x1, j1, p1, nrows1, ncols1);
    let mat2 = CsrSlots::from_r(x2, j2, p2, nrows2, ncols2);
    let names1 = strings_to_str_vec(mat1_rownames);
    let names2 = strings_to_str_vec(mat2_rownames);
    let all_names = strings_to_str_vec(all_rownames);
    row_merge_matrices_impl(mat1, mat2, &names1, &names2, &all_names).to_r_list()
}

#[extendr]
fn replace_cols(
    x: Doubles,
    i: Integers,
    p: Integers,
    nrows: i32,
    ncols: i32,
    col_idx: Doubles,
    rx: Doubles,
    ri: Integers,
    rp: Integers,
    rnrows: i32,
    rncols: i32,
) -> List {
    let mat = CscSlots::from_r(x, i, p, nrows, ncols);
    let replacement = CscSlots::from_r(rx, ri, rp, rnrows, rncols);
    let cols: Vec<i32> = vec_from_doubles(&col_idx)
        .into_iter()
        .map(|v| v as i32)
        .collect();
    replace_cols_impl(mat, &cols, replacement).to_r_list()
}

#[extendr]
fn graph_to_neighbor_helper(
    x: Doubles,
    i: Integers,
    p: Integers,
    nrows: i32,
    ncols: i32,
) -> Robj {
    let mat = CscSlots::from_r(x, i, p, nrows, ncols);
    graph_to_neighbor_helper_impl(mat)
}

#[extendr]
fn fast_dist(x: RMatrix<f64>, y: RMatrix<f64>, n: List) -> Robj {
    fast_dist_impl(&x, &y, &n)
}

#[extendr]
fn find_weights(
    cells2: Doubles,
    distances: RMatrix<f64>,
    anchor_cells2: Strings,
    integration_matrix_rownames: Strings,
    cell_index: RMatrix<f64>,
    anchor_score: Doubles,
    min_dist: f64,
    sd: f64,
    display_progress: bool,
) -> List {
    let cells: Vec<i32> = vec_from_doubles(&cells2)
        .into_iter()
        .map(|v| v as i32)
        .collect();
    find_weights_impl(
        &cells,
        &distances,
        &strings_to_str_vec(anchor_cells2),
        &strings_to_str_vec(integration_matrix_rownames),
        &cell_index,
        &anchor_score.iter().map(|v| v.0).collect::<Vec<_>>(),
        min_dist,
        sd,
        display_progress,
    )
    .to_r_list()
}

#[extendr]
fn integrate_data(
    ix: Doubles,
    ii: Integers,
    ip: Integers,
    inrows: i32,
    incols: i32,
    wx: Doubles,
    wi: Integers,
    wp: Integers,
    wnrows: i32,
    wncols: i32,
    ex: Doubles,
    ei: Integers,
    ep: Integers,
    enrows: i32,
    encols: i32,
) -> List {
    let integration_matrix = CscSlots::from_r(ix, ii, ip, inrows, incols);
    let weights = CscSlots::from_r(wx, wi, wp, wnrows, wncols);
    let expression = CscSlots::from_r(ex, ei, ep, enrows, encols);
    integrate_data_impl(integration_matrix, weights, expression).to_r_list()
}

#[extendr]
fn score_helper(
    x: Doubles,
    i: Integers,
    p: Integers,
    nrows: i32,
    ncols: i32,
    query_pca: RMatrix<f64>,
    query_dists: RMatrix<f64>,
    corrected_nns: RMatrix<f64>,
    k_snn: i32,
    subtract_first_nn: bool,
    display_progress: bool,
) -> Doubles {
    let snn = CscSlots::from_r(x, i, p, nrows, ncols);
    score_helper_impl(
        snn,
        &query_pca,
        &query_dists,
        &corrected_nns,
        k_snn,
        subtract_first_nn,
        display_progress,
    )
}

#[extendr]
fn compute_snn(nn_ranked: RMatrix<f64>, prune: f64) -> extendr_api::Result<Robj> {
    compute_snn_to_r_impl(&nn_ranked, prune)
}

#[extendr]
fn write_edge_file(
    x: Doubles,
    i: Integers,
    p: Integers,
    nrows: i32,
    ncols: i32,
    filename: &str,
    display_progress: bool,
) {
    let snn = CscSlots::from_r(x, i, p, nrows, ncols);
    write_edge_file_impl(&snn, filename, display_progress);
}

#[extendr]
fn direct_snn_to_file(
    nn_ranked: RMatrix<f64>,
    prune: f64,
    display_progress: bool,
    filename: &str,
) -> List {
    direct_snn_to_file_impl(&nn_ranked, prune, display_progress, filename).into_r_list()
}

#[extendr]
fn run_modularity_clustering(
    x: Doubles,
    i: Integers,
    p: Integers,
    nrows: i32,
    ncols: i32,
    modularity_function: i32,
    resolution: f64,
    algorithm: i32,
    n_random_starts: i32,
    n_iterations: i32,
    random_seed: i32,
    print_output: bool,
    edgefilename: &str,
) -> Result<Integers, extendr_api::Error> {
    let clusters = run_modularity_clustering_impl(
        x.as_robj().as_real_slice().expect("numeric x"),
        i.as_robj().as_integer_slice().expect("integer i"),
        p.as_robj().as_integer_slice().expect("integer p"),
        nrows,
        ncols,
        modularity_function,
        resolution,
        algorithm,
        n_random_starts,
        n_iterations,
        random_seed,
        print_output,
        edgefilename,
    )
    .map_err(|msg| extendr_api::Error::Other(msg.into()))?;

    let mut out = Integers::new(clusters.len());
    out.as_robj_mut()
        .as_integer_slice_mut()
        .expect("cluster labels")
        .copy_from_slice(&clusters);
    Ok(out)
}

#[extendr]
fn snn_smallest_nonzero_dist(
    x: Doubles,
    i: Integers,
    p: Integers,
    nrows: i32,
    ncols: i32,
    mat: RMatrix<f64>,
    n: i32,
    nearest_dist: Doubles,
) -> Doubles {
    let snn = CscSlots::from_r(x, i, p, nrows, ncols);
    let nd = vec_from_doubles(&nearest_dist);
    snn_smallest_nonzero_dist_impl(snn, &mat, n, &nd)
}

extendr_module! {
    mod RSeurat;
    fn row_sum_dgcmatrix;
    fn row_mean_dgcmatrix;
    fn row_var_dgcmatrix;
    fn log_norm;
    fn standardize;
    fn fast_cov;
    fn fast_cov_mats;
    fn fast_rbind;
    fn row_var;
    fn fast_exp_mean;
    fn sparse_row_var;
    fn sparse_row_var2;
    fn sparse_row_var_std;
    fn fast_log_vmr;
    fn fast_sparse_row_scale;
    fn fast_sparse_row_scale_with_known_stats;
    fn run_umi_sampling;
    fn run_umi_sampling_per_cell;
    fn row_merge_matrices;
    fn replace_cols;
    fn graph_to_neighbor_helper;
    fn fast_dist;
    fn find_weights;
    fn integrate_data;
    fn score_helper;
    fn compute_snn;
    fn write_edge_file;
    fn direct_snn_to_file;
    fn snn_smallest_nonzero_dist;
    fn run_modularity_clustering;
}
