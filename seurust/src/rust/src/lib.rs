#![allow(dead_code, non_snake_case)]

mod data_manipulation;
mod fast_nn_dist;
mod integration;
mod markers;
mod modularity;
mod pca;
mod sct;
mod snn;
mod sparse;
mod stats;
mod utils;

use extendr_api::prelude::*;

use data_manipulation::{
    fast_cov_impl, fast_cov_mats_impl, fast_dense_row_scale_impl, fast_exp_mean_impl,
    fast_log_vmr_impl, fast_rbind_impl, fast_sparse_row_scale_impl,
    fast_sparse_row_scale_with_known_stats_impl, graph_to_neighbor_helper_impl, log_norm_impl,
    replace_cols_impl, row_merge_matrices_impl, row_var_impl, run_umi_sampling_impl,
    run_umi_sampling_per_cell_impl, sparse_row_mean_var_impl, sparse_row_mean_var_legacy_impl,
    sparse_row_var2_impl, sparse_row_var_impl, sparse_row_var_std_impl,
    sparse_row_var_std_legacy_rows_impl, standardize_impl,
};
use fast_nn_dist::fast_dist_impl;
use integration::{
    cca_crossprod_multiply_impl, count_anchor_shared_neighbors_impl, find_weights_impl,
    integrate_data_impl, score_helper_impl,
};
use markers::find_all_markers_sparse_fold_change_stats_impl;
use modularity::{run_modularity_clustering_impl, run_modularity_clustering_multi_impl};
use pca::eigen_gram_pca_impl;
use sct::{sct_pearson_residual_matrix_impl, sct_residual_stats_and_corrected_impl};
use snn::{
    compute_snn_to_r_impl, direct_snn_to_file_impl, snn_smallest_nonzero_dist_impl,
    write_edge_file_impl,
};
use sparse::{
    doubles_from_slice, strings_to_str_vec, vec_from_doubles, vec_from_integers, CscSlots, CscView,
    CsrView,
};
use stats::{row_mean_dgcmatrix_impl, row_sum_dgcmatrix_impl, row_var_dgcmatrix_impl};

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
    let mut out = doubles_from_slice(x_in);
    let x_out = out
        .as_robj_mut()
        .as_real_slice_mut()
        .expect("numeric output");
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
fn standardize(mat: RMatrix<f64>, display_progress: bool, _nthreads: i32) -> RMatrix<f64> {
    standardize_impl(&mat, display_progress)
}

#[extendr]
fn fast_cov(mat: RMatrix<f64>, center: bool) -> RMatrix<f64> {
    fast_cov_impl(&mat, center)
}

#[extendr]
fn fast_cov_mats(
    mat1: RMatrix<f64>,
    mat2: RMatrix<f64>,
    center: bool,
) -> extendr_api::Result<RMatrix<f64>> {
    fast_cov_mats_impl(&mat1, &mat2, center)
}

#[extendr]
fn fast_rbind(mat1: RMatrix<f64>, mat2: RMatrix<f64>) -> extendr_api::Result<RMatrix<f64>> {
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
    let mu_slice = mu.as_robj().as_real_slice().expect("numeric mu");
    sparse_row_var2_impl(view, mu_slice, display_progress)
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
    features: Integers,
    scale: bool,
    center: bool,
    scale_max: f64,
    display_progress: bool,
) -> RMatrix<f64> {
    let view = CscView::from_slots(&x, &i, &p, nrows, ncols);
    let features = vec_from_integers(&features);
    fast_sparse_row_scale_impl(view, &features, scale, center, scale_max, display_progress)
}

#[extendr]
fn fast_dense_row_scale(
    mat: RMatrix<f64>,
    features: Integers,
    scale: bool,
    center: bool,
    scale_max: f64,
    _nthreads: i32,
    display_progress: bool,
) -> RMatrix<f64> {
    let features = vec_from_integers(&features);
    fast_dense_row_scale_impl(&mat, &features, scale, center, scale_max, display_progress)
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
) -> extendr_api::Result<Robj> {
    let mat1 = CsrView::from_slots(&x1, &j1, &p1, nrows1, ncols1);
    let mat2 = CsrView::from_slots(&x2, &j2, &p2, nrows2, ncols2);
    row_merge_matrices_impl(mat1, mat2, &mat1_rownames, &mat2_rownames, &all_rownames)
        .into_r_dgcmatrix()
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
) -> extendr_api::Result<Robj> {
    let mat = CscView::from_slots(&x, &i, &p, nrows, ncols);
    let replacement = CscView::from_slots(&rx, &ri, &rp, rnrows, rncols);
    let cols: Vec<i32> = vec_from_doubles(&col_idx)
        .into_iter()
        .map(|v| v as i32)
        .collect();
    replace_cols_impl(mat, &cols, replacement)
}

#[extendr]
fn graph_to_neighbor_helper(
    x: Doubles,
    i: Integers,
    p: Integers,
    nrows: i32,
    ncols: i32,
) -> extendr_api::Result<Robj> {
    let mat = CscSlots::from_r(x, i, p, nrows, ncols);
    graph_to_neighbor_helper_impl(mat)
}

#[extendr]
fn fast_dist(x: RMatrix<f64>, y: RMatrix<f64>, n: List) -> extendr_api::Result<Robj> {
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
    _nthreads: i32,
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
) -> extendr_api::Result<Robj> {
    let integration_matrix = CscView::from_slots(&ix, &ii, &ip, inrows, incols);
    let weights = CscView::from_slots(&wx, &wi, &wp, wnrows, wncols);
    let expression = CscView::from_slots(&ex, &ei, &ep, enrows, encols);
    integrate_data_impl(integration_matrix, weights, expression).into_r_dgcmatrix()
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
    _nthreads: i32,
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
fn compute_snn(nn_ranked: RMatrix<f64>, prune: f64, _nthreads: i32) -> extendr_api::Result<Robj> {
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
) -> extendr_api::Result<()> {
    let snn = CscView::from_slots(&x, &i, &p, nrows, ncols);
    write_edge_file_impl(&snn, filename, display_progress)
}

#[extendr]
fn direct_snn_to_file(
    nn_ranked: RMatrix<f64>,
    prune: f64,
    display_progress: bool,
    filename: &str,
) -> extendr_api::Result<Robj> {
    direct_snn_to_file_impl(&nn_ranked, prune, display_progress, filename)?.into_r_dgcmatrix()
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
    n_threads: i32,
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
        n_threads,
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

#[extendr]
fn sparse_row_mean_var(
    x: Doubles,
    i: Integers,
    p: Integers,
    nrows: i32,
    ncols: i32,
    _nthreads: i32,
    _display_progress: bool,
) -> List {
    let view = CscView::from_slots(&x, &i, &p, nrows, ncols);
    let (mean, variance) = sparse_row_mean_var_impl(view);
    list!(
        mean = Doubles::from_values(mean),
        variance = Doubles::from_values(variance)
    )
}

#[extendr]
fn sparse_row_mean_var_legacy(
    x: Doubles,
    i: Integers,
    p: Integers,
    nrows: i32,
    ncols: i32,
    _display_progress: bool,
) -> List {
    let view = CscView::from_slots(&x, &i, &p, nrows, ncols);
    let (mean, variance) = sparse_row_mean_var_legacy_impl(view);
    list!(
        mean = Doubles::from_values(mean),
        variance = Doubles::from_values(variance)
    )
}

#[extendr]
fn sparse_row_var_std_legacy_rows(
    x: Doubles,
    i: Integers,
    p: Integers,
    mu: Doubles,
    sd: Doubles,
    vmax: f64,
    _rows: i32,
    _cols: i32,
    rows_use: Integers,
    _nthreads: i32,
    _display_progress: bool,
) -> Doubles {
    let view = CscView::from_slots(&x, &i, &p, _rows, _cols);
    let mu = vec_from_doubles(&mu);
    let sd = vec_from_doubles(&sd);
    let rows_use = vec_from_integers(&rows_use);
    Doubles::from_values(sparse_row_var_std_legacy_rows_impl(
        view, &mu, &sd, vmax, &rows_use,
    ))
}

#[extendr]
fn find_all_markers_sparse_fold_change_stats(
    x: Doubles,
    i: Integers,
    p: Integers,
    rows: i32,
    cols: i32,
    groups: Integers,
    n_groups: i32,
    log_normalize: bool,
    _nthreads: i32,
) -> extendr_api::Result<Robj> {
    let view = CscView::from_slots(&x, &i, &p, rows, cols);
    let groups = vec_from_integers(&groups);
    find_all_markers_sparse_fold_change_stats_impl(view, &groups, n_groups, log_normalize)
}

#[extendr]
fn cca_crossprod_multiply(
    left: RMatrix<f64>,
    right: RMatrix<f64>,
    x: Doubles,
) -> extendr_api::Result<Doubles> {
    let x = vec_from_doubles(&x);
    cca_crossprod_multiply_impl(&left, &right, &x)
}

#[extendr]
fn count_anchor_shared_neighbors(
    indices_aa: RMatrix<f64>,
    indices_ab: RMatrix<f64>,
    indices_ba: RMatrix<f64>,
    indices_bb: RMatrix<f64>,
    anchor_cell1: Integers,
    anchor_cell2: Integers,
    offset: i32,
    k_score: i32,
    _nthreads: i32,
) -> Integers {
    let cell1 = vec_from_integers(&anchor_cell1);
    let cell2 = vec_from_integers(&anchor_cell2);
    count_anchor_shared_neighbors_impl(
        &indices_aa,
        &indices_ab,
        &indices_ba,
        &indices_bb,
        &cell1,
        &cell2,
        offset,
        k_score,
    )
}

#[extendr]
fn eigen_gram_pca(
    object: RMatrix<f64>,
    npcs: i32,
    weight_by_var: bool,
    nthreads: i32,
) -> extendr_api::Result<Robj> {
    eigen_gram_pca_impl(&object, npcs, weight_by_var, nthreads)
}

#[extendr]
fn sct_residual_stats_and_corrected(
    x: Doubles,
    i: Integers,
    p: Integers,
    rows: i32,
    cols: i32,
    theta: Doubles,
    intercept: Doubles,
    slope: Doubles,
    log_umi: Doubles,
    target_log_umi: f64,
    min_var: f64,
    residual_clip_min: f64,
    residual_clip_max: f64,
    _n_threads: i32,
    compute_corrected: bool,
) -> extendr_api::Result<Robj> {
    let view = CscView::from_slots(&x, &i, &p, rows, cols);
    sct_residual_stats_and_corrected_impl(
        view,
        &vec_from_doubles(&theta),
        &vec_from_doubles(&intercept),
        &vec_from_doubles(&slope),
        &vec_from_doubles(&log_umi),
        target_log_umi,
        min_var,
        residual_clip_min,
        residual_clip_max,
        compute_corrected,
    )
}

#[extendr]
fn sct_pearson_residual_matrix(
    x: Doubles,
    i: Integers,
    p: Integers,
    rows: i32,
    cols: i32,
    theta: Doubles,
    intercept: Doubles,
    slope: Doubles,
    log_umi: Doubles,
    feature_index: Integers,
    min_var: Doubles,
    clip_min: f64,
    clip_max: f64,
    do_center: bool,
    _n_threads: i32,
    _display_progress: bool,
) -> extendr_api::Result<RMatrix<f64>> {
    let view = CscView::from_slots(&x, &i, &p, rows, cols);
    sct_pearson_residual_matrix_impl(
        view,
        &vec_from_doubles(&theta),
        &vec_from_doubles(&intercept),
        &vec_from_doubles(&slope),
        &vec_from_doubles(&log_umi),
        &vec_from_integers(&feature_index),
        &vec_from_doubles(&min_var),
        clip_min,
        clip_max,
        do_center,
    )
}

#[extendr]
fn run_modularity_clustering_multi(
    x: Doubles,
    i: Integers,
    p: Integers,
    nrows: i32,
    ncols: i32,
    modularity_function: i32,
    resolutions: Doubles,
    algorithm: i32,
    n_random_starts: i32,
    n_iterations: i32,
    random_seed: i32,
    print_output: bool,
    edgefilename: &str,
    n_threads: i32,
) -> Result<List, extendr_api::Error> {
    let clusters = run_modularity_clustering_multi_impl(
        x.as_robj().as_real_slice().expect("numeric x"),
        i.as_robj().as_integer_slice().expect("integer i"),
        p.as_robj().as_integer_slice().expect("integer p"),
        nrows,
        ncols,
        modularity_function,
        &vec_from_doubles(&resolutions),
        algorithm,
        n_random_starts,
        n_iterations,
        random_seed,
        print_output,
        edgefilename,
        n_threads,
    )
    .map_err(|msg| extendr_api::Error::Other(msg.into()))?;

    let items: Vec<Robj> = clusters
        .into_iter()
        .map(|labels| {
            let mut out = Integers::new(labels.len());
            out.as_robj_mut()
                .as_integer_slice_mut()
                .expect("cluster labels")
                .copy_from_slice(&labels);
            Robj::from(out)
        })
        .collect();
    Ok(List::from_values(items))
}

extendr_module! {
    mod seurust;
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
    fn fast_dense_row_scale;
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
    fn run_modularity_clustering_multi;
    fn sparse_row_mean_var;
    fn sparse_row_mean_var_legacy;
    fn sparse_row_var_std_legacy_rows;
    fn find_all_markers_sparse_fold_change_stats;
    fn cca_crossprod_multiply;
    fn count_anchor_shared_neighbors;
    fn eigen_gram_pca;
    fn sct_residual_stats_and_corrected;
    fn sct_pearson_residual_matrix;
}
