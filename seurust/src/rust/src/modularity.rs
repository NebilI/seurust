//! Pure Rust modularity optimizer (ported from ModularityOptimizer C++).

#[path = "modularity/clustering.rs"]
mod clustering;
#[path = "modularity/java_random.rs"]
mod java_random;
#[path = "modularity/network.rs"]
mod network;
#[path = "modularity/runner.rs"]
mod runner;
#[path = "modularity/vos.rs"]
mod vos;

use runner::{cluster_network, load_network, run_modularity_clustering};

pub fn run_modularity_clustering_impl(
    x: &[f64],
    i: &[i32],
    p: &[i32],
    nrows: i32,
    ncols: i32,
    modularity_function: i32,
    resolution: f64,
    algorithm: i32,
    n_random_starts: i32,
    n_iterations: i32,
    random_seed: i32,
    print_output: bool,
    edge_filename: &str,
    n_threads: i32,
) -> Result<Vec<i32>, String> {
    run_modularity_clustering(
        x,
        i,
        p,
        nrows,
        ncols,
        modularity_function,
        resolution,
        algorithm,
        n_random_starts,
        n_iterations,
        random_seed,
        print_output,
        edge_filename,
        n_threads,
    )
}

pub fn run_modularity_clustering_multi_impl(
    x: &[f64],
    i: &[i32],
    p: &[i32],
    nrows: i32,
    ncols: i32,
    modularity_function: i32,
    resolutions: &[f64],
    algorithm: i32,
    n_random_starts: i32,
    n_iterations: i32,
    random_seed: i32,
    print_output: bool,
    edge_filename: &str,
    n_threads: i32,
) -> Result<Vec<Vec<i32>>, String> {
    let _ = print_output;
    let network = load_network(
        x,
        i,
        p,
        nrows,
        ncols,
        modularity_function,
        algorithm,
        n_random_starts,
        n_iterations,
        edge_filename,
    )?;
    let mut results = Vec::with_capacity(resolutions.len());
    for &resolution in resolutions {
        results.push(cluster_network(
            std::sync::Arc::clone(&network),
            modularity_function,
            resolution,
            algorithm,
            n_random_starts,
            n_iterations,
            random_seed,
            n_threads,
        )?);
    }
    Ok(results)
}
