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

use rayon::prelude::*;
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
    // Each resolution owns a fresh JavaRandom(random_seed) and only reads
    // `network`, so the labels stay bit-identical to a serial loop.
    if resolutions.len() <= 1 {
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
        return Ok(results);
    }
    resolutions
        .par_iter()
        .map(|&resolution| {
            cluster_network(
                std::sync::Arc::clone(&network),
                modularity_function,
                resolution,
                algorithm,
                n_random_starts,
                n_iterations,
                random_seed,
                n_threads,
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn multi_resolutions_match_repeated_single_runs() {
        // 4-cycle stored as the upper triangle of a symmetric CSC.
        let x = [1.0, 1.0, 1.0, 1.0];
        let i = [1, 3, 2, 3];
        let p = [0, 2, 3, 4, 4];
        let resolutions = [0.4, 0.8, 1.2];
        let multi = run_modularity_clustering_multi_impl(
            &x,
            &i,
            &p,
            4,
            4,
            1,
            &resolutions,
            1,
            1,
            5,
            7,
            false,
            "",
            1,
        )
        .unwrap();
        assert_eq!(multi.len(), 3);
        for (idx, &res) in resolutions.iter().enumerate() {
            let single =
                run_modularity_clustering_impl(&x, &i, &p, 4, 4, 1, res, 1, 1, 5, 7, false, "", 1)
                    .unwrap();
            assert_eq!(multi[idx], single);
        }
    }
}
