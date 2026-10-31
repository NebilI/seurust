use super::clustering::Clustering;
use super::java_random::JavaRandom;
use super::network::{network_from_snn_csc, read_input_file, Network};
use super::vos::VOSClusteringTechnique;
use rayon::prelude::*;
use std::sync::Arc;

fn run_single_start(
    network: &Arc<Network>,
    resolution2: f64,
    algorithm: i32,
    n_iterations: i32,
    random: &mut JavaRandom,
) -> (Clustering, f64) {
    let mut vos = VOSClusteringTechnique::new(Arc::clone(network), resolution2);
    let mut j = 0;
    let mut update = true;
    let mut modularity = 0.0;

    while j < n_iterations && update {
        match algorithm {
            1 => update = vos.run_louvain_algorithm(random),
            2 => update = vos.run_louvain_algorithm_with_multilevel_refinement(random),
            3 => {
                vos.run_smart_local_moving_algorithm(random);
            }
            _ => {}
        }
        j += 1;
        modularity = vos.calc_quality_function();
    }

    (vos.clustering, modularity)
}

pub(crate) fn load_network(
    x: &[f64],
    i: &[i32],
    p: &[i32],
    nrows: i32,
    ncols: i32,
    modularity_function: i32,
    algorithm: i32,
    n_random_starts: i32,
    n_iterations: i32,
    edge_filename: &str,
) -> Result<Arc<Network>, String> {
    if modularity_function != 1 && modularity_function != 2 {
        return Err("Modularity parameter must be equal to 1 or 2.".to_string());
    }
    if algorithm != 1 && algorithm != 2 && algorithm != 3 && algorithm != 4 {
        return Err("Algorithm for modularity optimization must be 1, 2, 3, or 4".to_string());
    }
    if n_random_starts < 1 {
        return Err("Have to have at least one start".to_string());
    }
    if n_iterations < 1 {
        return Err("Need at least one interation".to_string());
    }
    Ok(Arc::new(if !edge_filename.is_empty() {
        read_input_file(edge_filename, modularity_function)?
    } else {
        build_network_from_snn(x, i, p, nrows, ncols, modularity_function)?
    }))
}

/// Run modularity clustering on an already-built network (Seurat 5.6 multi-resolution).
pub(crate) fn cluster_network(
    network: Arc<Network>,
    modularity_function: i32,
    resolution: f64,
    algorithm: i32,
    n_random_starts: i32,
    n_iterations: i32,
    random_seed: i32,
    n_threads: i32,
) -> Result<Vec<i32>, String> {
    if modularity_function == 2 && resolution > 1.0 {
        return Err("error: resolution<1 for alternative modularity".to_string());
    }

    let resolution2 = if modularity_function == 1 {
        resolution / (2.0 * network.get_total_edge_weight() + network.total_edge_weight_self_links)
    } else {
        resolution
    };

    let mut best_clustering: Option<Clustering> = None;
    let mut max_modularity = f64::NEG_INFINITY;
    let serial = n_threads <= 1 || n_random_starts <= 1;
    if serial {
        let mut random = JavaRandom::new(random_seed as u64);
        for _start in 0..n_random_starts {
            let (clustering, modularity) =
                run_single_start(&network, resolution2, algorithm, n_iterations, &mut random);

            if modularity > max_modularity {
                best_clustering = Some(clustering);
                max_modularity = modularity;
            }
        }
    } else {
        let starts: Vec<(Clustering, f64)> = (0..n_random_starts)
            .into_par_iter()
            .map(|s| {
                let mut random = JavaRandom::new(random_seed as u64 + s as u64);
                run_single_start(&network, resolution2, algorithm, n_iterations, &mut random)
            })
            .collect();
        for (clustering, modularity) in starts {
            if modularity > max_modularity {
                best_clustering = Some(clustering);
                max_modularity = modularity;
            }
        }
    }

    let mut clustering = best_clustering.ok_or_else(|| "Clustering step failed.".to_string())?;
    clustering.order_clusters_by_n_nodes();
    Ok(clustering.cluster)
}

/// Run modularity clustering; mirrors Seurat 5.6 `RunModularityClusteringCpp`.
///
/// `n_threads <= 1` keeps a shared `JavaRandom` stream across starts (exact
/// serial C++ parity). `n_threads > 1` assigns `seed + start` to each restart.
pub fn run_modularity_clustering(
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
    _print_output: bool,
    edge_filename: &str,
    n_threads: i32,
) -> Result<Vec<i32>, String> {
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
    cluster_network(
        network,
        modularity_function,
        resolution,
        algorithm,
        n_random_starts,
        n_iterations,
        random_seed,
        n_threads,
    )
}

fn build_network_from_snn(
    snn_x: &[f64],
    snn_i: &[i32],
    snn_p: &[i32],
    snn_nrows: i32,
    snn_ncols: i32,
    modularity_function: i32,
) -> Result<Network, String> {
    network_from_snn_csc(
        snn_x,
        snn_i,
        snn_p,
        snn_nrows,
        snn_ncols,
        modularity_function,
    )
}
