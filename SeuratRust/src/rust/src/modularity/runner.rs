use super::clustering::Clustering;
use super::java_random::JavaRandom;
use super::network::{matrix_to_network, read_input_file, Network};
use super::vos::VOSClusteringTechnique;
use std::sync::Arc;

/// Run modularity clustering; mirrors the former `modularity_bridge.cpp` logic.
///
/// Random starts run sequentially with a shared `JavaRandom` for exact C++ parity.
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
) -> Result<Vec<i32>, String> {
    if modularity_function != 1 && modularity_function != 2 {
        return Err("Modularity parameter must be equal to 1 or 2.".to_string());
    }
    if algorithm != 1 && algorithm != 2 && algorithm != 3 && algorithm != 4 {
        return Err(
            "Algorithm for modularity optimization must be 1, 2, 3, or 4".to_string(),
        );
    }
    if n_random_starts < 1 {
        return Err("Have to have at least one start".to_string());
    }
    if n_iterations < 1 {
        return Err("Need at least one interation".to_string());
    }
    if modularity_function == 2 && resolution > 1.0 {
        return Err("error: resolution<1 for alternative modularity".to_string());
    }

    let network = Arc::new(if !edge_filename.is_empty() {
        read_input_file(edge_filename, modularity_function)?
    } else {
        build_network_from_snn(x, i, p, nrows, ncols, modularity_function)?
    });

    let resolution2 = if modularity_function == 1 {
        resolution
            / (2.0 * network.get_total_edge_weight() + network.total_edge_weight_self_links)
    } else {
        resolution
    };

    let mut best_clustering: Option<Clustering> = None;
    let mut max_modularity = f64::NEG_INFINITY;
    let mut random = JavaRandom::new(random_seed as u64);

    for _start in 0..n_random_starts {
        let mut vos = VOSClusteringTechnique::new(Arc::clone(&network), resolution2);
        let mut j = 0;
        let mut update = true;
        let mut modularity = 0.0;

        while j < n_iterations && update {
            match algorithm {
                1 => update = vos.run_louvain_algorithm(&mut random),
                2 => update = vos.run_louvain_algorithm_with_multilevel_refinement(&mut random),
                3 => {
                    vos.run_smart_local_moving_algorithm(&mut random);
                }
                _ => {}
            }
            j += 1;
            modularity = vos.calc_quality_function();
        }

        if modularity > max_modularity {
            best_clustering = Some(vos.clustering);
            max_modularity = modularity;
        }
    }

    let mut clustering = best_clustering.ok_or_else(|| "Clustering step failed.".to_string())?;
    clustering.order_clusters_by_n_nodes();
    Ok(clustering.cluster)
}

fn build_network_from_snn(
    snn_x: &[f64],
    snn_i: &[i32],
    snn_p: &[i32],
    snn_nrows: i32,
    snn_ncols: i32,
    modularity_function: i32,
) -> Result<Network, String> {
    let ncols = snn_ncols;
    let mut node1 = Vec::with_capacity(snn_x.len());
    let mut node2 = Vec::with_capacity(snn_x.len());
    let mut edgeweights = Vec::with_capacity(snn_x.len());

    for col in 0..ncols {
        let start = snn_p[col as usize] as usize;
        let end = snn_p[(col + 1) as usize] as usize;
        for idx in start..end {
            let row = snn_i[idx];
            if col >= row {
                continue;
            }
            node1.push(col);
            node2.push(row);
            edgeweights.push(snn_x[idx]);
        }
    }

    if node1.is_empty() {
        return Err("Matrix contained no network data.  Check format.".to_string());
    }

    let n_nodes = snn_nrows.max(snn_ncols);
    Ok(matrix_to_network(
        &node1,
        &node2,
        &edgeweights,
        modularity_function,
        n_nodes,
    ))
}
