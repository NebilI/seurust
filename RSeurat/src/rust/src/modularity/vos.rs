use super::clustering::Clustering;
use super::java_random::{generate_random_permutation, JavaRandom};
use super::network::Network;
use std::sync::Arc;

pub struct VOSClusteringTechnique {
    pub network: Arc<Network>,
    pub clustering: Clustering,
    pub resolution: f64,
}

impl VOSClusteringTechnique {
    pub fn new(network: Arc<Network>, resolution: f64) -> Self {
        let n_nodes = network.n_nodes();
        let mut clustering = Clustering::new(n_nodes);
        clustering.init_singleton_clusters();
        Self {
            network,
            clustering,
            resolution,
        }
    }

    pub fn with_clustering(network: Arc<Network>, clustering: Clustering, resolution: f64) -> Self {
        Self {
            network,
            clustering,
            resolution,
        }
    }

    pub fn calc_quality_function(&self) -> f64 {
        let mut quality_function = 0.0;
        for i in 0..self.network.n_nodes {
            let j = self.clustering.cluster[i as usize];
            let start = self.network.first_neighbor_index[i as usize] as usize;
            let end = self.network.first_neighbor_index[(i + 1) as usize] as usize;
            for k in start..end {
                if self.clustering.cluster[self.network.neighbor[k] as usize] == j {
                    quality_function += self.network.edge_weight[k];
                }
            }
        }
        quality_function += self.network.total_edge_weight_self_links;

        let mut cluster_weight = vec![0.0; self.clustering.n_clusters as usize];
        for i in 0..self.network.n_nodes as usize {
            cluster_weight[self.clustering.cluster[i] as usize] += self.network.node_weight[i];
        }
        for i in 0..self.clustering.n_clusters as usize {
            quality_function -= cluster_weight[i] * cluster_weight[i] * self.resolution;
        }

        quality_function
            / (2.0 * self.network.get_total_edge_weight()
                + self.network.total_edge_weight_self_links)
    }

    pub fn run_local_moving_algorithm(&mut self, random: &mut JavaRandom) -> bool {
        let mut update = false;
        let n_nodes = self.network.n_nodes;
        if n_nodes == 1 {
            return false;
        }

        let mut cluster_weight = vec![0.0; n_nodes as usize];
        let mut n_nodes_per_cluster = vec![0; n_nodes as usize];

        for i in 0..n_nodes {
            let c = self.clustering.cluster[i as usize] as usize;
            cluster_weight[c] += self.network.node_weight[i as usize];
            n_nodes_per_cluster[c] += 1;
        }

        let mut n_unused_clusters = 0i32;
        let mut unused_cluster = vec![0; n_nodes as usize];
        for i in 0..n_nodes {
            if n_nodes_per_cluster[i as usize] == 0 {
                unused_cluster[n_unused_clusters as usize] = i;
                n_unused_clusters += 1;
            }
        }

        let node_permutation = generate_random_permutation(n_nodes, random);
        let mut edge_weight_per_cluster = vec![0.0; n_nodes as usize];
        let mut neighboring_cluster = vec![0; (n_nodes - 1) as usize];
        let mut n_stable_nodes = 0i32;
        let mut i = 0i32;

        loop {
            let j = node_permutation[i as usize];
            let mut n_neighboring_clusters = 0i32;

            let start = self.network.first_neighbor_index[j as usize] as usize;
            let end = self.network.first_neighbor_index[(j + 1) as usize] as usize;
            for k in start..end {
                let l = self.clustering.cluster[self.network.neighbor[k] as usize];
                if edge_weight_per_cluster[l as usize] == 0.0 {
                    neighboring_cluster[n_neighboring_clusters as usize] = l;
                    n_neighboring_clusters += 1;
                }
                edge_weight_per_cluster[l as usize] += self.network.edge_weight[k];
            }

            let old_cluster = self.clustering.cluster[j as usize] as usize;
            cluster_weight[old_cluster] -= self.network.node_weight[j as usize];
            n_nodes_per_cluster[old_cluster] -= 1;
            if n_nodes_per_cluster[old_cluster] == 0 {
                unused_cluster[n_unused_clusters as usize] = old_cluster as i32;
                n_unused_clusters += 1;
            }

            let mut best_cluster = -1i32;
            let mut max_quality_function = 0.0;
            for k in 0..n_neighboring_clusters {
                let l = neighboring_cluster[k as usize];
                let quality_function = edge_weight_per_cluster[l as usize]
                    - self.network.node_weight[j as usize]
                        * cluster_weight[l as usize]
                        * self.resolution;
                if (quality_function > max_quality_function)
                    || ((quality_function == max_quality_function) && (l < best_cluster))
                {
                    best_cluster = l;
                    max_quality_function = quality_function;
                }
                edge_weight_per_cluster[l as usize] = 0.0;
            }

            if max_quality_function == 0.0 {
                best_cluster = unused_cluster[(n_unused_clusters - 1) as usize];
                n_unused_clusters -= 1;
            }

            let best = best_cluster as usize;
            cluster_weight[best] += self.network.node_weight[j as usize];
            n_nodes_per_cluster[best] += 1;
            if best_cluster == self.clustering.cluster[j as usize] {
                n_stable_nodes += 1;
            } else {
                self.clustering.cluster[j as usize] = best_cluster;
                n_stable_nodes = 1;
                update = true;
            }

            i = if i < n_nodes - 1 { i + 1 } else { 0 };

            if n_stable_nodes >= n_nodes {
                break;
            }
        }

        let mut new_cluster = vec![0; n_nodes as usize];
        self.clustering.n_clusters = 0;
        for i in 0..n_nodes as usize {
            if n_nodes_per_cluster[i] > 0 {
                new_cluster[i] = self.clustering.n_clusters;
                self.clustering.n_clusters += 1;
            }
        }
        for i in 0..n_nodes as usize {
            self.clustering.cluster[i] = new_cluster[self.clustering.cluster[i] as usize];
        }

        update
    }

    pub fn run_louvain_algorithm(&mut self, random: &mut JavaRandom) -> bool {
        if self.network.n_nodes == 1 {
            return false;
        }
        let mut update = self.run_local_moving_algorithm(random);
        if self.clustering.n_clusters < self.network.n_nodes {
            let reduced = self.network.create_reduced_network(&self.clustering);
            let mut vos = VOSClusteringTechnique::new(Arc::new(reduced), self.resolution);
            let update2 = vos.run_louvain_algorithm(random);
            if update2 {
                update = true;
                self.clustering.merge_clusters(&vos.clustering);
            }
        }
        update
    }

    pub fn run_iterated_louvain_algorithm(
        &mut self,
        max_n_iterations: i32,
        random: &mut JavaRandom,
    ) -> bool {
        let mut update;
        let mut i = 0;
        loop {
            update = self.run_louvain_algorithm(random);
            i += 1;
            if i >= max_n_iterations || !update {
                break;
            }
        }
        (i > 1) || update
    }

    pub fn run_louvain_algorithm_with_multilevel_refinement(
        &mut self,
        random: &mut JavaRandom,
    ) -> bool {
        if self.network.n_nodes == 1 {
            return false;
        }

        let mut update = self.run_local_moving_algorithm(random);
        if self.clustering.n_clusters < self.network.n_nodes {
            let reduced = self.network.create_reduced_network(&self.clustering);
            let mut vos = VOSClusteringTechnique::new(Arc::new(reduced), self.resolution);
            let update2 = vos.run_louvain_algorithm_with_multilevel_refinement(random);
            if update2 {
                update = true;
                self.clustering.merge_clusters(&vos.clustering);
                self.run_local_moving_algorithm(random);
            }
        }
        update
    }

    pub fn run_iterated_louvain_algorithm_with_multilevel_refinement(
        &mut self,
        max_n_iterations: i32,
        random: &mut JavaRandom,
    ) -> bool {
        let mut update;
        let mut i = 0;
        loop {
            update = self.run_louvain_algorithm_with_multilevel_refinement(random);
            i += 1;
            if i >= max_n_iterations || !update {
                break;
            }
        }
        (i > 1) || update
    }

    pub fn run_smart_local_moving_algorithm(&mut self, random: &mut JavaRandom) -> bool {
        if self.network.n_nodes == 1 {
            return false;
        }

        let mut update = self.run_local_moving_algorithm(random);
        if self.clustering.n_clusters < self.network.n_nodes {
            let subnetworks = self.network.create_subnetworks(&self.clustering);
            let node_per_cluster = self.clustering.get_nodes_per_cluster();
            self.clustering.n_clusters = 0;
            let mut n_nodes_per_cluster_reduced_network = vec![0; subnetworks.len()];

            for (idx, sub) in subnetworks.into_iter().enumerate() {
                let mut vos = VOSClusteringTechnique::new(Arc::new(sub), self.resolution);
                vos.run_local_moving_algorithm(random);
                for j in 0..vos.network.n_nodes {
                    self.clustering.cluster[node_per_cluster[idx][j as usize] as usize] =
                        self.clustering.n_clusters + vos.clustering.cluster[j as usize];
                }
                self.clustering.n_clusters += vos.clustering.n_clusters;
                n_nodes_per_cluster_reduced_network[idx] = vos.clustering.n_clusters;
            }

            let reduced = self.network.create_reduced_network(&self.clustering);
            let mut vos2 = VOSClusteringTechnique::new(Arc::new(reduced), self.resolution);

            let mut idx = 0i32;
            for (j, &n) in n_nodes_per_cluster_reduced_network.iter().enumerate() {
                for _ in 0..n {
                    vos2.clustering.cluster[idx as usize] = j as i32;
                    idx += 1;
                }
            }
            vos2.clustering.n_clusters = n_nodes_per_cluster_reduced_network.len() as i32;

            update |= vos2.run_smart_local_moving_algorithm(random);
            self.clustering.merge_clusters(&vos2.clustering);
        }
        update
    }

    pub fn run_iterated_smart_local_moving_algorithm(
        &mut self,
        n_iterations: i32,
        random: &mut JavaRandom,
    ) -> bool {
        let mut update = false;
        for _ in 0..n_iterations {
            update |= self.run_smart_local_moving_algorithm(random);
        }
        update
    }

    pub fn remove_cluster(&mut self, cluster: i32) -> i32 {
        let mut cluster_weight = vec![0.0; self.clustering.n_clusters as usize];
        let mut total_edge_weight_per_cluster = vec![0.0; self.clustering.n_clusters as usize];

        for i in 0..self.network.n_nodes {
            let iu = i as usize;
            cluster_weight[self.clustering.cluster[iu] as usize] += self.network.node_weight[iu];
            if self.clustering.cluster[iu] == cluster {
                let start = self.network.first_neighbor_index[i as usize] as usize;
                let end = self.network.first_neighbor_index[(i + 1) as usize] as usize;
                for j in start..end {
                    total_edge_weight_per_cluster
                        [self.clustering.cluster[self.network.neighbor[j] as usize] as usize] +=
                        self.network.edge_weight[j];
                }
            }
        }

        let mut best = -1i32;
        let mut max_quality_function = 0.0;
        for j in 0..self.clustering.n_clusters {
            if (j != cluster) && (cluster_weight[j as usize] > 0.0) {
                let quality_function =
                    total_edge_weight_per_cluster[j as usize] / cluster_weight[j as usize];
                if quality_function > max_quality_function {
                    best = j;
                    max_quality_function = quality_function;
                }
            }
        }

        if best >= 0 {
            for j in 0..self.network.n_nodes {
                let ju = j as usize;
                if self.clustering.cluster[ju] == cluster {
                    self.clustering.cluster[ju] = best;
                }
            }
            if cluster == self.clustering.n_clusters - 1 {
                self.clustering.n_clusters =
                    *self.clustering.cluster.iter().max().unwrap_or(&0) + 1;
            }
        }
        best
    }

    pub fn remove_small_clusters(&mut self, min_n_nodes_per_cluster: i32) {
        let reduced = self.network.create_reduced_network(&self.clustering);
        let mut vos = VOSClusteringTechnique::new(Arc::new(reduced), self.resolution);
        let mut n_nodes_per_cluster = self.clustering.get_n_nodes_per_cluster();

        loop {
            let mut best = -1i32;
            let mut min_count = min_n_nodes_per_cluster;
            for k in 0..vos.clustering.n_clusters {
                if (n_nodes_per_cluster[k as usize] > 0)
                    && (n_nodes_per_cluster[k as usize] < min_count)
                {
                    best = k;
                    min_count = n_nodes_per_cluster[k as usize];
                }
            }

            if best < 0 {
                break;
            }

            let merged = vos.remove_cluster(best);
            if merged >= 0 {
                n_nodes_per_cluster[merged as usize] += n_nodes_per_cluster[best as usize];
            }
            n_nodes_per_cluster[best as usize] = 0;
        }

        self.clustering.merge_clusters(&vos.clustering);
    }
}
