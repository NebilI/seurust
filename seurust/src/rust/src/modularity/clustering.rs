pub struct Clustering {
    n_nodes: i32,
    pub n_clusters: i32,
    pub cluster: Vec<i32>,
}

impl Clustering {
    pub fn new(n_nodes: i32) -> Self {
        Self {
            n_nodes,
            n_clusters: 1,
            cluster: vec![0; n_nodes as usize],
        }
    }

    pub fn from_cluster_vec(cluster: Vec<i32>) -> Self {
        let n_nodes = cluster.len() as i32;
        let n_clusters = *cluster.iter().max().unwrap_or(&0) + 1;
        Self {
            n_nodes,
            n_clusters,
            cluster,
        }
    }

    pub fn n_nodes(&self) -> i32 {
        self.n_nodes
    }

    pub fn n_clusters(&self) -> i32 {
        self.n_clusters
    }

    pub fn get_n_nodes_per_cluster(&self) -> Vec<i32> {
        let mut n_nodes_per_cluster = vec![0; self.n_clusters as usize];
        for &clust in &self.cluster {
            n_nodes_per_cluster[clust as usize] += 1;
        }
        n_nodes_per_cluster
    }

    pub fn get_nodes_per_cluster(&self) -> Vec<Vec<i32>> {
        let n_nodes_per_cluster = self.get_n_nodes_per_cluster();
        let mut nodes_per_cluster: Vec<Vec<i32>> = n_nodes_per_cluster
            .iter()
            .map(|&cnt| Vec::with_capacity(cnt as usize))
            .collect();
        for i in 0..self.n_nodes {
            nodes_per_cluster[self.cluster[i as usize] as usize].push(i);
        }
        nodes_per_cluster
    }

    pub fn set_cluster(&mut self, node: i32, cluster: i32) {
        self.cluster[node as usize] = cluster;
        self.n_clusters = self.n_clusters.max(cluster + 1);
    }

    pub fn init_singleton_clusters(&mut self) {
        for i in 0..self.n_nodes {
            self.cluster[i as usize] = i;
        }
        self.n_clusters = self.n_nodes;
    }

    pub fn order_clusters_by_n_nodes(&mut self) {
        let n_nodes_per_cluster = self.get_n_nodes_per_cluster();
        let mut cluster_n_nodes: Vec<(i32, i32)> = (0..self.n_clusters)
            .map(|i| (n_nodes_per_cluster[i as usize], i))
            .collect();

        // Rust's sort is stable, matching C++ stable_sort with descending node count.
        cluster_n_nodes.sort_by(|a, b| b.0.cmp(&a.0));

        let mut new_cluster = vec![0; self.n_clusters as usize];
        let mut i = 0;
        loop {
            new_cluster[cluster_n_nodes[i as usize].1 as usize] = i;
            i += 1;
            if i >= self.n_clusters || cluster_n_nodes[i as usize].0 <= 0 {
                break;
            }
        }
        self.n_clusters = i;
        for node in 0..self.n_nodes as usize {
            self.cluster[node] = new_cluster[self.cluster[node] as usize];
        }
    }

    pub fn merge_clusters(&mut self, clustering: &Clustering) {
        for i in 0..self.n_nodes as usize {
            self.cluster[i] = clustering.cluster[self.cluster[i] as usize];
        }
        self.n_clusters = clustering.n_clusters;
    }
}
