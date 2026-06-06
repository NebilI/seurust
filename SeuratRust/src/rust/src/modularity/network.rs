use super::clustering::Clustering;

pub struct Network {
    pub n_nodes: i32,
    pub n_edges: i32,
    pub node_weight: Vec<f64>,
    pub first_neighbor_index: Vec<i32>,
    pub neighbor: Vec<i32>,
    pub edge_weight: Vec<f64>,
    pub total_edge_weight_self_links: f64,
}

impl Network {
    pub fn new() -> Self {
        Self {
            n_nodes: 0,
            n_edges: 0,
            node_weight: Vec::new(),
            first_neighbor_index: Vec::new(),
            neighbor: Vec::new(),
            edge_weight: Vec::new(),
            total_edge_weight_self_links: 0.0,
        }
    }

    pub fn from_adjacency(
        n_nodes: i32,
        node_weight: Option<&[f64]>,
        first_neighbor_index: Vec<i32>,
        neighbor: Vec<i32>,
        edge_weight: Option<&[f64]>,
    ) -> Self {
        let n_edges = neighbor.len() as i32;
        let mut edge_weight_vec = vec![1.0; n_edges as usize];
        if let Some(ew) = edge_weight {
            edge_weight_vec.copy_from_slice(ew);
        }

        let node_weight_vec = if let Some(nw) = node_weight {
            nw.to_vec()
        } else {
            let mut net = Self {
                n_nodes,
                n_edges,
                node_weight: Vec::new(),
                first_neighbor_index: first_neighbor_index.clone(),
                neighbor: neighbor.clone(),
                edge_weight: edge_weight_vec.clone(),
                total_edge_weight_self_links: 0.0,
            };
            net.get_total_edge_weight_per_node()
        };

        Self {
            n_nodes,
            n_edges,
            node_weight: node_weight_vec,
            first_neighbor_index,
            neighbor,
            edge_weight: edge_weight_vec,
            total_edge_weight_self_links: 0.0,
        }
    }

    pub fn from_edge_list(
        n_nodes: i32,
        node_weight: Option<&[f64]>,
        edge: &[Vec<i32>],
        edge_weight: Option<&[f64]>,
    ) -> Self {
        if edge.len() != 2 || edge[0].len() != edge[1].len() {
            panic!("Edge was supposed to be an array with 2 columns of equal size.");
        }

        let mut neighbor = vec![0; edge[0].len()];
        let mut edge_weight2 = vec![0.0; edge[0].len()];
        let mut first_neighbor_index = vec![0; (n_nodes + 1) as usize];
        let mut n_edges = 0i32;
        let mut total_edge_weight_self_links = 0.0;
        let mut i = 1i32;

        for j in 0..edge[0].len() {
            if edge[0][j] != edge[1][j] {
                if edge[0][j] >= i {
                    while i <= edge[0][j] {
                        first_neighbor_index[i as usize] = n_edges;
                        i += 1;
                    }
                }
                neighbor[n_edges as usize] = edge[1][j];
                edge_weight2[n_edges as usize] =
                    edge_weight.map(|ew| ew[j]).unwrap_or(1.0);
                n_edges += 1;
            } else {
                total_edge_weight_self_links +=
                    edge_weight.map(|ew| ew[j]).unwrap_or(1.0);
            }
        }
        for idx in i..=n_nodes {
            first_neighbor_index[idx as usize] = n_edges;
        }

        neighbor.truncate(n_edges as usize);
        edge_weight2.truncate(n_edges as usize);

        let node_weight_vec = if let Some(nw) = node_weight {
            nw.to_vec()
        } else {
            let mut net = Self {
                n_nodes,
                n_edges,
                node_weight: Vec::new(),
                first_neighbor_index: first_neighbor_index.clone(),
                neighbor: neighbor.clone(),
                edge_weight: edge_weight2.clone(),
                total_edge_weight_self_links,
            };
            net.get_total_edge_weight_per_node()
        };

        Self {
            n_nodes,
            n_edges,
            node_weight: node_weight_vec,
            first_neighbor_index,
            neighbor,
            edge_weight: edge_weight2,
            total_edge_weight_self_links,
        }
    }

    pub fn n_nodes(&self) -> i32 {
        self.n_nodes
    }

    pub fn n_edges(&self) -> i32 {
        self.n_edges / 2
    }

    pub fn get_total_node_weight(&self) -> f64 {
        self.node_weight.iter().sum()
    }

    pub fn get_total_edge_weight(&self) -> f64 {
        self.edge_weight.iter().sum::<f64>() / 2.0
    }

    pub fn get_total_edge_weight_node(&self, node: i32) -> f64 {
        let start = self.first_neighbor_index[node as usize] as usize;
        let end = self.first_neighbor_index[(node + 1) as usize] as usize;
        self.edge_weight[start..end].iter().sum()
    }

    pub fn get_total_edge_weight_per_node(&self) -> Vec<f64> {
        (0..self.n_nodes)
            .map(|i| self.get_total_edge_weight_node(i))
            .collect()
    }

    pub fn create_subnetworks(&self, clustering: &Clustering) -> Vec<Network> {
        let node_per_cluster = clustering.get_nodes_per_cluster();
        let mut subnetwork_node = vec![0; self.n_nodes as usize];
        let mut subnetwork_neighbor = vec![0; self.n_edges as usize];
        let mut subnetwork_edge_weight = vec![0.0; self.n_edges as usize];

        node_per_cluster
            .iter()
            .enumerate()
            .map(|(i, nodes)| {
                self.create_subnetwork(
                    clustering,
                    i as i32,
                    nodes,
                    &mut subnetwork_node,
                    &mut subnetwork_neighbor,
                    &mut subnetwork_edge_weight,
                )
            })
            .collect()
    }

    pub fn create_reduced_network(&self, clustering: &Clustering) -> Network {
        let mut reduced_network = Network::new();
        reduced_network.n_nodes = clustering.n_clusters;
        reduced_network.n_edges = 0;
        reduced_network.node_weight = vec![0.0; clustering.n_clusters as usize];
        reduced_network.first_neighbor_index = vec![0; (clustering.n_clusters + 1) as usize];
        reduced_network.total_edge_weight_self_links = self.total_edge_weight_self_links;

        let mut reduced_network_neighbor1 = vec![0; self.n_edges as usize];
        let mut reduced_network_edge_weight1 = vec![0.0; self.n_edges as usize];
        let mut reduced_network_neighbor2 = vec![0; (clustering.n_clusters - 1).max(0) as usize];
        let mut reduced_network_edge_weight2 = vec![0.0; clustering.n_clusters as usize];

        let node_per_cluster = clustering.get_nodes_per_cluster();

        for i in 0..clustering.n_clusters {
            let mut j = 0i32;
            for &l in &node_per_cluster[i as usize] {
                reduced_network.node_weight[i as usize] += self.node_weight[l as usize];

                let start = self.first_neighbor_index[l as usize] as usize;
                let end = self.first_neighbor_index[(l + 1) as usize] as usize;
                for m in start..end {
                    let n = clustering.cluster[self.neighbor[m] as usize];
                    if n != i {
                        if reduced_network_edge_weight2[n as usize] == 0.0 {
                            reduced_network_neighbor2[j as usize] = n;
                            j += 1;
                        }
                        reduced_network_edge_weight2[n as usize] += self.edge_weight[m];
                    } else {
                        reduced_network.total_edge_weight_self_links += self.edge_weight[m];
                    }
                }
            }

            for k in 0..j {
                let n = reduced_network_neighbor2[k as usize] as usize;
                reduced_network_neighbor1[(reduced_network.n_edges + k) as usize] = n as i32;
                reduced_network_edge_weight1[(reduced_network.n_edges + k) as usize] =
                    reduced_network_edge_weight2[n];
                reduced_network_edge_weight2[n] = 0.0;
            }
            reduced_network.n_edges += j;
            reduced_network.first_neighbor_index[(i + 1) as usize] = reduced_network.n_edges;
        }

        reduced_network.neighbor = reduced_network_neighbor1
            [..reduced_network.n_edges as usize]
            .to_vec();
        reduced_network.edge_weight = reduced_network_edge_weight1
            [..reduced_network.n_edges as usize]
            .to_vec();
        reduced_network
    }

    /// Clone network state for a fresh VOS random start (same topology/weights).
    pub fn clone_for_run(&self) -> Network {
        Network {
            n_nodes: self.n_nodes,
            n_edges: self.n_edges,
            node_weight: self.node_weight.clone(),
            first_neighbor_index: self.first_neighbor_index.clone(),
            neighbor: self.neighbor.clone(),
            edge_weight: self.edge_weight.clone(),
            total_edge_weight_self_links: self.total_edge_weight_self_links,
        }
    }

    fn create_subnetwork(
        &self,
        clustering: &Clustering,
        cluster: i32,
        node: &[i32],
        subnetwork_node: &mut [i32],
        subnetwork_neighbor: &mut [i32],
        subnetwork_edge_weight: &mut [f64],
    ) -> Network {
        let mut subnetwork = Network::new();
        subnetwork.n_nodes = node.len() as i32;

        if subnetwork.n_nodes == 1 {
            subnetwork.n_edges = 0;
            subnetwork.node_weight = vec![self.node_weight[node[0] as usize]];
            subnetwork.first_neighbor_index = vec![0, 0];
            subnetwork.neighbor = Vec::new();
            subnetwork.edge_weight = Vec::new();
        } else {
            for (i, &n) in node.iter().enumerate() {
                subnetwork_node[n as usize] = i as i32;
            }

            subnetwork.n_edges = 0;
            subnetwork.node_weight = vec![0.0; subnetwork.n_nodes as usize];
            subnetwork.first_neighbor_index = vec![0; (subnetwork.n_nodes + 1) as usize];

            for i in 0..subnetwork.n_nodes {
                let j = node[i as usize];
                subnetwork.node_weight[i as usize] = self.node_weight[j as usize];
                let start = self.first_neighbor_index[j as usize] as usize;
                let end = self.first_neighbor_index[(j + 1) as usize] as usize;
                for k in start..end {
                    if clustering.cluster[self.neighbor[k] as usize] == cluster {
                        subnetwork_neighbor[subnetwork.n_edges as usize] =
                            subnetwork_node[self.neighbor[k] as usize];
                        subnetwork_edge_weight[subnetwork.n_edges as usize] =
                            self.edge_weight[k];
                        subnetwork.n_edges += 1;
                    }
                }
                subnetwork.first_neighbor_index[(i + 1) as usize] = subnetwork.n_edges;
            }

            subnetwork.neighbor =
                subnetwork_neighbor[..subnetwork.n_edges as usize].to_vec();
            subnetwork.edge_weight =
                subnetwork_edge_weight[..subnetwork.n_edges as usize].to_vec();
        }

        subnetwork.total_edge_weight_self_links = 0.0;
        subnetwork
    }
}

pub fn matrix_to_network(
    node1: &[i32],
    node2: &[i32],
    edge_weight1: &[f64],
    modularity_function: i32,
    n_nodes: i32,
) -> Network {
    let mut n_neighbors = vec![0; n_nodes as usize];
    for i in 0..node1.len() {
        if node1[i] < node2[i] {
            n_neighbors[node1[i] as usize] += 1;
            n_neighbors[node2[i] as usize] += 1;
        }
    }

    let mut first_neighbor_index = vec![0; (n_nodes + 1) as usize];
    let mut n_edges = 0i32;
    for i in 0..n_nodes {
        first_neighbor_index[i as usize] = n_edges;
        n_edges += n_neighbors[i as usize];
    }
    first_neighbor_index[n_nodes as usize] = n_edges;

    let mut neighbor = vec![0; n_edges as usize];
    let mut edge_weight2 = vec![0.0; n_edges as usize];
    n_neighbors.fill(0);

    for i in 0..node1.len() {
        if node1[i] < node2[i] {
            let mut j =
                first_neighbor_index[node1[i] as usize] + n_neighbors[node1[i] as usize];
            neighbor[j as usize] = node2[i];
            edge_weight2[j as usize] = edge_weight1[i];
            n_neighbors[node1[i] as usize] += 1;

            j = first_neighbor_index[node2[i] as usize] + n_neighbors[node2[i] as usize];
            neighbor[j as usize] = node1[i];
            edge_weight2[j as usize] = edge_weight1[i];
            n_neighbors[node2[i] as usize] += 1;
        }
    }

    if modularity_function == 1 {
        Network::from_adjacency(
            n_nodes,
            None,
            first_neighbor_index,
            neighbor,
            Some(&edge_weight2),
        )
    } else {
        let node_weight = vec![1.0; n_nodes as usize];
        Network::from_adjacency(
            n_nodes,
            Some(&node_weight),
            first_neighbor_index,
            neighbor,
            Some(&edge_weight2),
        )
    }
}

pub fn split(s: &str, delimiter: char) -> Vec<String> {
    s.split(delimiter).map(|t| t.to_string()).collect()
}

pub fn read_input_file(fname: &str, modularity_function: i32) -> Result<Network, String> {
    let content = std::fs::read_to_string(fname)
        .map_err(|_| "File could not be opened.".to_string())?;

    let lines: Vec<&str> = content.lines().collect();
    let n_lines = lines.len();

    let mut node1 = vec![0; n_lines];
    let mut node2 = vec![0; n_lines];
    let mut edge_weight1 = vec![1.0; n_lines];

    for (j, line) in lines.iter().enumerate() {
        let splitted_line = split(line, '\t');
        node1[j] = splitted_line[0]
            .parse()
            .map_err(|e: std::num::ParseIntError| e.to_string())?;
        node2[j] = splitted_line[1]
            .parse()
            .map_err(|e: std::num::ParseIntError| e.to_string())?;
        if splitted_line.len() > 2 {
            edge_weight1[j] = splitted_line[2]
                .parse()
                .map_err(|e: std::num::ParseFloatError| e.to_string())?;
        }
    }

    let n_nodes = node1.iter().chain(node2.iter()).max().unwrap_or(&0) + 1;
    Ok(matrix_to_network(
        &node1,
        &node2,
        &edge_weight1,
        modularity_function,
        n_nodes,
    ))
}
