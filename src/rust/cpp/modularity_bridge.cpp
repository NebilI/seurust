// C bridge so the Rust/extendr crate can call the existing ModularityOptimizer C++.
// Logic mirrors RModularityOptimizer.cpp without Rcpp dependencies.

#include "ModularityOptimizer.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <string>

using namespace ModularityOptimizer;

static void set_error(char* error_msg, int error_msg_len, const char* msg) {
  if (error_msg == nullptr || error_msg_len <= 0) {
    return;
  }
  std::strncpy(error_msg, msg, static_cast<size_t>(error_msg_len - 1));
  error_msg[error_msg_len - 1] = '\0';
}

extern "C" int* modularity_cluster_from_snn(
    const double* snn_x,
    int snn_x_len,
    const int* snn_i,
    int snn_i_len,
    const int* snn_p,
    int snn_p_len,
    int snn_nrows,
    int snn_ncols,
    int modularity_function,
    double resolution,
    int algorithm,
    int n_random_starts,
    int n_iterations,
    int64_t random_seed,
    int print_output,
    const char* edge_filename,
    int* out_len,
    char* error_msg,
    int error_msg_len) {

  (void)snn_x_len;
  (void)snn_i_len;
  (void)snn_p_len;
  (void)print_output;

  try {
    if (modularity_function != 1 && modularity_function != 2) {
      set_error(error_msg, error_msg_len, "Modularity parameter must be equal to 1 or 2.");
      return nullptr;
    }
    if (algorithm != 1 && algorithm != 2 && algorithm != 3 && algorithm != 4) {
      set_error(error_msg, error_msg_len, "Algorithm for modularity optimization must be 1, 2, 3, or 4");
      return nullptr;
    }
    if (n_random_starts < 1) {
      set_error(error_msg, error_msg_len, "Have to have at least one start");
      return nullptr;
    }
    if (n_iterations < 1) {
      set_error(error_msg, error_msg_len, "Need at least one interation");
      return nullptr;
    }
    if (modularity_function == 2 && resolution > 1.0) {
      set_error(error_msg, error_msg_len, "error: resolution<1 for alternative modularity");
      return nullptr;
    }

    std::shared_ptr<Network> network;
    std::string edgefile = edge_filename == nullptr ? "" : edge_filename;
    if (!edgefile.empty()) {
      network = readInputFile(edgefile, modularity_function);
    } else {
      if (snn_x == nullptr || snn_i == nullptr || snn_p == nullptr) {
        set_error(error_msg, error_msg_len, "Missing SNN matrix slots.");
        return nullptr;
      }

      int ncols = snn_ncols;
      IVector node1;
      IVector node2;
      DVector edgeweights;
      node1.reserve(snn_x_len);
      node2.reserve(snn_x_len);
      edgeweights.reserve(snn_x_len);

      for (int col = 0; col < ncols; ++col) {
        int start = snn_p[col];
        int end = snn_p[col + 1];
        for (int idx = start; idx < end; ++idx) {
          int row = snn_i[idx];
          if (col >= row) {
            continue;
          }
          node1.push_back(col);
          node2.push_back(row);
          edgeweights.push_back(snn_x[idx]);
        }
      }

      if (node1.empty()) {
        set_error(error_msg, error_msg_len, "Matrix contained no network data.  Check format.");
        return nullptr;
      }

      int n_nodes = std::max(snn_nrows, snn_ncols);
      network = matrixToNetwork(node1, node2, edgeweights, modularity_function, n_nodes);
    }

    double resolution2 = (modularity_function == 1)
                             ? (resolution / (2 * network->getTotalEdgeWeight() +
                                              network->getTotalEdgeWeightSelfLinks()))
                             : resolution;

    std::shared_ptr<Clustering> clustering;
    double max_modularity = -std::numeric_limits<double>::infinity();
    JavaRandom random(static_cast<uint64_t>(random_seed));

    for (int start = 0; start < n_random_starts; ++start) {
      VOSClusteringTechnique vos_clustering_technique(network, resolution2);
      int j = 0;
      bool update = true;
      double modularity = 0.0;
      do {
        if (algorithm == 1) {
          update = vos_clustering_technique.runLouvainAlgorithm(random);
        } else if (algorithm == 2) {
          update = vos_clustering_technique.runLouvainAlgorithmWithMultilevelRefinement(random);
        } else if (algorithm == 3) {
          vos_clustering_technique.runSmartLocalMovingAlgorithm(random);
        }
        j++;
        modularity = vos_clustering_technique.calcQualityFunction();
      } while ((j < n_iterations) && update);

      if (modularity > max_modularity) {
        clustering = vos_clustering_technique.getClustering();
        max_modularity = modularity;
      }
    }

    if (clustering == nullptr) {
      set_error(error_msg, error_msg_len, "Clustering step failed.");
      return nullptr;
    }

    clustering->orderClustersByNNodes();
    int* result = new int[clustering->cluster.size()];
    std::copy(clustering->cluster.begin(), clustering->cluster.end(), result);
    *out_len = static_cast<int>(clustering->cluster.size());
    return result;
  } catch (const std::exception& ex) {
    set_error(error_msg, error_msg_len, ex.what());
    return nullptr;
  } catch (...) {
    set_error(error_msg, error_msg_len, "Unknown modularity optimizer error.");
    return nullptr;
  }
}

extern "C" void modularity_free(int* ptr) {
  delete[] ptr;
}
