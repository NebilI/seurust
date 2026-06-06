// C bridge: Eigen ComputeSNN (mirrors src/snn.cpp) for the Rust/extendr crate.

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>

typedef Eigen::Triplet<double> Triplet;

static void set_error(char* error_msg, int error_msg_len, const char* msg) {
  if (error_msg == nullptr || error_msg_len <= 0) {
    return;
  }
  std::strncpy(error_msg, msg, static_cast<size_t>(error_msg_len - 1));
  error_msg[error_msg_len - 1] = '\0';
}

static int fill_snn_csc_prealloc(
    const Eigen::SparseMatrix<double>& snn,
    int k,
    double prune,
    double* out_x,
    int* out_i,
    int* out_p,
    int out_x_len,
    int out_i_len,
    int out_p_len,
    int* out_nnz_required) {
  const double k_f = static_cast<double>(k);
  const int n = static_cast<int>(snn.cols());

  if (out_p_len < n + 1) {
    return -2;
  }

  int nnz = 0;
  for (int col = 0; col < snn.outerSize(); ++col) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(snn, col); it; ++it) {
      const double scaled = it.value() / (k_f + (k_f - it.value()));
      if (scaled >= prune) {
        ++nnz;
      }
    }
  }

  if (out_nnz_required != nullptr) {
    *out_nnz_required = nnz;
  }

  if (out_x_len < nnz || out_i_len < nnz) {
    return -3;
  }

  int nz = 0;
  for (int col = 0; col < n; ++col) {
    out_p[col] = nz;
    for (Eigen::SparseMatrix<double>::InnerIterator it(snn, col); it; ++it) {
      const double scaled = it.value() / (k_f + (k_f - it.value()));
      if (scaled >= prune) {
        out_i[nz] = static_cast<int>(it.row());
        out_x[nz] = scaled;
        ++nz;
      }
    }
  }
  out_p[n] = nz;
  return nnz;
}

static int fill_snn_csc_malloc(
    const Eigen::SparseMatrix<double>& snn,
    int k,
    double prune,
    double** out_x,
    int** out_i,
    int** out_p,
    int* out_nnz) {
  const double k_f = static_cast<double>(k);
  const int n = static_cast<int>(snn.cols());

  int nnz = 0;
  for (int col = 0; col < snn.outerSize(); ++col) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(snn, col); it; ++it) {
      const double scaled = it.value() / (k_f + (k_f - it.value()));
      if (scaled >= prune) {
        ++nnz;
      }
    }
  }

  double* x = static_cast<double*>(std::malloc(sizeof(double) * static_cast<size_t>(nnz)));
  int* i = static_cast<int*>(std::malloc(sizeof(int) * static_cast<size_t>(nnz)));
  int* p = static_cast<int*>(std::malloc(sizeof(int) * static_cast<size_t>(n + 1)));
  if (x == nullptr || i == nullptr || p == nullptr) {
    std::free(x);
    std::free(i);
    std::free(p);
    return -2;
  }

  int nz = 0;
  for (int col = 0; col < n; ++col) {
    p[col] = nz;
    for (Eigen::SparseMatrix<double>::InnerIterator it(snn, col); it; ++it) {
      const double scaled = it.value() / (k_f + (k_f - it.value()));
      if (scaled >= prune) {
        i[nz] = static_cast<int>(it.row());
        x[nz] = scaled;
        ++nz;
      }
    }
  }
  p[n] = nz;

  *out_x = x;
  *out_i = i;
  *out_p = p;
  *out_nnz = nnz;
  return 0;
}

extern "C" int compute_snn_csc(
    const double* nn_ranked,
    int nrows,
    int ncols,
    double prune,
    double** out_x,
    int** out_i,
    int** out_p,
    int* out_nnz,
    char* error_msg,
    int error_msg_len) {
  if (out_x == nullptr || out_i == nullptr || out_p == nullptr || out_nnz == nullptr) {
    set_error(error_msg, error_msg_len, "Null output pointer.");
    return -1;
  }
  if (nn_ranked == nullptr || nrows < 1 || ncols < 1) {
    set_error(error_msg, error_msg_len, "Invalid nn_ranked input.");
    return -1;
  }

  try {
    const int k = ncols;
    std::vector<Triplet> triplet_list;
    triplet_list.reserve(static_cast<size_t>(nrows) * static_cast<size_t>(ncols));

    for (int col = 0; col < ncols; ++col) {
      for (int row = 0; row < nrows; ++row) {
        const double neighbor = nn_ranked[col * nrows + row];
        triplet_list.emplace_back(row, static_cast<int>(neighbor) - 1, 1.0);
      }
    }

    Eigen::SparseMatrix<double> snn(nrows, nrows);
    snn.setFromTriplets(triplet_list.begin(), triplet_list.end());
    snn = snn * snn.transpose();

    return fill_snn_csc_malloc(snn, k, prune, out_x, out_i, out_p, out_nnz);
  } catch (const std::exception& ex) {
    set_error(error_msg, error_msg_len, ex.what());
    return -1;
  } catch (...) {
    set_error(error_msg, error_msg_len, "Unknown error in compute_snn_csc.");
    return -1;
  }
}

extern "C" void compute_snn_csc_free(double* x, int* i, int* p) {
  std::free(x);
  std::free(i);
  std::free(p);
}

extern "C" int compute_snn_csc_into(
    const double* nn_ranked,
    int nrows,
    int ncols,
    double prune,
    double* out_x,
    int out_x_len,
    int* out_i,
    int out_i_len,
    int* out_p,
    int out_p_len,
    int* out_nnz_required,
    char* error_msg,
    int error_msg_len) {
  if (out_x == nullptr || out_i == nullptr || out_p == nullptr) {
    set_error(error_msg, error_msg_len, "Null output pointer.");
    return -1;
  }
  if (nn_ranked == nullptr || nrows < 1 || ncols < 1) {
    set_error(error_msg, error_msg_len, "Invalid nn_ranked input.");
    return -1;
  }

  try {
    const int k = ncols;
    std::vector<Triplet> triplet_list;
    triplet_list.reserve(static_cast<size_t>(nrows) * static_cast<size_t>(ncols));

    for (int col = 0; col < ncols; ++col) {
      for (int row = 0; row < nrows; ++row) {
        const double neighbor = nn_ranked[col * nrows + row];
        triplet_list.emplace_back(row, static_cast<int>(neighbor) - 1, 1.0);
      }
    }

    Eigen::SparseMatrix<double> snn(nrows, nrows);
    snn.setFromTriplets(triplet_list.begin(), triplet_list.end());
    snn = snn * snn.transpose();

    return fill_snn_csc_prealloc(
        snn,
        k,
        prune,
        out_x,
        out_i,
        out_p,
        out_x_len,
        out_i_len,
        out_p_len,
        out_nnz_required);
  } catch (const std::exception& ex) {
    set_error(error_msg, error_msg_len, ex.what());
    return -1;
  } catch (...) {
    set_error(error_msg, error_msg_len, "Unknown error in compute_snn_csc_into.");
    return -1;
  }
}
