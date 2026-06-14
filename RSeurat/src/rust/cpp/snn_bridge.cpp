// C bridge: Eigen ComputeSNN (mirrors src/snn.cpp) for the Rust/extendr crate.

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <cstdlib>
#include <cstring>
#include <vector>

#ifdef SNN_BRIDGE_RCPP
#include <Rcpp.h>
#include <RcppEigen.h>
#endif

typedef Eigen::Triplet<double> Triplet;

static void set_error(char* error_msg, int error_msg_len, const char* msg) {
  if (error_msg == nullptr || error_msg_len <= 0) {
    return;
  }
  std::strncpy(error_msg, msg, static_cast<size_t>(error_msg_len - 1));
  error_msg[error_msg_len - 1] = '\0';
}

static Eigen::SparseMatrix<double> build_neighbor_matrix(
    const double* nn_ranked,
    int nrows,
    int ncols) {
  const int n = nrows;
  const int k = ncols;

  std::vector<Triplet> triplets;
  triplets.reserve(static_cast<size_t>(n) * static_cast<size_t>(k));
  for (int rank = 0; rank < k; ++rank) {
    const int base = rank * n;
    for (int row = 0; row < n; ++row) {
      triplets.emplace_back(
          row,
          static_cast<int>(nn_ranked[base + row]) - 1,
          1.0);
    }
  }

  Eigen::SparseMatrix<double> neighbor(n, n);
  neighbor.setFromTriplets(triplets.begin(), triplets.end());
  return neighbor;
}

static Eigen::SparseMatrix<double> build_snn_matrix(
    const double* nn_ranked,
    int nrows,
    int ncols) {
  const Eigen::SparseMatrix<double> neighbor =
      build_neighbor_matrix(nn_ranked, nrows, ncols);
  return neighbor * neighbor.transpose();
}

static void scale_prune_snn(Eigen::SparseMatrix<double>& snn, int k, double prune) {
  const double k_f = static_cast<double>(k);
  for (int col = 0; col < snn.outerSize(); ++col) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(snn, col); it; ++it) {
      const double scaled = it.value() / (k_f + (k_f - it.value()));
      if (scaled < prune) {
        it.valueRef() = 0.0;
      } else {
        it.valueRef() = scaled;
      }
    }
  }
  snn.prune(0.0);
}

#ifdef SNN_BRIDGE_RCPP
static Eigen::SparseMatrix<double> compute_snn_from_matrix(
    const Eigen::Ref<const Eigen::MatrixXd>& nn_ranked,
    double prune) {
  const int k = static_cast<int>(nn_ranked.cols());
  Eigen::SparseMatrix<double> snn = build_snn_matrix(
      nn_ranked.data(),
      static_cast<int>(nn_ranked.rows()),
      k);
  scale_prune_snn(snn, k, prune);
  return snn;
}

extern "C" SEXP compute_snn_rcpp_fast(
    const double* nn_ranked,
    int nrows,
    int ncols,
    double prune) {
  try {
    Eigen::SparseMatrix<double> snn = build_snn_matrix(nn_ranked, nrows, ncols);
    scale_prune_snn(snn, ncols, prune);
    return Rcpp::wrap(snn);
  } catch (const std::exception& ex) {
    Rf_error("%s", ex.what());
  } catch (...) {
    Rf_error("Unknown error in compute_snn_rcpp_fast.");
  }
}

extern "C" SEXP compute_snn_rcpp(SEXP nn_ranked_sexp, double prune) {
  if (!Rf_isMatrix(nn_ranked_sexp) || TYPEOF(nn_ranked_sexp) != REALSXP) {
    Rf_error("nn_ranked must be a numeric matrix.");
  }
  SEXP dim_sexp = Rf_getAttrib(nn_ranked_sexp, R_DimSymbol);
  if (dim_sexp == R_NilValue || LENGTH(dim_sexp) != 2) {
    Rf_error("nn_ranked must be a matrix with Dim attribute.");
  }
  const int* dims = INTEGER(dim_sexp);
  const Eigen::Map<const Eigen::MatrixXd> nn_map(
      REAL(nn_ranked_sexp),
      dims[0],
      dims[1]);
  try {
    return Rcpp::wrap(compute_snn_from_matrix(nn_map, prune));
  } catch (const std::exception& ex) {
    Rf_error("%s", ex.what());
  } catch (...) {
    Rf_error("Unknown error in compute_snn_rcpp.");
  }
}
#endif

static int export_snn_csc(
    const Eigen::SparseMatrix<double>& snn,
    double** out_x,
    int** out_i,
    int** out_p,
    int* out_nnz) {
  const int n = static_cast<int>(snn.cols());
  const int nnz = static_cast<int>(snn.nonZeros());

  double* x_out =
      static_cast<double*>(std::malloc(sizeof(double) * static_cast<size_t>(nnz)));
  int* i_out =
      static_cast<int*>(std::malloc(sizeof(int) * static_cast<size_t>(nnz)));
  int* p_out =
      static_cast<int*>(std::malloc(sizeof(int) * static_cast<size_t>(n + 1)));
  if (x_out == nullptr || i_out == nullptr || p_out == nullptr) {
    std::free(x_out);
    std::free(i_out);
    std::free(p_out);
    return -2;
  }

  int nz = 0;
  for (int col = 0; col < n; ++col) {
    p_out[col] = nz;
    for (Eigen::SparseMatrix<double>::InnerIterator it(snn, col); it; ++it) {
      i_out[nz] = static_cast<int>(it.row());
      x_out[nz] = it.value();
      ++nz;
    }
  }
  p_out[n] = nz;

  *out_x = x_out;
  *out_i = i_out;
  *out_p = p_out;
  *out_nnz = nnz;
  return 0;
}

static int export_snn_csc_into(
    const Eigen::SparseMatrix<double>& snn,
    double* out_x,
    int out_x_len,
    int* out_i,
    int out_i_len,
    int* out_p) {
  const int n = static_cast<int>(snn.cols());
  const int nnz = static_cast<int>(snn.nonZeros());
  if (out_x_len < nnz || out_i_len < nnz) {
    return -2;
  }

  int nz = 0;
  for (int col = 0; col < n; ++col) {
    out_p[col] = nz;
    for (Eigen::SparseMatrix<double>::InnerIterator it(snn, col); it; ++it) {
      out_i[nz] = static_cast<int>(it.row());
      out_x[nz] = it.value();
      ++nz;
    }
  }
  out_p[n] = nz;
  return nz;
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
    Eigen::SparseMatrix<double> snn = build_snn_matrix(nn_ranked, nrows, ncols);
    scale_prune_snn(snn, ncols, prune);
    return export_snn_csc(snn, out_x, out_i, out_p, out_nnz);
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
  if (nn_ranked == nullptr || nrows < 1 || ncols < 1) {
    set_error(error_msg, error_msg_len, "Invalid nn_ranked input.");
    return -1;
  }

  try {
    static thread_local Eigen::SparseMatrix<double> tls_snn;
    static thread_local bool tls_ready = false;

    if (!tls_ready) {
      tls_snn = build_snn_matrix(nn_ranked, nrows, ncols);
      scale_prune_snn(tls_snn, ncols, prune);
      tls_ready = true;
    }

    const bool sizing_only = out_x == nullptr || out_i == nullptr;
    if (sizing_only) {
      if (out_nnz_required != nullptr) {
        *out_nnz_required = static_cast<int>(tls_snn.nonZeros());
      }
      return -3;
    }

    if (out_p == nullptr || out_p_len < nrows + 1) {
      set_error(error_msg, error_msg_len, "Null output pointer.");
      tls_ready = false;
      tls_snn.resize(0, 0);
      return -1;
    }

    const int rc = export_snn_csc_into(
        tls_snn, out_x, out_x_len, out_i, out_i_len, out_p);
    tls_ready = false;
    tls_snn.resize(0, 0);
    return rc;
  } catch (const std::exception& ex) {
    set_error(error_msg, error_msg_len, ex.what());
    return -1;
  } catch (...) {
    set_error(error_msg, error_msg_len, "Unknown error in compute_snn_csc_into.");
    return -1;
  }
}

extern "C" void compute_snn_csc_clear_cache() {}
