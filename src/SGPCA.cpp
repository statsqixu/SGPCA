#include <RcppArmadillo.h>
using namespace Rcpp;
#include <iostream>
#include <string>

double soft_thresholding_elementwise(double x, double lambda) {
  return std::copysign(std::max(0.0, std::abs(x) - lambda), x);
}

void printProgressBar(int current, int total) {
  int barWidth = 70;
  float progress = (float)current / total;
  Rcpp::Rcout << "[";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos) Rcpp::Rcout << "=";
    else if (i == pos) Rcpp::Rcout << ">";
    else Rcpp::Rcout << " ";
  }
  Rcpp::Rcout << "] " << int(progress * 100.0) << " %\r";
  Rcpp::Rcout.flush();
}

// [[Rcpp::export]]
arma::vec SGPCA_cpp(const arma::mat& X, const arma::ivec& group_label,
                       double tau_reg, double eta_reg, int max_iter, double tol,
                       const arma::vec& init_v) {

  int G = arma::max(group_label);

  // Initialize v from sparse PCA
  arma::vec v = init_v;
  arma::vec v_old = v;
  arma::vec gamma = X.t() * X * v / X.n_rows;
  arma::vec gamma_old = gamma;

  for (int iter = 0; iter < max_iter; ++iter) {

    v_old = v;
    gamma_old = gamma;
    gamma = X * v;
    gamma = X.t() * gamma / X.n_rows;

    for (int g = 1; g <= G; ++g) {
      arma::uvec index_g = arma::find(group_label == g);
      arma::vec v_g = v.elem(index_g);
      arma::vec gamma_g = gamma.elem(index_g);

      // Calculate the scaling factor based on group size
      double group_size = index_g.n_elem;
      double scaled_penalty = eta_reg * std::sqrt(group_size);

      // Step 1: Group-level thresholding
      double gamma_g_norm = arma::norm(gamma_g, 2);
      if (gamma_g_norm > 0) {

        double group_shrinkage = 1 - scaled_penalty / gamma_g_norm;

        if (group_shrinkage > 0) {
          // Group survives, apply soft thresholding
          gamma_g = group_shrinkage * gamma_g;

          // Step 2: Element-wise soft thresholding for non-zero groups
          for (size_t i = 0; i < v_g.n_elem; ++i) {
            v_g[i] = soft_thresholding_elementwise(gamma_g[i], tau_reg);
          }
        } else {
          gamma_g.zeros();
        }
      }

      gamma.elem(index_g) = gamma_g;
    }

    // Normalization
    if (arma::norm(gamma, 2) > 0) {
      v = gamma / arma::norm(gamma, 2);
    } else {
      v.zeros();
      break;
    }

    // Stopping rule
    if (arma::norm(v - v_old, 2) < tol) {
      break;
    }
  }

  return v;
}

// [[Rcpp::export]]
Rcpp::List SGPCA_rs_cpp(const arma::mat& X, const arma::ivec& group_label, int B, double rho,
                        const arma::vec& tau_range, const arma::vec& eta_range, int max_iter, double tol, std::string mode,
                        const arma::vec& init_v){

  int N = X.n_rows;
  int P = X.n_cols;

  arma::cube rs_results(tau_range.n_elem * eta_range.n_elem, B, P, arma::fill::zeros);

  for (size_t param_index = 0; param_index < tau_range.n_elem * eta_range.n_elem; ++param_index) {
    double tau = tau_range(param_index % tau_range.n_elem);
    double eta = eta_range(param_index / eta_range.n_elem);

    for (int rs_index = 0; rs_index < B; ++rs_index) {
      arma::uvec resample_indices = arma::randperm(N, floor(rho * N));
      arma::mat X_resample = X.rows(resample_indices);

      // Call to sgPCA_cpp_v3 function
      arma::vec v = SGPCA_cpp(X_resample, group_label, tau, eta, max_iter, tol, init_v);

      if (arma::norm(v, 2) > 0) {
        rs_results.tube(param_index, rs_index) = v / arma::norm(v, 2);
      } else {
        rs_results.tube(param_index, rs_index).zeros();
      }
    }
    // Print progress bar
    printProgressBar(param_index + 1, tau_range.n_elem * eta_range.n_elem);
  }

  Rcpp::Rcout << std::endl;

  arma::vec alignment_results(tau_range.n_elem * eta_range.n_elem, arma::fill::zeros);
  arma::vec mean_support_sizes(tau_range.n_elem * eta_range.n_elem, arma::fill::zeros);

  for (size_t param_index = 0; param_index < tau_range.n_elem * eta_range.n_elem; ++param_index) {
    // Calculate alignment
    arma::mat loadings_rs = rs_results.tube(param_index, 0, param_index, B - 1);
    arma::mat inner_product_matrix = arma::abs(loadings_rs * loadings_rs.t());
    alignment_results(param_index) = (arma::accu(inner_product_matrix) - B) / (B * (B - 1));

    // Calculate mean support size for features
    double total_support_size = 0.0;
    for (int rs_index = 0; rs_index < B; ++rs_index) {
      arma::uvec support_indices = arma::find(arma::abs(loadings_rs.row(rs_index)) > 1e-10);
      total_support_size += support_indices.n_elem;
    }
    mean_support_sizes(param_index) = total_support_size / B;
  }

  return Rcpp::List::create(
    Rcpp::Named("rs_results") = rs_results,
    Rcpp::Named("alignment") = alignment_results,
    Rcpp::Named("mean_support_size") = mean_support_sizes
  );

}
