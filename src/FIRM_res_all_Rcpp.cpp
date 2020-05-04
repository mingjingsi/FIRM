// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]
#include <RcppArmadillo.h>
#include "FIRM_res_aux.hpp"
// #include "Update_alpha.hpp"
// #include "Update_beta.hpp"
// #include "functions.hpp"
#include "MixingMetric.hpp"
#include <stdio.h>
#include <math.h>
// using namespace Rcpp;
using namespace arma;
// using namespace std;


// [[Rcpp::export]]
RcppExport SEXP FIRM_res_all(arma::mat& SS2, arma::uvec& hvg_ind_SS2, arma::mat& SS2_FindClusters,
                             arma::mat& tenx, arma::uvec& hvg_ind_tenx, arma::mat& tenx_FindClusters,
                             const int dims, const int gene_all_num, arma::uvec gene_all_hvg_ind,
                             arma::uvec gene_all_ind_SS2, arma::uvec gene_all_ind_tenx,
                             arma::uvec dataset_list, double quantile_default, const int rept_ds = 50,
                             const int k = 5, const int max_k = 300, const int coreNum = 1){

  uword nSS2 = SS2.n_cols;    // No. of cells in SS2
  uword ntenx = tenx.n_cols;  // No. of cells in 10X

  mat integrated_PCA = zeros<mat>(gene_all_num, nSS2 + ntenx);
  mat integrated_PCA_SS2 = zeros<mat>(gene_all_num, nSS2);
  integrated_PCA_SS2.rows(gene_all_ind_SS2-1) = SS2;
  mat integrated_PCA_tenx = zeros<mat>(gene_all_num, ntenx);
  integrated_PCA_tenx.rows(gene_all_ind_tenx-1) = tenx;
  integrated_PCA.cols(0, nSS2-1) = integrated_PCA_SS2;
  integrated_PCA.cols(nSS2, nSS2+ntenx-1) = integrated_PCA_tenx;

  mat integrated_PCA_scaled = integrated_PCA.each_col()/stddev(integrated_PCA, 0, 1).as_col();
  integrated_PCA_scaled(find(integrated_PCA_scaled > 10)).fill(10);

  mat integrated_PCA_hvg = integrated_PCA_scaled.rows(gene_all_hvg_ind-1);
  mat U_integrated_PCA;
  vec s_integrated_PCA;
  mat V_integrated_PCA;
  svd(U_integrated_PCA, s_integrated_PCA, V_integrated_PCA, integrated_PCA_hvg.t());
  mat integrated_PCA_embedding = U_integrated_PCA.cols(0, dims-1) * diagmat(s_integrated_PCA.subvec(0, dims-1));

  double Metric_PCA = mean(MixingMetric(integrated_PCA_embedding, dataset_list, k, max_k));

  mat Metric = ones<mat>(SS2_FindClusters.n_cols, tenx_FindClusters.n_cols)*Metric_PCA;

  FIRM_res FIRM_resObj (Metric, SS2, hvg_ind_SS2, SS2_FindClusters, tenx, hvg_ind_tenx, tenx_FindClusters,
                        dims, gene_all_num, gene_all_hvg_ind, gene_all_ind_SS2, gene_all_ind_tenx,
                        dataset_list, quantile_default, rept_ds, k, max_k);

  const int n_thread = coreNum;
  std::vector<std::thread> threads(n_thread);
  for (int i_thread = 0; i_thread < n_thread; i_thread++){
    threads[i_thread] = std::thread(&FIRM_res::fit_by_thread, &FIRM_resObj, i_thread);
  }

  for (int i = 0; i < n_thread; i++){
    threads[i].join();
  }

  Rcpp::List ret;
  ret["Metric_FIRM"] = FIRM_resObj.Metric;
  ret["Metric_PCA"] = Metric_PCA;

  return ret;
}

