// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "FIRM_res_aux.hpp"
#include "functions.hpp"
#include <stdio.h>
#include <math.h>
using namespace arma;

// [[Rcpp::export]]
RcppExport SEXP FIRM_res_all(arma::mat& SS2, arma::mat& SS2_FindClusters,
                             arma::mat& tenx, arma::mat& tenx_FindClusters,
                             const int dims,
                             double quantile_default, const int rept_ds = 50,
                             const int k = 5, const int max_k = 300, const int coreNum = 1){

  uword nSS2 = SS2.n_cols;    // No. of cells in SS2
  uword ntenx = tenx.n_cols;  // No. of cells in 10X

  mat integrated_PCA_embedding = integrated_fill_hvg_PCA(SS2, tenx, dims);

  cube Embedding = zeros<cube>(integrated_PCA_embedding.n_rows, integrated_PCA_embedding.n_cols, SS2_FindClusters.n_cols*tenx_FindClusters.n_cols);

  FIRM_res FIRM_resObj(Embedding, integrated_PCA_embedding, SS2, SS2_FindClusters, tenx, tenx_FindClusters,
                       dims, quantile_default, rept_ds);

  const int n_thread = coreNum;
  std::vector<std::thread> threads(n_thread);
  for (int i_thread = 0; i_thread < n_thread; i_thread++){
    threads[i_thread] = std::thread(&FIRM_res::fit_by_thread, &FIRM_resObj, i_thread);
  }

  for (int i = 0; i < n_thread; i++){
    threads[i].join();
  }

  Rcpp::List ret;

  ret["Embedding_FIRM"] = FIRM_resObj.Embedding;
  ret["Embedding_PCA"] = integrated_PCA_embedding;

  return ret;
}

