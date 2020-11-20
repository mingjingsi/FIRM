#ifndef FIRM_res_aux_hpp
#define FIRM_res_aux_hpp

#include <stdio.h>
#include <RcppArmadillo.h>
#include <thread>
#include <mutex>

using namespace std;
using namespace arma;

class FIRM_res{
public:
  uword current_idx = 0;

  cube Embedding;
  mat integrated_PCA_embedding;
  mat SS2;
  uvec hvg_ind_SS2;
  mat SS2_FindClusters;
  mat tenx;
  uvec hvg_ind_tenx;
  mat tenx_FindClusters;
  int dims;
  int gene_all_num;
  uvec gene_all_hvg_ind;
  uvec gene_all_ind_SS2;
  uvec gene_all_ind_tenx;
  double quantile_default;
  int rept_ds;

  FIRM_res(cube& Embedding, const mat& integrated_PCA_embedding, const mat& SS2, const uvec& hvg_ind_SS2, const mat& SS2_FindClusters  ,
           const mat& tenx, const uvec& hvg_ind_tenx, const mat& tenx_FindClusters,
           const int dims, const int gene_all_num, const uvec gene_all_hvg_ind,
           const uvec gene_all_ind_SS2, const uvec gene_all_ind_tenx,
           double quantile_default, const int rept_ds = 50){

    this -> Embedding = Embedding;
    this -> integrated_PCA_embedding = integrated_PCA_embedding;
    this -> SS2 = SS2;
    this -> hvg_ind_SS2 = hvg_ind_SS2;
    this -> SS2_FindClusters = SS2_FindClusters;
    this -> tenx = tenx;
    this -> hvg_ind_tenx = hvg_ind_tenx;
    this -> tenx_FindClusters = tenx_FindClusters;
    this -> dims = dims;
    this -> gene_all_num = gene_all_num;
    this -> gene_all_hvg_ind = gene_all_hvg_ind;
    this -> gene_all_ind_SS2 = gene_all_ind_SS2;
    this -> gene_all_ind_tenx = gene_all_ind_tenx;
    this -> quantile_default = quantile_default;
    this -> rept_ds = rept_ds;

  }

  void loop_by_thread(int i, int j, int k);
  void fit_by_thread(int thread_id);
  int next();
};

#endif
