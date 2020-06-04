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

  // mat Metric;
  cube Embedding;
  // cube integrated_all;
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
  // uvec dataset_list;
  double quantile_default;
  int rept_ds;
  // int k;
  // int max_k;

  // FIRM_res(mat& Metric, const mat& SS2, const uvec& hvg_ind_SS2, const mat& SS2_FindClusters  ,
  //          const mat& tenx, const uvec& hvg_ind_tenx, const mat& tenx_FindClusters,
  //          const int dims, const int gene_all_num, const uvec gene_all_hvg_ind,
  //          const uvec gene_all_ind_SS2, const uvec gene_all_ind_tenx,
  //          const uvec dataset_list, double quantile_default, const int rept_ds = 50,
  //          const int k = 5, const int max_k = 300){
  FIRM_res(cube& Embedding, const mat& integrated_PCA_embedding, const mat& SS2, const uvec& hvg_ind_SS2, const mat& SS2_FindClusters  ,
           const mat& tenx, const uvec& hvg_ind_tenx, const mat& tenx_FindClusters,
           const int dims, const int gene_all_num, const uvec gene_all_hvg_ind,
           const uvec gene_all_ind_SS2, const uvec gene_all_ind_tenx,
           double quantile_default, const int rept_ds = 50){
  // FIRM_res(cube& integrated_all, const mat& SS2, const uvec& hvg_ind_SS2, const mat& SS2_FindClusters  ,
  //          const mat& tenx, const uvec& hvg_ind_tenx, const mat& tenx_FindClusters,
  //          const int dims, const int gene_all_num, const uvec gene_all_hvg_ind,
  //          const uvec gene_all_ind_SS2, const uvec gene_all_ind_tenx,
  //          const uvec dataset_list, double quantile_default, const int rept_ds = 50){

    // this -> Metric = Metric;
    this -> Embedding = Embedding;
    // this -> integrated_all = integrated_all;
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
    // this -> dataset_list = dataset_list;
    this -> quantile_default = quantile_default;
    this -> rept_ds = rept_ds;
    // this -> k = k;
    // this -> max_k = max_k;

  }

  // void loop_by_thread(int i, int j);
  void loop_by_thread(int i, int j, int k);
  void fit_by_thread(int thread_id);
  int next();
};

#endif
