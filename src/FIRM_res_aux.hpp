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
  mat SS2_FindClusters;
  mat tenx;
  mat tenx_FindClusters;
  int dims;
  double quantile_default;
  int rept_ds;

  FIRM_res(cube& Embedding, const mat& integrated_PCA_embedding, const mat& SS2, const mat& SS2_FindClusters  ,
           const mat& tenx, const mat& tenx_FindClusters,
           const int dims, double quantile_default, const int rept_ds = 50){

    this -> Embedding = Embedding;
    this -> integrated_PCA_embedding = integrated_PCA_embedding;
    this -> SS2 = SS2;
    this -> SS2_FindClusters = SS2_FindClusters;
    this -> tenx = tenx;
    this -> tenx_FindClusters = tenx_FindClusters;
    this -> dims = dims;
    this -> quantile_default = quantile_default;
    this -> rept_ds = rept_ds;

  }

  void loop_by_thread(int i, int j, int k);
  void fit_by_thread(int thread_id);
  int next();
};

#endif
