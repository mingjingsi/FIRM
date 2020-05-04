// #include <RcppArmadillo.h>
// #include <Rcpp.h>
// #include <ANN/ANN.h>
// #include "nn.hpp"
#include "MixingMetric.hpp"
// using namespace arma;

vec MixingMetric(mat embedding, uvec dataset_list, const int k = 5, const int max_k = 300) {

  umat nn = nn2_cpp(embedding, max_k+1);
  nn.shed_col(0);

  uvec dataset = unique(dataset_list);
  umat dataset_rank = ones<umat>(dataset.n_elem, nn.n_rows)*(max_k-1);

  for (int i = 0; i < nn.n_rows; i++){
  // int i = 2;
    int nn_SS2_tmp = 0;
    int nn_tenx_tmp = 0;

    for (int k_tmp = 0; k_tmp < max_k; k_tmp++){
      if(dataset_list(nn(i, k_tmp)) == dataset(0)){
        // printf("if1 %d ", k_tmp);
        if (nn_SS2_tmp == k-1){
          dataset_rank(0, i) = k_tmp;
          // printf("if end %d ", k_tmp);
          break;
        }
        nn_SS2_tmp += 1;
      }
      if(k_tmp == max_k-1){
        dataset_rank(0, i) = max_k-1;
      }
    }

    for (int k_tmp = 0; k_tmp < max_k; k_tmp++){
      if(dataset_list(nn(i, k_tmp)) == dataset(1)){
        if(nn_tenx_tmp == k - 1){
          dataset_rank(1, i) = k_tmp;
          break;
        }
        nn_tenx_tmp += 1;
      }
      if(k_tmp == max_k-1){
        dataset_rank(1, i) = max_k-1;
      }
    }
  }

  mat dataset_rank_double = conv_to<mat>::from(dataset_rank+1);

  vec metric = median(dataset_rank_double, 0).as_col();

  return metric;
}

