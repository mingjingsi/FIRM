#include "FIRM_res_aux.hpp"
#include "functions.hpp"
#include <stdio.h>
#include <math.h>
using namespace arma;

void FIRM_res::loop_by_thread(int idx, int res_SS2_ind, int res_tenx_ind){

  uword nSS2 = SS2.n_cols;    // No. of cells in SS2
  uword ntenx = tenx.n_cols;  // No. of cells in 10X

  mat SS2_scale = SS2.rows(hvg_ind_SS2 - 1);
  mat tenx_scale = tenx.rows(hvg_ind_tenx - 1);

  vec SS2_FindClusters_current = SS2_FindClusters.col(res_SS2_ind);
  vec tenx_FindClusters_current = tenx_FindClusters.col(res_tenx_ind);

  uword nSS2_cluster = SS2_FindClusters_current.max() + 1;   // No. of clusters in SS2
  uword ntenx_cluster = tenx_FindClusters_current.max() + 1; // No. of clusters in 10X

  uvec SS2_cluster_num_ini = hist(SS2_FindClusters_current, nSS2_cluster);
  uvec tenx_cluster_num_ini = hist(tenx_FindClusters_current, ntenx_cluster);

  uvec SS2_cluster_name_ini = zeros<uvec>(nSS2_cluster);
  for (int i = 0; i < nSS2_cluster; i++){
    SS2_cluster_name_ini(i) = i;
  }

  ///// find same clusters
  // calculate the center of each cluster based on cell embeddings of PCA
  mat all_ini_SS2 = zeros<mat>(nSS2, dims);
  mat all_ini_tenx = zeros<mat>(ntenx, dims);
  join_scale_PCA_sep(all_ini_SS2, all_ini_tenx, SS2_scale, tenx_scale, dims);

  mat SS2_center_ini = zeros<mat>(nSS2_cluster, dims);
  for (int i = 0; i < nSS2_cluster; i++){
    SS2_center_ini.row(i) = mean(all_ini_SS2.rows(find(SS2_FindClusters_current == i)), 0);
  }
  mat tenx_center_ini = zeros<mat>(ntenx_cluster, dims);
  for (int i = 0; i < ntenx_cluster; i++){
    tenx_center_ini.row(i) = mean(all_ini_tenx.rows(find(tenx_FindClusters_current == i)), 0);
  }

  uvec remove_ind = ones<uvec>(nSS2_cluster)*999;
  for (int i = 0; i < nSS2_cluster; i++){
    mat SS2_center_ini_i = SS2_center_ini;
    SS2_center_ini_i.shed_row(i);
    uword ind = (sum(square(SS2_center_ini_i.each_row() - SS2_center_ini.row(i)), 1)).index_min();
    if (ind >= i){
      ind = ind + 1;
    }
    mat all_ini_i = all_ini_SS2.rows(find(SS2_FindClusters_current == i));
    vec dist_ini_i = (sum(square(all_ini_i.each_row() - SS2_center_ini.row(i)), 1)).as_col();
    mat all_ini_ind = all_ini_SS2.rows(find(SS2_FindClusters_current == ind));
    vec dist_ini_ind = (sum(square(all_ini_ind.each_row() - SS2_center_ini.row(i)), 1)).as_col();

    if (mean(dist_ini_i) > 1.8*mean(dist_ini_ind)){
      remove_ind(i) = i;
    }
  }

  uvec SS2_cluster_num;
  uvec SS2_cluster_name;
  mat SS2_center;
  if (min(remove_ind) != 999){
    SS2_cluster_num = SS2_cluster_num_ini(find(remove_ind == 999));
    SS2_cluster_name = SS2_cluster_name_ini(find(remove_ind == 999));
    nSS2_cluster = SS2_cluster_name.n_elem;
    SS2_center = SS2_center_ini.rows(find(remove_ind == 999));
  } else {
    SS2_cluster_num = SS2_cluster_num_ini;
    SS2_cluster_name = SS2_cluster_name_ini;
    SS2_center = SS2_center_ini;
  }

  // calculate quantile of the distance from cells to their center within each cluster
  vec SS2_dist_quantile = zeros<vec>(nSS2_cluster);
  vec quantile_def = ones<vec>(1)*quantile_default;
  for (int i = 0; i < nSS2_cluster; i++){
    mat all_SS2_i = all_ini_SS2.rows(find(SS2_FindClusters_current == SS2_cluster_name(i)));
    SS2_dist_quantile(i) = as_scalar(quantile(sum(square(all_SS2_i.each_row() - SS2_center.row(i)), 1).as_col(), quantile_def));
  }

  // find which cluster in SS2 can be merge with each 10X cluster
  // (whether the distance from the center of the j-th 10X clusters to the center of i-th SS2 cluster is less than the quantile for the i-th SS2 cluster)
  // (if two SS2 clusters are detected for one 10X cluster, choose the closer one)
  uvec merge_pair = ones<uvec>(ntenx_cluster)*999;
  uword n_NN = 5;
  uword n_NN_min = std::min(n_NN, nSS2_cluster);
  umat NNC_tenx = zeros<umat>(ntenx_cluster, n_NN_min);
  for (int j = 0; j < ntenx_cluster; j++){
    vec distance = sum(square(SS2_center.each_row() - tenx_center_ini.row(j)), 1).as_col();
    uvec NNC_tenx_all = sort_index(distance);
    uvec NNC_tenx_all_name = SS2_cluster_name(NNC_tenx_all);

    for (int i = 0; i < n_NN_min; i++){
      NNC_tenx(j, i) = NNC_tenx_all_name(i);
    }

    if (sum(distance < SS2_dist_quantile) !=  0){
      for (int i = 0; i < n_NN_min; i++){
        if (distance(NNC_tenx_all(i)) < SS2_dist_quantile(NNC_tenx_all(i))){
          merge_pair(j) = NNC_tenx(j, i);
          break;
        }
      }
    }
  }

  // unpaired cluster id
  uvec tenx_unpaired_ind = find(merge_pair == 999);

  // check the unpaired clusters by subsampling
  if (tenx_unpaired_ind.n_elem != 0){ // if there is unpaired cluster for either SS2 or 10X
    for (int j = 0; j < tenx_unpaired_ind.n_elem; j++){
      uvec prop_tenx_tmp = zeros<uvec>(1);
      prop_tenx_tmp(0) = sum(tenx_FindClusters_current == tenx_unpaired_ind(j));
      double prop_tenx = as_scalar(conv_to<vec>::from(prop_tenx_tmp))/ntenx; //the j-th cluster proportion in 10X

      for (int i = 0; i < NNC_tenx.n_cols; i++){
        uvec prop_SS2_tmp = zeros<uvec>(1);
        prop_SS2_tmp(0) = sum(SS2_FindClusters_current == NNC_tenx(tenx_unpaired_ind(j), i));
        double prop_SS2 = as_scalar(conv_to<vec>::from(prop_SS2_tmp))/nSS2; //the i-th cluster proportion in SS2

        if (prop_tenx < prop_SS2){ // If the proportion in 10X is less than that in SS2, subsample cells in the i-th cluster in SS2
          float num_tmp = sum(SS2_FindClusters_current != NNC_tenx(tenx_unpaired_ind(j), i))*prop_tenx/(1 - prop_tenx);
          fvec num_float = zeros<fvec>(1);
          num_float(0) = std::round(num_tmp);
          uword num = as_scalar(conv_to<uvec>::from(num_float));

          if (num > 20){ // subsample when the No. of cells is not very small
            // calculate the standard deviation based on the cells after subsampling and scale the data
            mat SS2_cluster_tmp = SS2_scale.cols(find(SS2_FindClusters_current == NNC_tenx(tenx_unpaired_ind(j), i)));
            float rept_tmp =  std::round(SS2_cluster_tmp.n_cols/num) + 1;
            float rept_min = 50;
            double rept = std::min(rept_tmp, rept_min);
            vec sd_SS2_tmp = zeros<vec>(SS2_cluster_tmp.n_rows);
            for (int iter = 0; iter < rept; iter++){
              mat SS2_cluster_tmp_new = SS2_cluster_tmp.cols(randperm(SS2_cluster_tmp.n_cols, num));
              mat SS2_all_tmp_new = join_rows(SS2_scale.cols(find(SS2_FindClusters_current != NNC_tenx(tenx_unpaired_ind(j), i))), SS2_cluster_tmp_new);
              sd_SS2_tmp += stddev(SS2_all_tmp_new, 0, 1);
            }
            vec sd_SS2_tmp_new = sd_SS2_tmp/rept;

            int check = check_merge(SS2_scale, sd_SS2_tmp_new, tenx_scale, dims, SS2_FindClusters_current, tenx_FindClusters_current,
                                    NNC_tenx(tenx_unpaired_ind(j), i), tenx_unpaired_ind(j), quantile_def);

            if (check >= 1){
              merge_pair(tenx_unpaired_ind(j)) = NNC_tenx(tenx_unpaired_ind(j), i);
            }

          } else {
            break;
          }
        } else { //If the proportion in 10X is greater than that in SS2, subsample cells in the j-th cluster in 10X
          float num_tmp = sum(tenx_FindClusters_current != tenx_unpaired_ind(j))*prop_SS2/(1 - prop_SS2); //No. of cells in the j-th cluster in 10X after downsampling
          fvec num_float = zeros<fvec>(1);
          num_float(0) = std::round(num_tmp);
          uword num = as_scalar(conv_to<uvec>::from(num_float));

          if (num > 20){ // subsample when the No. of cells is not very small
            // calculate the standard deviation based on the cells after subsampling and scale the data
            mat tenx_cluster_tmp = tenx_scale.cols(find(tenx_FindClusters_current == tenx_unpaired_ind(j)));
            float rept_tmp =  std::round(tenx_cluster_tmp.n_cols/num) + 1;
            float rept_min = 50;
            double rept = std::min(rept_tmp, rept_min);
            vec sd_tenx_tmp = zeros<vec>(tenx_cluster_tmp.n_rows);
            for (int iter = 0; iter < rept; iter++){
              mat tenx_cluster_tmp_new = tenx_cluster_tmp.cols(randperm(tenx_cluster_tmp.n_cols, num));
              mat tenx_all_tmp_new = join_rows(tenx_scale.cols(find(tenx_FindClusters_current != tenx_unpaired_ind(j))), tenx_cluster_tmp_new);
              sd_tenx_tmp += stddev(tenx_all_tmp_new, 0, 1);
            }
            vec sd_tenx_tmp_new = sd_tenx_tmp/rept;

            int check = check_merge(SS2_scale, tenx_scale, sd_tenx_tmp_new, dims, SS2_FindClusters_current, tenx_FindClusters_current,
                                    NNC_tenx(tenx_unpaired_ind(j), i), tenx_unpaired_ind(j), quantile_def);

            if (check >= 1){
              merge_pair(tenx_unpaired_ind(j)) = NNC_tenx(tenx_unpaired_ind(j), i);
            }

          } else {
            break;
          }
        }
      }
    }
  }

  if (arma::all(merge_pair == 999)){

    Embedding.slice(idx) = integrated_PCA_embedding;

  } else {
    uvec SS2_paired_name_ini = unique(merge_pair(find(merge_pair != 999))); // SS2 paired name
    uword n_paired = SS2_paired_name_ini.n_elem; // No. of paired clusters
    uvec tenx_paired_name_ini = zeros<uvec>(n_paired); // corresponding paired id of 10X
    for (int i = 0; i < n_paired; i++){
      uword check = sum(merge_pair == SS2_paired_name_ini(i));
      if (check == 1){
        tenx_paired_name_ini(i) = as_scalar(find(merge_pair == SS2_paired_name_ini(i)));
      } else {
        uvec check_tenx_ind = find(merge_pair == SS2_paired_name_ini(i));
        for (int j = 1; j < check_tenx_ind.n_elem; j++){
          tenx_FindClusters_current(find(tenx_FindClusters_current == check_tenx_ind(j))).fill(check_tenx_ind(0));
        }
        tenx_paired_name_ini(i) = check_tenx_ind(0);
      }
    }

    if (n_paired == 1){
      mat SS2_cluster_tmp = SS2.cols(find(SS2_FindClusters_current == SS2_paired_name_ini(0)));
      vec sd_SS2 = stddev(SS2_cluster_tmp, 0, 1);

      mat tenx_cluster_tmp = tenx.cols(find(tenx_FindClusters_current == tenx_paired_name_ini(0)));
      vec sd_tenx = stddev(tenx_cluster_tmp, 0, 1);

      Embedding.slice(idx) = integrated_scale_fill_hvg_PCA(SS2, tenx, sd_SS2, sd_tenx, gene_all_num, nSS2, ntenx, gene_all_ind_SS2, gene_all_ind_tenx,
                      gene_all_hvg_ind, dims);


    } else{
      ///// compute number of cells in each paired cluster after subsampling
      // the proportion of the corresponding clusters should be the same
      uvec num_SS2_ini = zeros<uvec>(n_paired);
      uvec num_tenx_ini = zeros<uvec>(n_paired);
      for (int i = 0; i < n_paired; i++){
        num_SS2_ini(i) = sum(SS2_FindClusters_current == SS2_paired_name_ini(i));
        num_tenx_ini(i) = sum(tenx_FindClusters_current == tenx_paired_name_ini(i));
      }

      uvec SS2_paired_name;
      uvec tenx_paired_name;
      uvec num_SS2;
      uvec num_tenx;
      pair_info(SS2_paired_name, tenx_paired_name, num_SS2, num_tenx, num_SS2_ini, num_tenx_ini, SS2_paired_name_ini, tenx_paired_name_ini);

      uword n_paired_new = SS2_paired_name.n_elem;

      ///// calculate standard deviation for scaling based on subsampling
      // arma_rng::set_seed(0);
      vec sd_SS2_sum = zeros<vec>(SS2.n_rows);
      vec sd_tenx_sum = zeros<vec>(tenx.n_rows);
      for (int iter = 0; iter < rept_ds; iter++){
        mat SS2_matrix = zeros<mat>(SS2.n_rows, sum(num_SS2));
        mat tenx_matrix = zeros<mat>(tenx.n_rows, sum(num_tenx));
        uword ind_begin_SS2 = 0;
        uword ind_begin_tenx = 0;
        for (int i = 0; i < n_paired_new; i++){
          mat SS2_cluster = SS2.cols(find(SS2_FindClusters_current == SS2_paired_name(i)));
          mat tenx_cluster = tenx.cols(find(tenx_FindClusters_current == tenx_paired_name(i)));

          SS2_matrix.cols(ind_begin_SS2, ind_begin_SS2+num_SS2(i)-1) = SS2_cluster.cols(randperm(SS2_cluster.n_cols, num_SS2(i)));
          tenx_matrix.cols(ind_begin_tenx, ind_begin_tenx+num_tenx(i)-1) = tenx_cluster.cols(randperm(tenx_cluster.n_cols, num_tenx(i)));

          ind_begin_SS2 += num_SS2(i);
          ind_begin_tenx += num_tenx(i);
        }
        sd_SS2_sum += stddev(SS2_matrix, 0, 1);
        sd_tenx_sum += stddev(tenx_matrix, 0, 1);
      }

      vec sd_SS2 = sd_SS2_sum/rept_ds;
      vec sd_tenx = sd_tenx_sum/rept_ds;

      ///// do scaling and obtain new expression matrix
      Embedding.slice(idx) = integrated_scale_fill_hvg_PCA(SS2, tenx, sd_SS2, sd_tenx, gene_all_num, nSS2, ntenx, gene_all_ind_SS2, gene_all_ind_tenx,
                      gene_all_hvg_ind, dims);

    }
  }
}


std::mutex mtx;
int FIRM_res::next(){
  std::lock_guard<std::mutex> lock(mtx);
  if (current_idx >= SS2_FindClusters.n_cols*tenx_FindClusters.n_cols){
    return -1;
  }
  current_idx++;

  return current_idx-1;
}


void FIRM_res::fit_by_thread(int thread_id){
  while(true){
    int idx = next();
    if (idx == -1){
      break;
    }
    int current_idx_SS2 = idx/tenx_FindClusters.n_cols;
    int current_idx_tenx = idx%tenx_FindClusters.n_cols;
    printf("thread_id = %d, idx = %d, SS2_idx = %d, tenx_idx = %d, SS2_idx_num = %d, tenx_idx_num = %d \n",
           thread_id, idx, current_idx_SS2, current_idx_tenx, SS2_FindClusters.n_cols, tenx_FindClusters.n_cols);
    loop_by_thread(idx, current_idx_SS2, current_idx_tenx);
  }
}
