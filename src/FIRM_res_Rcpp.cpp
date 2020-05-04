// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>
using namespace arma;

// [[Rcpp::export]]
RcppExport SEXP FIRM_res(arma::mat& SS2, arma::uvec& hvg_ind_SS2, arma::mat& SS2_FindClusters  ,
                         arma::mat& tenx, arma::uvec& hvg_ind_tenx, arma::mat& tenx_FindClusters,
                         const int dims, const int gene_all_num, arma::uvec gene_all_hvg_ind,
                         arma::uvec gene_all_ind_SS2, arma::uvec gene_all_ind_tenx,
                         arma::uvec dataset_list, double quantile_default, const int rept_ds=50,
                         const int k = 5, const int max_k = 300){

  uword nSS2 = SS2.n_cols;    // No. of cells in SS2
  uword ntenx = tenx.n_cols;  // No. of cells in 10X

  mat SS2_scale = SS2.rows(hvg_ind_SS2-1);
  mat tenx_scale = tenx.rows(hvg_ind_tenx-1);

  vec SS2_FindClusters_current = SS2_FindClusters.col(4);
  vec tenx_FindClusters_current = tenx_FindClusters.col(4);

  uword nSS2_cluster = SS2_FindClusters_current.max()+1;   // No. of clusters in SS2
  uword ntenx_cluster = tenx_FindClusters_current.max()+1; // No. of clusters in 10X

  uvec SS2_cluster_num_ini = hist(SS2_FindClusters_current, nSS2_cluster);
  uvec tenx_cluster_num_ini = hist(tenx_FindClusters_current, ntenx_cluster);

  uvec SS2_cluster_name_ini = zeros<uvec>(nSS2_cluster);
  for (int i = 0; i < nSS2_cluster; i++){
    SS2_cluster_name_ini(i) = i;
  }

  ///// find same clusters
  // calculate the center of each cluster based on cell embeddings of PCA
  mat all_matrix_ini = join_rows(SS2_scale, tenx_scale);
  mat all_matrix_ini_scaled = all_matrix_ini.each_col()/stddev(all_matrix_ini, 0, 1).as_col();
  mat U;
  vec s;
  mat V;
  svd(U, s, V, all_matrix_ini_scaled.t());
  mat all_ini = U.cols(0, dims-1) * diagmat(s.subvec(0, dims-1));
  mat all_ini_SS2 = all_ini.rows(0, nSS2-1);
  mat all_ini_tenx = all_ini.rows(nSS2, nSS2+ntenx-1);

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
    uword ind = sum(square(SS2_center_ini_i.each_row() - SS2_center_ini.row(i)), 1).index_min();
    if (ind >= i){
      ind = ind + 1;
    }
    mat all_ini_i = all_ini_SS2.rows(find(SS2_FindClusters_current == i));
    vec dist_ini_i = sum(square(all_ini_i.each_row() - SS2_center_ini.row(i)), 1).as_col();
    mat all_ini_ind = all_ini_SS2.rows(find(SS2_FindClusters_current == ind));
    vec dist_ini_ind = sum(square(all_ini_ind.each_row() - SS2_center_ini.row(i)), 1).as_col();

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
      double prop_tenx = sum(tenx_FindClusters_current == tenx_unpaired_ind(j))/ntenx; //the j-th cluster proportion in 10X

      for (int i = 0; i < NNC_tenx.n_cols; i++){
        double prop_SS2 = sum(SS2_FindClusters_current == NNC_tenx(j, i))/nSS2; //the i-th cluster proportion in SS2

        if (prop_tenx < prop_SS2){ // If the proportion in 10X is less than that in SS2, subsample cells in the i-th cluster in SS2
          float num_tmp = sum(SS2_FindClusters_current != NNC_tenx(j, i))*prop_tenx/(1 - prop_tenx);
          fvec num_float = zeros<fvec>(1);
          num_float(0) = std::round(num_tmp);
          uword num = as_scalar(conv_to<uvec>::from(num_float));

          if (num > 20){ // subsample when the No. of cells is not very small
            // calculate the standard deviation based on the cells after subsampling and scale the data
            mat SS2_cluster_tmp = SS2_scale.cols(find(SS2_FindClusters_current == NNC_tenx(j, i)));
            float rept_tmp =  std::round(SS2_cluster_tmp.n_cols/num) + 1;
            float rept_min = 50;
            double rept = std::min(rept_tmp, rept_min);
            mat sd_SS2_tmp = zeros<mat>(rept, SS2_cluster_tmp.n_rows);
            // arma_rng::set_seed(0);
            for (int iter = 0; iter < rept; iter++){
              mat SS2_cluster_tmp_new = SS2_cluster_tmp.cols(randperm(SS2_cluster_tmp.n_cols, num));
              mat SS2_all_tmp_new = join_rows(SS2_scale.cols(find(SS2_FindClusters_current != NNC_tenx(j, i))), SS2_cluster_tmp_new);
              sd_SS2_tmp.row(iter) = stddev(SS2_all_tmp_new, 0, 1);
            }
            vec sd_SS2_tmp_new = mean(sd_SS2_tmp, 1).as_col();
            mat SS2_scale_tmp_new = SS2_scale.each_col()/(sd_SS2_tmp_new + (sd_SS2_tmp_new == 0));
            mat all_matrix_tmp_new = join_rows(SS2_scale_tmp_new, tenx_scale);

            // calculate the new center of each cluster based on cell embeddings of PCA
            mat all_matrix_tmp_new_scaled = all_matrix_tmp_new.each_col()/stddev(all_matrix_tmp_new, 0, 1).as_col();
            mat U_tmp;
            vec s_tmp;
            mat V_tmp;
            svd(U_tmp, s_tmp, V_tmp, all_matrix_tmp_new_scaled.t());
            mat all_tmp = U_tmp.cols(0, dims-1) * diagmat(s_tmp.subvec(0, dims-1));
            mat all_SS2_tmp = all_tmp.rows(0, SS2_scale_tmp_new.n_cols-1);
            mat all_tenx_tmp = all_tmp.rows(SS2_scale_tmp_new.n_cols, all_tmp.n_cols-1);

            vec SS2_cluster_center = mean(all_SS2_tmp.rows(find(SS2_FindClusters_current == NNC_tenx(j, i))), 0);
            vec tenx_cluster_center = mean(all_tenx_tmp.rows(find(tenx_FindClusters_current == tenx_unpaired_ind(j))), 0);

            mat all_SS2_tmp_i = all_SS2_tmp.rows(find(SS2_FindClusters_current == NNC_tenx(j, i)));
            double SS2_cluster_dist_quantile = as_scalar(quantile(sum(square(all_SS2_tmp_i.each_row() - SS2_cluster_center), 1).as_col(), quantile_def));

            mat all_tenx_tmp_i = all_tenx_tmp.rows(find(tenx_FindClusters_current == tenx_unpaired_ind(j)));
            double tenx_cluster_dist_quantile = as_scalar(quantile(sum(square(all_tenx_tmp_i.each_row() - tenx_cluster_center), 1).as_col(), quantile_def));

            if ((sum(square(SS2_cluster_center - tenx_cluster_center)) < SS2_cluster_dist_quantile) ||
                (sum(square(SS2_cluster_center - tenx_cluster_center)) < tenx_cluster_dist_quantile)){
              merge_pair(tenx_unpaired_ind(j)) = NNC_tenx(j, i);
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
            mat sd_tenx_tmp = zeros<mat>(rept, tenx_cluster_tmp.n_rows);
            // arma_rng::set_seed(0);
            for (int iter = 0; iter < rept; iter++){
              mat tenx_cluster_tmp_new = tenx_cluster_tmp.cols(randperm(tenx_cluster_tmp.n_cols, num));
              mat tenx_all_tmp = join_rows(tenx_scale.cols(find(tenx_FindClusters_current != tenx_unpaired_ind(j))), tenx_cluster_tmp_new);
              sd_tenx_tmp.row(iter) = stddev(tenx_all_tmp, 0, 1);
            }
            vec sd_tenx_tmp_new = mean(sd_tenx_tmp, 1).as_col();
            vec tenx_scale_tmp_new = tenx_scale.each_col()/(sd_tenx_tmp_new + (sd_tenx_tmp_new == 0));
            mat all_matrix_tmp_new = join_rows(SS2_scale, tenx_scale_tmp_new);

            // calculate the new center of each cluster based on cell embeddings of PCA
            mat all_matrix_tmp_new_scaled = all_matrix_tmp_new.each_col()/stddev(all_matrix_tmp_new, 0, 1).as_col();
            mat U_tmp;
            vec s_tmp;
            mat V_tmp;
            svd(U_tmp, s_tmp, V_tmp, all_matrix_tmp_new_scaled.t());
            mat all_tmp = U_tmp.cols(0, dims-1) * diagmat(s_tmp.subvec(0, dims-1));
            mat all_SS2_tmp = all_tmp.rows(0, SS2_scale.n_cols-1);
            mat all_tenx_tmp = all_tmp.rows(SS2_scale.n_cols, all_tmp.n_cols-1);

            vec SS2_cluster_center = mean(all_SS2_tmp.rows(find(SS2_FindClusters_current == NNC_tenx(j, i))), 0);
            vec tenx_cluster_center = mean(all_tenx_tmp.rows(find(tenx_FindClusters_current == tenx_unpaired_ind(j))), 0);

            // vec quantile_def = ones<vec>(1)*quantile_default;
            mat all_SS2_tmp_i = all_SS2_tmp.rows(find(SS2_FindClusters_current == NNC_tenx(j, i)));
            double SS2_cluster_dist_quantile = as_scalar(quantile(sum(square(all_SS2_tmp_i.each_row() - SS2_cluster_center), 1).as_col(), quantile_def));

            mat all_tenx_tmp_i = all_tenx_tmp.rows(find(tenx_FindClusters_current == tenx_unpaired_ind(j)));
            double tenx_cluster_dist_quantile = as_scalar(quantile(sum(square(all_tenx_tmp_i.each_row() - tenx_cluster_center), 1).as_col(), quantile_def));

            if ((sum(square(SS2_cluster_center - tenx_cluster_center)) < SS2_cluster_dist_quantile) ||
                (sum(square(SS2_cluster_center - tenx_cluster_center)) < tenx_cluster_dist_quantile)){
              merge_pair(tenx_unpaired_ind(j)) = NNC_tenx(j, i);
            }

          } else {
            break;
          }
        }
      }
    }
  }

  if (arma::all(merge_pair == 999)){
    Rcpp::List ret;
    ret["SS2"] = SS2;
    ret["tenx"] = tenx;

    return ret;
  }

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
  // if the smallest clusters in SS2 and 10X are paired
  if (num_SS2_ini.index_min() == num_tenx_ini.index_min()){
    uword ind = num_SS2_ini.index_min();
    mat prop_join = zeros<mat>(num_SS2_ini.n_elem, 2);
    prop_join.col(0) = conv_to<vec>::from(num_SS2_ini)/num_SS2_ini(ind);
    prop_join.col(1) = conv_to<vec>::from(num_tenx_ini)/num_tenx_ini(ind);
    vec prop_paired = min(prop_join, 1).as_col();

    // If more than 75% cells are removed, we don't consider this pair in subsampling
    if((sum(round(prop_paired*num_SS2_ini(ind)))/sum(num_SS2_ini) < 0.25) ||
       (sum(round(prop_paired*num_tenx_ini(ind)))/sum(num_tenx_ini) < 0.25)){
      umat SS2_paired_name_mat = reshape(SS2_paired_name_ini, 1, SS2_paired_name_ini.n_elem);
      SS2_paired_name_mat.shed_col(ind);
      SS2_paired_name = SS2_paired_name_mat.as_col();

      umat tenx_paired_name_mat = reshape(tenx_paired_name_ini, 1, tenx_paired_name_ini.n_elem);
      tenx_paired_name_mat.shed_col(ind);
      tenx_paired_name = tenx_paired_name_mat.as_col();

      umat num_SS2_mat = reshape(num_SS2_ini, 1, num_SS2_ini.n_elem);
      num_SS2_mat.shed_col(ind);
      num_SS2 = num_SS2_mat.as_col();

      umat num_tenx_mat = reshape(num_tenx, 1, num_tenx.n_elem);
      num_tenx_mat.shed_col(ind);
      num_tenx = num_tenx_mat.as_col();
    } else {
      num_SS2 = conv_to<uvec>::from(round(prop_paired*num_SS2_ini(ind)));
      num_tenx = conv_to<uvec>::from(round(prop_paired*num_tenx_ini(ind)));
      SS2_paired_name = SS2_paired_name_ini;
      tenx_paired_name = tenx_paired_name_ini;
    }
  } else {
    // if the smallest clusters in SS2 and 10X are not paired
    //choose which one should be the smallest cluster after subsampling based on the No. of cells kept
    uword ind_SS2 = num_SS2_ini.index_min();
    uword ind_tenx = num_tenx_ini.index_min();
    uvec num_tenx1_tmp = num_tenx_ini;
    num_tenx1_tmp(ind_SS2) = num_tenx_ini(ind_tenx)*num_SS2_ini(ind_SS2)/num_SS2_ini(ind_tenx);
    mat prop_join1 = zeros<mat>(num_SS2_ini.n_elem, 2);
    prop_join1.col(0) = conv_to<vec>::from(num_SS2_ini)/num_SS2_ini(ind_SS2);
    prop_join1.col(1) = conv_to<vec>::from(num_tenx1_tmp)/num_tenx1_tmp(ind_SS2);
    vec prop_paired1 = min(prop_join1, 1).as_col();
    uvec num_SS21 = conv_to<uvec>::from(round(prop_paired1*num_SS2_ini(ind_SS2)));
    uvec num_tenx1 = conv_to<uvec>::from(round(prop_paired1*num_tenx1_tmp(ind_SS2)));

    uvec num_SS22_tmp = num_SS2_ini;
    num_SS22_tmp(ind_tenx) = num_SS2_ini(ind_SS2)*num_tenx_ini(ind_tenx)/num_tenx_ini(ind_SS2);
    mat prop_join2 = zeros<mat>(num_tenx_ini.n_elem, 2);
    prop_join2.col(0) = conv_to<vec>::from(num_SS22_tmp)/num_SS22_tmp(ind_tenx);
    prop_join2.col(1) = conv_to<vec>::from(num_tenx_ini)/num_tenx_ini(ind_tenx);
    vec prop_paired2 = min(prop_join2, 1).as_col();
    uvec num_SS22 = conv_to<uvec>::from(round(prop_paired2*num_SS22_tmp(ind_tenx)));
    uvec num_tenx2 = conv_to<uvec>::from(round(prop_paired2*num_tenx_ini(ind_tenx)));

    if (sum(num_SS21 + num_tenx1) > sum(num_SS22 + num_tenx2)){
      if ((sum(num_SS21)/sum(num_SS2) < 0.5) || (sum(num_tenx1)/sum(num_tenx) < 0.25)){
        umat SS2_paired_name_mat = reshape(SS2_paired_name_ini, 1, SS2_paired_name_ini.n_elem);
        SS2_paired_name_mat.shed_col(ind_SS2);
        SS2_paired_name = SS2_paired_name_mat.as_col();

        umat tenx_paired_name_mat = reshape(tenx_paired_name_ini, 1, tenx_paired_name_ini.n_elem);
        tenx_paired_name_mat.shed_col(ind_SS2);
        tenx_paired_name = tenx_paired_name_mat.as_col();

        umat num_SS2_mat = reshape(num_SS2_ini, 1, num_SS2_ini.n_elem);
        num_SS2_mat.shed_col(ind_SS2);
        num_SS2 = num_SS2_mat.as_col();

        umat num_tenx_mat = reshape(num_tenx_ini, 1, num_tenx_ini.n_elem);
        num_tenx_mat.shed_col(ind_SS2);
        num_tenx = num_tenx_mat.as_col();
      } else{
        num_SS2 = num_SS21;
        num_tenx = num_tenx1;
        SS2_paired_name = SS2_paired_name_ini;
        tenx_paired_name = tenx_paired_name_ini;
      }
    } else{
      if ((sum(num_SS22)/sum(num_SS2) < 0.5) || (sum(num_tenx2)/sum(num_tenx) < 0.25)){
        umat SS2_paired_name_mat = reshape(SS2_paired_name_ini, 1, SS2_paired_name_ini.n_elem);
        SS2_paired_name_mat.shed_col(ind_tenx);
        SS2_paired_name = SS2_paired_name_mat.as_col();

        umat tenx_paired_name_mat = reshape(tenx_paired_name_ini, 1, tenx_paired_name_ini.n_elem);
        tenx_paired_name_mat.shed_col(ind_tenx);
        tenx_paired_name = tenx_paired_name_mat.as_col();

        umat num_SS2_mat = reshape(num_SS2_ini, 1, num_SS2_ini.n_elem);
        num_SS2_mat.shed_col(ind_tenx);
        num_SS2 = num_SS2_mat.as_col();

        umat num_tenx_mat = reshape(num_tenx_ini, 1, num_tenx_ini.n_elem);
        num_tenx_mat.shed_col(ind_tenx);
        num_tenx = num_tenx_mat.as_col();
      } else{
        num_SS2 = num_SS22;
        num_tenx = num_tenx2;
        SS2_paired_name = SS2_paired_name_ini;
        tenx_paired_name = tenx_paired_name_ini;
      }
    }
  }

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
  mat SS2_scale_new = SS2.each_col()/(sd_SS2 + (sd_SS2 == 0));
  mat tenx_scale_new = tenx.each_col()/(sd_tenx + (sd_tenx == 0));

  ///// compute the Mixing Metric
  mat integrated = zeros<mat>(gene_all_num, nSS2 + ntenx);
  mat integrated_SS2 = zeros<mat>(gene_all_num, nSS2);
  integrated_SS2.rows(gene_all_ind_SS2-1) = SS2_scale_new;
  mat integrated_tenx = zeros<mat>(gene_all_num, ntenx);
  integrated_tenx.rows(gene_all_ind_tenx-1) = tenx_scale_new;
  integrated.cols(0, nSS2-1) = integrated_SS2;
  integrated.cols(nSS2, nSS2+ntenx-1) = integrated_tenx;

  mat integrated_scaled = integrated.each_col()/stddev(integrated, 0, 1).as_col();
  integrated_scaled(find(integrated_scaled > 10)).fill(10);

  Rcpp::List ret;

  ret["SS2"] = SS2_scale_new;
  ret["tenx"] = tenx_scale_new;
  ret["integrated"] = integrated_scaled;

  return ret;

}

