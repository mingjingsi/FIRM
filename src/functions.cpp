#include "functions.hpp"
#include <RcppEigen.h>

#include "randompca.h"

mat scale_PCA(mat x, const int dims){

  mat x_scaled = x.each_col()/stddev(x, 0, 1).as_col();
    
  mat xT = x_scaled.t();
    
  Eigen::MatrixXd Xm = Eigen::Map<Eigen::MatrixXd>(xT.memptr(), xT.n_rows, xT.n_cols);
    
  RandomPCA rpca;
  rpca.stand_method_x = 0;
  rpca.divisor = 0;
  rpca.verbose = FALSE;
    
  rpca.pca_fast(Xm, 0, dims, 1e2, 1e-4, 1, FALSE);
    
  Rcpp::NumericMatrix P(Rcpp::wrap(rpca.Px));
  mat embedding = Rcpp::as<arma::mat>(P);

  return embedding;
}

void join_scale_PCA_sep(mat& embedding_SS2, mat& embedding_tenx, mat SS2, mat tenx, const int dims){

  mat all_matrix = join_rows(SS2, tenx);
  mat embedding = scale_PCA(all_matrix, dims);
  embedding_SS2 = embedding.rows(0, SS2.n_cols-1);
  embedding_tenx = embedding.rows(SS2.n_cols, SS2.n_cols+tenx.n_cols-1);

}

mat integrated_scale_PCA(mat x, const int dims){

  mat x_scaled = x.each_col()/stddev(x, 0, 1).as_col();
  x_scaled(find(x_scaled > 10)).fill(10);

  mat xT = x_scaled.t();
    
  Eigen::MatrixXd Xm = Eigen::Map<Eigen::MatrixXd>(xT.memptr(), xT.n_rows, xT.n_cols);
    
  RandomPCA rpca;
  rpca.stand_method_x = 0;
  rpca.divisor = 0;
  rpca.verbose = FALSE;
    
  rpca.pca_fast(Xm, 0, dims, 1e2, 1e-4, 1, FALSE);
    
  Rcpp::NumericMatrix P(Rcpp::wrap(rpca.Px));
  mat embedding = Rcpp::as<arma::mat>(P);

  return embedding;
}

mat integrated_fill(mat SS2, mat tenx, uword gene_all_num, uword nSS2, uword ntenx, uvec gene_all_ind_SS2, uvec gene_all_ind_tenx){
  mat integrated = zeros<mat>(gene_all_num, nSS2 + ntenx);
  mat integrated_SS2 = zeros<mat>(gene_all_num, nSS2);
  mat integrated_tenx = zeros<mat>(gene_all_num, ntenx);

  integrated_SS2.rows(gene_all_ind_SS2-1) = SS2;
  integrated_tenx.rows(gene_all_ind_tenx-1) = tenx;
  integrated.cols(0, nSS2-1) = integrated_SS2;
  integrated.cols(nSS2, nSS2+ntenx-1) = integrated_tenx;

  return integrated;
}

mat integrated_scale_fill(mat SS2, mat tenx, vec sd_SS2, vec sd_tenx, uword gene_all_num, uword nSS2, uword ntenx,
                          uvec gene_all_ind_SS2, uvec gene_all_ind_tenx){
  mat integrated = zeros<mat>(gene_all_num, nSS2 + ntenx);
  mat integrated_SS2 = zeros<mat>(gene_all_num, nSS2);
  mat integrated_tenx = zeros<mat>(gene_all_num, ntenx);

  integrated_SS2.rows(gene_all_ind_SS2-1) = SS2.each_col()/(sd_SS2 + (sd_SS2 == 0));
  integrated_tenx.rows(gene_all_ind_tenx-1) = tenx.each_col()/(sd_tenx + (sd_tenx == 0));
  integrated.cols(0, nSS2-1) = integrated_SS2;
  integrated.cols(nSS2, nSS2+ntenx-1) = integrated_tenx;
    
  return integrated;
}

mat integrated_fill_hvg(mat SS2, mat tenx, uword gene_all_num, uword nSS2, uword ntenx, uvec gene_all_ind_SS2, uvec gene_all_ind_tenx,
                        uvec gene_all_hvg_ind){
  mat integrated = zeros<mat>(gene_all_num, nSS2 + ntenx);
  mat integrated_SS2 = zeros<mat>(gene_all_num, nSS2);
  mat integrated_tenx = zeros<mat>(gene_all_num, ntenx);

  integrated_SS2.rows(gene_all_ind_SS2-1) = SS2;
  integrated_tenx.rows(gene_all_ind_tenx-1) = tenx;
  integrated.cols(0, nSS2-1) = integrated_SS2;
  integrated.cols(nSS2, nSS2+ntenx-1) = integrated_tenx;

  return integrated.rows(gene_all_hvg_ind-1);
}

mat integrated_scale_fill_hvg(mat SS2, mat tenx, vec sd_SS2, vec sd_tenx, uword gene_all_num, uword nSS2, uword ntenx,
                              uvec gene_all_ind_SS2, uvec gene_all_ind_tenx, uvec gene_all_hvg_ind){
  mat integrated = zeros<mat>(gene_all_num, nSS2 + ntenx);
  mat integrated_SS2 = zeros<mat>(gene_all_num, nSS2);
  mat integrated_tenx = zeros<mat>(gene_all_num, ntenx);

  integrated_SS2.rows(gene_all_ind_SS2-1) = SS2.each_col()/(sd_SS2 + (sd_SS2 == 0));
  integrated_tenx.rows(gene_all_ind_tenx-1) = tenx.each_col()/(sd_tenx + (sd_tenx == 0));
  integrated.cols(0, nSS2-1) = integrated_SS2;
  integrated.cols(nSS2, nSS2+ntenx-1) = integrated_tenx;

  return integrated.rows(gene_all_hvg_ind-1);
}

mat integrated_fill_hvg_PCA(mat SS2, mat tenx, uword gene_all_num, uword nSS2, uword ntenx, uvec gene_all_ind_SS2, uvec gene_all_ind_tenx,
                            uvec gene_all_hvg_ind, int dims){

  mat integrated_hvg = integrated_fill_hvg(SS2, tenx, gene_all_num, nSS2, ntenx, gene_all_ind_SS2, gene_all_ind_tenx, gene_all_hvg_ind);

  mat integrated_embedding = integrated_scale_PCA(integrated_hvg, dims);

  return integrated_embedding;
}

mat integrated_scale_fill_hvg_PCA(mat SS2, mat tenx, vec sd_SS2, vec sd_tenx, uword gene_all_num, uword nSS2, uword ntenx, uvec gene_all_ind_SS2, uvec gene_all_ind_tenx,
                            uvec gene_all_hvg_ind, int dims){

  mat integrated_hvg = integrated_scale_fill_hvg(SS2, tenx, sd_SS2, sd_tenx, gene_all_num, nSS2, ntenx, gene_all_ind_SS2, gene_all_ind_tenx, gene_all_hvg_ind);

  mat integrated_embedding = integrated_scale_PCA(integrated_hvg, dims);

  return integrated_embedding;
}


void pair_info(uvec& SS2_paired_name, uvec& tenx_paired_name, uvec& num_SS2, uvec& num_tenx,
               uvec num_SS2_ini, uvec num_tenx_ini, uvec SS2_paired_name_ini, uvec tenx_paired_name_ini){
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

      umat num_tenx_mat = reshape(num_tenx_ini, 1, num_tenx_ini.n_elem);
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
      if ((sum(num_SS21)/sum(num_SS2_ini) < 0.5) || (sum(num_tenx1)/sum(num_tenx_ini) < 0.25)){
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
      if ((sum(num_SS22)/sum(num_SS2_ini) < 0.5) || (sum(num_tenx2)/sum(num_tenx_ini) < 0.25)){
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
}

int check_merge(mat SS2_scale, vec sd_SS2_tmp_new, mat tenx_scale, int dims, vec SS2_FindClusters_current, vec tenx_FindClusters_current,
                uword id_SS2, uword id_tenx, vec quantile_def){
  mat SS2_scale_tmp_new = SS2_scale.each_col()/(sd_SS2_tmp_new + (sd_SS2_tmp_new == 0));

  mat all_SS2_tmp = zeros<mat>(SS2_scale_tmp_new.n_cols, dims);
  mat all_tenx_tmp = zeros<mat>(tenx_scale.n_cols, dims);
  join_scale_PCA_sep(all_SS2_tmp, all_tenx_tmp, SS2_scale_tmp_new, tenx_scale, dims);

  // calculate the new center of each cluster based on cell embeddings of PC
  rowvec SS2_cluster_center = mean(all_SS2_tmp.rows(find(SS2_FindClusters_current == id_SS2)), 0);
  rowvec tenx_cluster_center = mean(all_tenx_tmp.rows(find(tenx_FindClusters_current == id_tenx)), 0);

  mat all_SS2_tmp_i = all_SS2_tmp.rows(find(SS2_FindClusters_current == id_SS2));
  double SS2_cluster_dist_quantile = as_scalar(quantile(sum(square(all_SS2_tmp_i.each_row() - SS2_cluster_center), 1).as_col(), quantile_def));

  mat all_tenx_tmp_i = all_tenx_tmp.rows(find(tenx_FindClusters_current == id_tenx));
  double tenx_cluster_dist_quantile = as_scalar(quantile(sum(square(all_tenx_tmp_i.each_row() - tenx_cluster_center), 1).as_col(), quantile_def));

  return (sum(square(SS2_cluster_center - tenx_cluster_center)) < SS2_cluster_dist_quantile)+
      (sum(square(SS2_cluster_center - tenx_cluster_center)) < tenx_cluster_dist_quantile);
}

int check_merge(mat SS2_scale, mat tenx_scale, vec sd_tenx_tmp_new, int dims, vec SS2_FindClusters_current, vec tenx_FindClusters_current,
                uword id_SS2, uword id_tenx, vec quantile_def){

  mat tenx_scale_tmp_new = tenx_scale.each_col()/(sd_tenx_tmp_new + (sd_tenx_tmp_new == 0));

  mat all_SS2_tmp = zeros<mat>(SS2_scale.n_cols, dims);
  mat all_tenx_tmp = zeros<mat>(tenx_scale_tmp_new.n_cols, dims);
  join_scale_PCA_sep(all_SS2_tmp, all_tenx_tmp, SS2_scale, tenx_scale_tmp_new, dims);

  // calculate the new center of each cluster based on cell embeddings of PC
  rowvec SS2_cluster_center = mean(all_SS2_tmp.rows(find(SS2_FindClusters_current == id_SS2)), 0);
  rowvec tenx_cluster_center = mean(all_tenx_tmp.rows(find(tenx_FindClusters_current == id_tenx)), 0);

  mat all_SS2_tmp_i = all_SS2_tmp.rows(find(SS2_FindClusters_current == id_SS2));
  double SS2_cluster_dist_quantile = as_scalar(quantile(sum(square(all_SS2_tmp_i.each_row() - SS2_cluster_center), 1).as_col(), quantile_def));

  mat all_tenx_tmp_i = all_tenx_tmp.rows(find(tenx_FindClusters_current == id_tenx));
  double tenx_cluster_dist_quantile = as_scalar(quantile(sum(square(all_tenx_tmp_i.each_row() - tenx_cluster_center), 1).as_col(), quantile_def));

  return (sum(square(SS2_cluster_center - tenx_cluster_center)) < SS2_cluster_dist_quantile)+
    (sum(square(SS2_cluster_center - tenx_cluster_center)) < tenx_cluster_dist_quantile);
}

