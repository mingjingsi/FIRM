#ifndef functions_hpp
#define functions_hpp

#include <RcppArmadillo.h>

using namespace arma;

mat scale_PCA(mat x, const int dims);

void join_scale_PCA_sep(mat& embedding_SS2, mat& embedding_tenx, mat SS2, mat tenx, const int dims);

mat integrated_scale_PCA(mat x, const int dims);

mat integrated_fill(mat SS2, mat tenx, uword gene_all_num, uword nSS2, uword ntenx, uvec gene_all_ind_SS2, uvec gene_all_ind_tenx);

mat integrated_scale_fill(mat SS2, mat tenx, vec sd_SS2, vec sd_tenx, uword gene_all_num, uword nSS2, uword ntenx,
                          uvec gene_all_ind_SS2, uvec gene_all_ind_tenx);

mat integrated_fill_hvg_PCA(mat SS2, mat tenx, uword gene_all_num, uword nSS2, uword ntenx, uvec gene_all_ind_SS2, uvec gene_all_ind_tenx,
                            uvec gene_all_hvg_ind, int dims);

mat integrated_scale_fill_hvg_PCA(mat SS2, mat tenx, vec sd_SS2, vec sd_tenx, uword gene_all_num, uword nSS2, uword ntenx, uvec gene_all_ind_SS2, uvec gene_all_ind_tenx,
                                  uvec gene_all_hvg_ind, int dims);

void pair_info(uvec& SS2_paired_name, uvec& tenx_paired_name, uvec& num_SS2, uvec& num_tenx,
               uvec num_SS2_ini, uvec num_tenx_ini, uvec SS2_paired_name_ini, uvec tenx_paired_name_ini);

int check_merge(mat SS2_scale, vec sd_SS2_tmp_new, mat tenx_scale, int dims, vec SS2_FindClusters_current, vec tenx_FindClusters_current,
                uword id_SS2, uword id_tenx, vec quantile_def);

int check_merge(mat SS2_scale, mat tenx_scale, vec sd_tenx_tmp_new, int dims, vec SS2_FindClusters_current, vec tenx_FindClusters_current,
                uword id_SS2, uword id_tenx, vec quantile_def);


#endif
