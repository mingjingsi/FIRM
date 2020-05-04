library(Seurat)
library(RANN)

FIRM <- function(SS2, tenx, hvg, dims, res_low_SS2 = 0.1, res_high_SS2 = 2,
                 res_low_tenx = 0.1, res_high_tenx = 2, quantile_default = 0.75, max.k = 300, rept = 50, coreNum = 1, verbose = FALSE){

  set.seed(0)
  gene_all <- union(rownames(SS2), rownames(tenx))

  SS2_embedding <- RunPCA(SS2[hvg, ], npc = dims, verbose = FALSE)@cell.embeddings
  tenx_embedding <- RunPCA(tenx[hvg, ], npc = dims, verbose = FALSE)@cell.embeddings

  SS2_snn <- FindNeighbors(SS2_embedding, force.recalc = TRUE, verbose = FALSE)$snn
  tenx_snn <- FindNeighbors(tenx_embedding, force.recalc = TRUE, verbose = FALSE)$snn

  res_seq_SS2 <- seq(res_low_SS2, res_high_SS2, 0.1)
  res_seq_tenx <- seq(res_low_tenx, res_high_tenx, 0.1)

  res_SS2 <- NULL
  res_tenx <- NULL

  for (i in 1:length(res_seq_SS2)){
    if (is.null(tryCatch(SS2_FindClusters <- FindClusters(SS2_snn, resolution = res_seq_SS2[i], verbose = FALSE)[, 1], error = function(e){}))){
      break
    }

    SS2_FindClusters <- FindClusters(SS2_snn, resolution = res_seq_SS2[i], verbose = FALSE)[, 1]
    num_cluster_SS2 <- length(unique(SS2_FindClusters))

    if (num_cluster_SS2 > 1){
      if (is.null(res_SS2)){
        res_SS2 <- c(res_SS2, res_seq_SS2[i])
        num_cluster_SS2_old <- length(unique(SS2_FindClusters))
      }
      else if (num_cluster_SS2 > num_cluster_SS2_old){
        res_SS2 <- c(res_SS2, res_seq_SS2[i])
      }
    }
    num_cluster_SS2_old <- num_cluster_SS2
  }

  dataset_list <- c(rep(1, ncol(SS2)), rep(2, ncol(tenx)))

  integrated_PCA <- matrix(0, length(gene_all), ncol(SS2) + ncol(tenx))
  rownames(integrated_PCA) <- gene_all
  colnames(integrated_PCA) <- c(colnames(SS2), colnames(tenx))
  integrated_PCA[rownames(SS2), 1:ncol(SS2)] <- SS2
  integrated_PCA[rownames(tenx), (ncol(SS2)+1):(ncol(SS2)+ncol(tenx))] <- tenx
  integrated_PCA <- ScaleData(integrated_PCA, do.center = FALSE, verbose = FALSE)
  # integrated_PCA_embedding <- RunPCA(integrated_PCA[hvg, ], npc = dims, verbose = FALSE)@cell.embeddings

  # Metric_PCA <- mean(MixingMetric(integrated_PCA_embedding, dataset_list, max.k = max.k))

  if (length(res_SS2) == 0){
    if (verbose == TRUE){
      return(list(integrated = integrated_PCA, Metric_PCA = Metric_PCA))
    }
    else{
      return(integrated_PCA)
    }
  }

  for (i in 1:length(res_seq_tenx)){
    tenx_FindClusters <- FindClusters(tenx_snn, resolution = res_seq_tenx[i], verbose = FALSE)[, 1]
    num_cluster_tenx <- length(unique(tenx_FindClusters))

    if (num_cluster_tenx > 1){
      if (is.null(res_tenx)){
        res_tenx <- c(res_tenx, res_seq_tenx[i])
        num_cluster_tenx_old <- length(unique(tenx_FindClusters))
      }
      else if (num_cluster_tenx > num_cluster_tenx_old){
        res_tenx <- c(res_tenx, res_seq_tenx[i])
      }
    }
    num_cluster_tenx_old <- num_cluster_tenx
  }

  res_SS2 <- res_SS2[length(res_SS2):1]
  res_tenx <- res_tenx[length(res_tenx):1]

  SS2_FindClusters <- FindClusters(SS2_snn, resolution = res_SS2, verbose = FALSE)
  tenx_FindClusters <- FindClusters(tenx_snn, resolution = res_tenx, verbose = FALSE)

  SS2_FindClusters <- matrix(as.numeric(as.matrix(SS2_FindClusters)), nrow(SS2_FindClusters), ncol(SS2_FindClusters))
  tenx_FindClusters <- matrix(as.numeric(as.matrix(tenx_FindClusters)), nrow(tenx_FindClusters), ncol(tenx_FindClusters))

  SS2_scale <- SS2[hvg, ]
  tenx_scale <- tenx[hvg, ]

  gene_all_num <- length(gene_all)

  tmp <- seq(1, nrow(SS2), 1)
  names(tmp) <- rownames(SS2)
  gene_all_ind_SS2 <- as.numeric(tmp[gene_all[which(gene_all %in% rownames(SS2))]])
  hvg_ind_SS2 <- as.numeric(tmp[hvg])

  tmp <- seq(1, nrow(tenx), 1)
  names(tmp) <- rownames(tenx)
  gene_all_ind_tenx <- as.numeric(tmp[gene_all[which(gene_all %in% rownames(tenx))]])
  hvg_ind_tenx <- as.numeric(tmp[hvg])

  gene_all_hvg <- which(gene_all %in% hvg)


  result_Metric <- FIRM_res_all(SS2, hvg_ind_SS2, SS2_FindClusters,
                         tenx, hvg_ind_tenx, tenx_FindClusters,
                         dims, gene_all_num, gene_all_hvg,
                         gene_all_ind_SS2, gene_all_ind_tenx,
                         dataset_list, quantile_default = quantile_default, rept = rept, coreNum = coreNum)

  if(min(result_Metric$Metric_FIRM) >= result_Metric$Metric_PCA){
    if (verbose == TRUE){
      return(list(integrated = integrated_PCA, Metric_PCA = Metric_PCA))
    }
    else{
      print("PCA")
      return(integrated_PCA)
    }
  } else{
    i <- which.min(result_Metric$Metric_FIRM) %% length(res_SS2)
    if (i == 0){
      i <- length(res_SS2)
    }
    j <- ceiling(which.min(result_Metric$Metric_FIRM)/length(res_SS2))

    result <- FIRM_res(SS2, hvg_ind_SS2, SS2_FindClusters[, i],
                       tenx, hvg_ind_tenx, tenx_FindClusters[, j],
                       dims, gene_all_num, gene_all_hvg,
                       gene_all_ind_SS2, gene_all_ind_tenx,
                       dataset_list, quantile_default = quantile_default, rept = rept)

    if (length(result) != 3){
      if (verbose == TRUE){
        return(list(integrated = integrated_PCA, Metric_PCA = Metric_PCA))
      }
      else{
        print("PCA")
        return(integrated_PCA)
      }
    }

    if (verbose == TRUE){
      return(list(integrated = result$integrated, Metric_FIRM = result$Metric_FIRM))
    }
    else{
      return(result$integrated)
    }
  }
}


SelectGene <- function(hvg_list, gene_all = NULL, num = 2000){
  K <- length(hvg_list)

  hvg_union <- NULL
  for (k in 1:K){
    hvg_union <- union(hvg_union, hvg_list[[k]])
    if (!is.null(gene_all)){
      hvg_union <- intersect(hvg_union, gene_all)
    }
  }

  if (length(hvg_union) < num){
    return(hvg_union)
  }

  hvg_union_rank <- matrix(0, length(hvg_union), K)
  rownames(hvg_union_rank) <- hvg_union
  for (i in 1:length(hvg_union)){
    for (k in 1:K){
      if (hvg_union[i] %in% hvg_list[[k]]) {
        hvg_union_rank[i, k] <- which(hvg_list[[k]] == hvg_union[i])
      }
    }
  }

  hvg <- hvg_union[which(rowSums(hvg_union_rank != 0) == K)]

  if (length(hvg) >= num){
    return(hvg)
  } else{
    for (k in (K-1):1){
      hvg_add <- hvg_union[which(rowSums(hvg_union_rank != 0) == k)]
      if (length(hvg_add) + length(hvg) >= num){
        hvg_add <- names(sort(apply(hvg_union_rank[hvg_add, ], 1, function(x) median(x[x != 0]))))[1:(num-length(hvg))]
        hvg <- c(hvg, hvg_add)
        return(hvg)
      } else {
        hvg <- c(hvg, hvg_add)
      }
    }
  }
}

