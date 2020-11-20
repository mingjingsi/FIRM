library(Seurat)
library(RANN)

FIRM <- function(SS2, tenx, hvg1, hvg2, dims, res_low_SS2 = 0.1, res_high_SS2 = 2,
                 res_low_tenx = 0.1, res_high_tenx = 2,
                 quantile_default = 0.75, max.k = 300, rept = 50, coreNum = 1, verbose = FALSE){

  set.seed(0)
  hvg <- intersect(hvg1, hvg2)
  gene_all <- union(rownames(SS2), rownames(tenx))

  SS2_embedding <- RunPCA(SS2[hvg, ], npc = dims, verbose = FALSE)@cell.embeddings
  tenx_embedding <- RunPCA(tenx[hvg, ], npc = dims, verbose = FALSE)@cell.embeddings
  gc()

  SS2_snn <- FindNeighbors(SS2_embedding, force.recalc = TRUE, verbose = FALSE)$snn
  tenx_snn <- FindNeighbors(tenx_embedding, force.recalc = TRUE, verbose = FALSE)$snn

  rm(SS2_embedding)
  rm(tenx_embedding)
  gc()

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

  if (length(res_SS2) == 0){

    integrated_PCA <- matrix(0, length(gene_all), ncol(SS2) + ncol(tenx))
    rownames(integrated_PCA) <- gene_all
    colnames(integrated_PCA) <- c(colnames(SS2), colnames(tenx))
    integrated_PCA[rownames(SS2), 1:ncol(SS2)] <- SS2
    integrated_PCA[rownames(tenx), (ncol(SS2)+1):(ncol(SS2)+ncol(tenx))] <- tenx
    integrated_PCA <- ScaleData(integrated_PCA, do.center = FALSE, verbose = FALSE)

    if (verbose == TRUE){
      return(list(integrated = integrated_PCA))
    }
    else{
      return(integrated_PCA)
    }
  }

  for (i in 1:length(res_seq_tenx)){
    if (is.null(tryCatch(tenx_FindClusters <- FindClusters(tenx_snn, resolution = res_seq_tenx[i], verbose = FALSE)[, 1], error = function(e){}))){
      break
    }

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

  rm(SS2_snn)
  rm(tenx_snn)
  gc()

  SS2_FindClusters <- matrix(as.numeric(as.matrix(SS2_FindClusters)), nrow(SS2_FindClusters), ncol(SS2_FindClusters))
  tenx_FindClusters <- matrix(as.numeric(as.matrix(tenx_FindClusters)), nrow(tenx_FindClusters), ncol(tenx_FindClusters))

  gene_all_num <- length(gene_all)

  tmp <- seq(1, gene_all_num, 1)
  names(tmp) <- gene_all
  gene_all_ind_SS2 <- as.numeric(tmp[rownames(SS2)[which(rownames(SS2) %in% gene_all)]])
  tmp <- seq(1, nrow(SS2), 1)
  names(tmp) <- rownames(SS2)
  hvg_ind_SS2 <- as.numeric(tmp[hvg])

  tmp <- seq(1, gene_all_num, 1)
  names(tmp) <- gene_all
  gene_all_ind_tenx <- as.numeric(tmp[rownames(tenx)[which(rownames(tenx) %in% gene_all)]])
  tmp <- seq(1, nrow(tenx), 1)
  names(tmp) <- rownames(tenx)
  hvg_ind_tenx <- as.numeric(tmp[hvg])

  gene_all_hvg <- which(gene_all %in% hvg)

  if (ncol(SS2) < ncol(tenx)){
    Dataset1 <- SS2
    Dataset2 <- tenx
    hvg_ind1 <- hvg_ind_SS2
    FindClusters1 <- SS2_FindClusters
    hvg_ind2 <- hvg_ind_tenx
    FindClusters2 <- tenx_FindClusters
    gene_all_ind1 <- gene_all_ind_SS2
    gene_all_ind2 <- gene_all_ind_tenx
    res1 <- res_SS2
    res2 <- res_tenx

  } else {
    Dataset1 <- tenx
    Dataset2 <- SS2
    hvg_ind1 <- hvg_ind_tenx
    FindClusters1 <- tenx_FindClusters
    hvg_ind2 <- hvg_ind_SS2
    FindClusters2 <- SS2_FindClusters
    gene_all_ind1 <- gene_all_ind_tenx
    gene_all_ind2 <- gene_all_ind_SS2
    res1 <- res_tenx
    res2 <- res_SS2
  }
  
  dataset_list <- c(rep(1, ncol(Dataset1)), rep(2, ncol(Dataset2)))
  dataset_list_PCA <- c(rep(1, ncol(SS2)), rep(2, ncol(tenx)))
  
  if ((length(res1)*length(res2)) == 1){
    result <- FIRM_res(Dataset1, hvg_ind1, FindClusters1,
                       Dataset2, hvg_ind2, FindClusters2,
                       dims, gene_all_num, gene_all_hvg,
                       gene_all_ind1, gene_all_ind2,
                       quantile_default = quantile_default, rept = rept)

    integrated_PCA <- matrix(0, length(gene_all), ncol(SS2) + ncol(tenx))
    rownames(integrated_PCA) <- gene_all
    colnames(integrated_PCA) <- c(colnames(SS2), colnames(tenx))
    integrated_PCA[rownames(SS2), 1:ncol(SS2)] <- SS2
    integrated_PCA[rownames(tenx), (ncol(SS2)+1):(ncol(SS2)+ncol(tenx))] <- tenx
    integrated_PCA <- ScaleData(integrated_PCA, do.center = FALSE, verbose = FALSE)
    integrated_PCA_embedding <- RunPCA(integrated_PCA[hvg, ], features = hvg, npcs = dims, verbose = FALSE)
    Metric_PCA <- mean(Mixing_Metric(integrated_PCA_embedding@cell.embeddings, dataset_list_PCA, max.k = max.k))
    
    if (all(result$integrated == 0)){
      # print("PCA")
      
      if (verbose == TRUE){
        return(list(integrated = integrated_PCA, Metric_PCA = Metric_PCA))
      }
      else{
        return(integrated_PCA)
      }
    }

    integrated_FIRM <- result$integrated
    rownames(integrated_FIRM) <- gene_all
    colnames(integrated_FIRM) <- c(colnames(Dataset1), colnames(Dataset2))
    integrated_FIRM <- ScaleData(integrated_FIRM, do.center = FALSE, verbose = FALSE)
    integrated_FIRM_embedding <- RunPCA(integrated_FIRM[hvg, ], features = hvg, npcs = dims, verbose = FALSE)
    Metric_FIRM <- mean(Mixing_Metric(integrated_FIRM_embedding@cell.embeddings, dataset_list, max.k = max.k))
    
    gc()

    if(min(Metric_FIRM) >= Metric_PCA){
      # print("PCA")
      
      if (verbose == TRUE){
        return(list(integrated = integrated_PCA, Metric_PCA = Metric_PCA, Metric_FIRM = Metric_FIRM))
      }
      else{
        return(integrated_PCA)
      }
    } else {
      if (verbose == TRUE){
        return(list(integrated = integrated_FIRM, Metric_FIRM = Metric_FIRM, Metric_PCA = Metric_PCA))
      }
      else{
        return(integrated_FIRM)
      }
    }
    
  }
  
  result <- FIRM_res_all(Dataset1, hvg_ind1, FindClusters1,
                         Dataset2, hvg_ind2, FindClusters2,
                         dims, gene_all_num, gene_all_hvg,
                         gene_all_ind1, gene_all_ind2,
                         quantile_default = quantile_default, rept = rept, coreNum = coreNum)

  Metric_PCA <- mean(Mixing_Metric(result$Embedding_PCA, dataset_list, max.k = max.k))

  Metric_FIRM <- matrix(0, length(res1), length(res2))
  rownames(Metric_FIRM) <- res1
  colnames(Metric_FIRM) <- res2

  idx <- 0;
  for (i in 1:length(res1)){
    for (j in 1:length(res2)){
      idx <- idx + 1
      if (all(result$Embedding_FIRM[, , idx] == 0)){
        Metric_FIRM[i, j] <- Metric_PCA
      } else {
        Metric_FIRM[i, j] <- mean(Mixing_Metric(result$Embedding_FIRM[, , idx], dataset_list, max.k = max.k))
      }
    }
  }

  if(min(Metric_FIRM) >= Metric_PCA){
    # print("PCA")
    integrated_PCA <- matrix(0, length(gene_all), ncol(SS2) + ncol(tenx))
    rownames(integrated_PCA) <- gene_all
    colnames(integrated_PCA) <- c(colnames(SS2), colnames(tenx))
    integrated_PCA[rownames(SS2), 1:ncol(SS2)] <- SS2
    integrated_PCA[rownames(tenx), (ncol(SS2)+1):(ncol(SS2)+ncol(tenx))] <- tenx
    integrated_PCA <- ScaleData(integrated_PCA, do.center = FALSE, verbose = FALSE)

    if (verbose == TRUE){
      return(list(integrated = integrated_PCA, Metric_PCA = Metric_PCA, Metric_FIRM = Metric_FIRM))
    }
    else{
      return(integrated_PCA)
    }
  } else{
    i <- which.min(Metric_FIRM) %% length(res1)
    if (i == 0){
      i <- length(res1)
    }
    j <- ceiling(which.min(Metric_FIRM)/length(res1))

    result <- FIRM_res(Dataset1, hvg_ind1, FindClusters1[, i],
                       Dataset2, hvg_ind2, FindClusters2[, j],
                       dims, gene_all_num, gene_all_hvg,
                       gene_all_ind1, gene_all_ind2,
                       quantile_default = quantile_default, rept = rept)

    if (all(result$integrated == 0)){
      # print("PCA")
      integrated_PCA <- matrix(0, length(gene_all), ncol(SS2) + ncol(tenx))
      rownames(integrated_PCA) <- gene_all
      colnames(integrated_PCA) <- c(colnames(SS2), colnames(tenx))
      integrated_PCA[rownames(SS2), 1:ncol(SS2)] <- SS2
      integrated_PCA[rownames(tenx), (ncol(SS2)+1):(ncol(SS2)+ncol(tenx))] <- tenx
      integrated_PCA <- ScaleData(integrated_PCA, do.center = FALSE, verbose = FALSE)

      if (verbose == TRUE){
        return(list(integrated = integrated_PCA, Metric_PCA = Metric_PCA, Metric_FIRM = Metric_FIRM))
      }
      else{
        return(integrated_PCA)
      }
    }

    integrated_FIRM <- result$integrated
    rownames(integrated_FIRM) <- gene_all
    colnames(integrated_FIRM) <- c(colnames(Dataset1), colnames(Dataset2))
    integrated_FIRM <- ScaleData(integrated_FIRM, do.center = FALSE, verbose = FALSE)

    gc()

    if (verbose == TRUE){
      return(list(integrated = integrated_FIRM, Metric_FIRM = Metric_FIRM, Metric_PCA = Metric_PCA))
    }
    else{
      return(integrated_FIRM)
    }
  }
}

