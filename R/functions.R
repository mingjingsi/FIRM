library(Seurat)

prep_data <- function(counts, file_path = NULL, file_name = NULL){

  Dataset <- CreateSeuratObject(counts)
  Dataset <- NormalizeData(Dataset, verbose = FALSE)
  Dataset <- FindVariableFeatures(Dataset, nfeatures = 4000, verbose = FALSE)
  Dataset <- ScaleData(Dataset, features = rownames(Dataset), do.center = FALSE, verbose = FALSE)

  hvg <- Dataset@assays$RNA@var.features
  Dataset <- Dataset@assays$RNA@scale.data

  if (is.null(file_name)){
    return(list(Dataset = Dataset, hvg = hvg))
  } else {
    save(hvg, file = paste(file_path, "/", file_name, "_hvg.RData", sep = ""))
    save(Dataset, file =  paste(file_path, "/", file_name, ".RData", sep = ""))
  }
}

SelectGene <- function(hvg_list, gene_all = NULL, num = 4000){
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

Select_hvg <- function(file_names, file_path, method){
  if (method == "all"){
    hvg_all <- NULL
    for (i in 1:length(file_names)){
      load(paste(file_path, "/", file_names[i], "_hvg.RData", sep = ""))
      hvg_all <- union(hvg_all, hvg)
    }
    return(hvg_all)
  }

  if (method == "top"){
    hvg_all <- NULL
    for (i in 1:length(file_names)){
      load(paste(file_path, "/", file_names[i], "_hvg.RData", sep = ""))
      hvg_all <- c(hvg_all, list(hvg))
    }
    hvg_top <- SelectGene(hvg_all, num = 4000)
    return(hvg_top)
  }
}
