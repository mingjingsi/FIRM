library(Seurat)
library(RANN)

Mixing_Metric <- function(embedding, dataset_list, k = 5, max.k = 300) {
  set.seed(0)
  nn <- nn2(embedding, k = max.k)$nn.idx[, -1]

  mixing <- sapply(
    X = 1:nrow(x = nn),
    FUN = function(x) {
      sapply(X = unique(dataset_list), FUN = function(y) {
        which(x = dataset_list[nn[x, ]] == y)[k]
      })
    }
  )

  mixing[is.na(x = mixing)] <- max.k
  mixing <- apply(
    X = mixing,
    MARGIN = 2,
    FUN = median
  )

  return(mixing)
}


Local_Struct <- function(SS2, tenx, integrated, dims = 30, neighbors = 20) {
  nn.SS2 <- nn2(SS2@reductions$pca@cell.embeddings[, 1:dims], k = neighbors)$nn.idx
  nn.tenx <- nn2(tenx@reductions$pca@cell.embeddings[, 1:dims], k = neighbors)$nn.idx

  nn.corrected.SS2 <- nn2(integrated@reductions$pca@cell.embeddings[which(integrated$dataset == unique(SS2$dataset)), 1:dims], k = neighbors)$nn.idx
  nn.corrected.tenx <- nn2(integrated@reductions$pca@cell.embeddings[which(integrated$dataset == unique(tenx$dataset)), 1:dims], k = neighbors)$nn.idx

  local.struct.SS2 <- sapply(X = 1:nrow(x = nn.SS2), FUN = function(x) {
    length(x = intersect(x = nn.SS2[x, ], y = nn.corrected.SS2[x, ])) / neighbors})
  local.struct.tenx <- sapply(X = 1:nrow(x = nn.tenx), FUN = function(x) {
    length(x = intersect(x = nn.tenx[x, ], y = nn.corrected.tenx[x, ])) / neighbors})

  local.struct <- list(SS2 = local.struct.SS2, tenx = local.struct.tenx)
  return(local.struct)
}

label_trans <- function(integrated, ref = "10X", query = "SS2", k = 10){
  set.seed(0)
  anno_tenx <- integrated$annotation[which(integrated$dataset == ref)]

  nn <- nn2(integrated@reductions$pca@cell.embeddings[which(integrated$dataset == ref), ],
            integrated@reductions$pca@cell.embeddings[which(integrated$dataset == query), ],
            k = k)

  nn_anno <- matrix(anno_tenx[as.vector(nn$nn.idx)], nrow = nrow(nn$nn.idx), ncol = ncol(nn$nn.idx))

  nn_table <- apply(nn_anno, 1, table)

  pred_anno_SS2 <- as.character(lapply(nn_table, function(x) names(x)[which.max(x)]))
  
  names(pred_anno_SS2) <- colnames(integrated)[which(integrated$dataset == query)]                                  

  return(pred_anno = pred_anno_SS2)
}

match_score <- function(integrated, ref = "10X", query = "SS2", k = 10){
  set.seed(0)
  anno_tenx <- integrated$annotation[which(integrated$dataset == ref)]

  nn <- nn2(integrated@reductions$pca@cell.embeddings[which(integrated$dataset == ref), ],
            integrated@reductions$pca@cell.embeddings[which(integrated$dataset == query), ],
            k = k)

  nn_SS2 <- nn2(integrated@reductions$pca@cell.embeddings[which(integrated$dataset == query), ],
                integrated@reductions$pca@cell.embeddings[which(integrated$dataset == query), ],
                k = k)

  metric <- rowMeans(nn$nn.dists)/rowMeans(nn_SS2$nn.dists[, -1])

  names(metric) <- rownames(integrated@reductions$pca@cell.embeddings[which(integrated$dataset == query), ])

  return(metric = metric)
}

label_trans_prob <- function(integrated, ref = "10X", query = "SS2", k = 10){
  set.seed(0)
  anno_tenx <- integrated$annotation[which(integrated$dataset == ref)]

  nn <- nn2(integrated@reductions$pca@cell.embeddings[which(integrated$dataset == ref), ],
            integrated@reductions$pca@cell.embeddings[which(integrated$dataset == query), ],
            k = k)

  nn_anno <- matrix(anno_tenx[as.vector(nn$nn.idx)], nrow = nrow(nn$nn.idx), ncol = ncol(nn$nn.idx))

  nn_table <- apply(nn_anno, 1, table)

  nn_prob <- lapply(nn_table, function(x) x/sum(x))

  anno_SS2_prob <- matrix(0, length(nn_prob), length(unique(anno_tenx)))
  rownames(anno_SS2_prob) <- rownames(integrated@reductions$pca@cell.embeddings[which(integrated$dataset == query), ])
  colnames(anno_SS2_prob) <- unique(anno_tenx)

  for (i in 1:length(nn_prob)){
    anno_SS2_prob[i, names(nn_prob[[i]])] <- nn_prob[[i]]
  }

  return(anno_SS2_prob = anno_SS2_prob)
}


