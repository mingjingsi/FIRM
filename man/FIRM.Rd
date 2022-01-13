\name{FIRM}
\alias{FIRM}
\title{
  FIRM
}
\description{
FIRM integration.
}
\usage{
FIRM(SS2, tenx, 
     hvg1, hvg2, 
     dims, 
     all_genes = FALSE, 
     res_seq_SS2 = seq(0.1, 2, 0.1), res_seq_tenx = seq(0.1, 2, 0.1),
     coreNum = 1, 
     verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SS2}{
     Scaled data for dataset 1, e.g., a SS2 dataset, after standard preprocessing workflow using function \code{prep_data}.
}
  \item{tenx}{
     Scaled data for dataset 2, e.g., a 10X dataset, after standard preprocessing workflow using function \code{prep_data}.
}
  \item{hvg1}{
     Highly variable genes for dataset 1 calculated using function \code{prep_data}.
}
  \item{hvg2}{
     Highly variable genes for dataset 2 calculated using function \code{prep_data}.
}
  \item{dims}{
     Number of PCs.
}
  \item{all_genes}{
     If \code{FALSE}, the integrated data contains only the overlapped hvgs in two datasets. If \code{TRUE}, the integrated data contains the union of all genes in two datasets. Default is \code{FALSE}.
}
  \item{res_seq_SS2}{
     The search search of resolution parameter for dataset 1. Default is \code{seq(0.1, 2, 0.1)}.
}
  \item{res_seq_tenx}{
     The search search of resolution parameter for dataset 2. Default is \code{seq(0.1, 2, 0.1)}.
}
  \item{coreNum}{
     Number of cores used. Default is 1.
  }
  \item{verbose}{
     Whether or not to output the mixing metrics in resolution parameter searching. Defalut is FALSE.
  }
  
}
\details{
  \code{FIRM} performs integration for two datasets. It requires to provide scaled data and hvgs for each dataset using \code{prep_data}.
}
\value{
Integrated data or list of integrated data and mixing metrics in resolution parameter searching.
  \item{integrated}{
     Integrated data.
}
  \item{Metric_FIRM}{
     Mixing metrics in resolution parameter searching.
}
  \item{Metric_PCA}{
     Mixing metrics for baseline model.
}
}

\references{
%% ~put references to the literature/web site here ~
}
\author{
Jingsi Ming
}
\seealso{
 \code{\link{prep_data}}, \code{\link{post}}, \code{\link{test_beta}}
}
\examples{
data("ExampleData")

prep_SS2 <- prep_data(SS2)
Dataset1 <- prep_SS2$Dataset
hvg1 <- prep_SS2$hvg

prep_tenx <- prep_data(tenx)
Dataset2 <- prep_tenx$Dataset
hvg2 <- prep_tenx$hvg

dims <- 15
coreNum <- 4
Dataset <- FIRM(Dataset1, Dataset2, hvg1, hvg2, dims = dims, coreNum = coreNum)
}