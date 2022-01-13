\name{prep_data}
\alias{FIRM}
\title{
  Data preprocessing for FIRM integration
}
\description{
Data preprocessing for FIRM integration. Standard pre-processing workflow in Seurat is performed, which includes normalization, scaling and feature selection. 
}
\usage{
prep_data(counts, file_path = NULL, file_name = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{counts}{
     Unnormalized data such as raw counts or UMIs.
}
  \item{file_path, file_name}{
     If \code{file_path} and \code{file_name} are provided, the ouput would be saved in the \code{file_path} with name \code{file_name}. Default is \code{NULL}.
}
}
\details{
  \code{prep_data} performs standard pre-processing workflow, which includes normalization, scaling and feature selection.
}
\value{
List of scaled data and highly variable genes.
  \item{Dataset}{
     Scaled data.
}
  \item{hvg}{
     Top 4000 highly variable genes.
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jingsi Ming
}
\seealso{
 \code{\link{FIRM}}
}
\examples{
data("ExampleData")

prep_SS2 <- prep_data(SS2)
Dataset1 <- prep_SS2$Dataset
hvg1 <- prep_SS2$hvg

prep_tenx <- prep_data(tenx)
Dataset2 <- prep_tenx$Dataset
hvg2 <- prep_tenx$hvg

}