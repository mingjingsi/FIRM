# FIRM

FIRM is an algorithm for accurate integration of heterogeneous scRNA-seq datasets across multiple platforms, which specifically accounts for the heterogeneity in cell type composition between SS2 and 10X datasets. The integrated datasets generated using FIRM show accurate mixing of shared cell type identities and superior preservation of original structure for each dataset.

# Installation

To install the development version of FIRM, it's easiest to use the 'devtools' package. Note that FIRM depends on the 'Seurat' pacakge, the 'RANN' package and the 'Rcpp' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.

```
#install.packages("devtools")
library(devtools)
install_github("mingjingsi/FIRM")
```

# Reference

Jingsi Ming, Zhixiang Lin, Xiang Wan, Can Yang, Angela Ruohao Wu, FIRM: Fast Integration of single-cell RNA-sequencing data across Multiple platforms
https://www.biorxiv.org/content/10.1101/2020.06.02.129031v1.full.pdf
