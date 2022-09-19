# FIRM

FIRM is an algorithm for flexible integration of heterogeneous scRNA-seq datasets across multiple tissue types, platforms and experimental batches. FIRM-integrated datasets show accurate mixing of shared cell type identities and superior preservation of original structure without overcorrection, generating robust integrated datasets for downstream exploration and analysis. 

# Installation

To install the development version of FIRM, it's easiest to use the 'devtools' package. Note that FIRM depends on the 'Seurat' pacakge, the 'RANN' package and the 'Rcpp' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.

```
#install.packages("devtools")
library(devtools)
install_github("mingjingsi/FIRM")
```

# Usage

[The 'FIRM' vignette](https://github.com/mingjingsi/FIRM/blob/master/inst/doc/FIRM_package.pdf?raw=true) will provide a good start point for the analysis using FIRM package. Two demos are also provided ([Demo 1](https://mingjingsi.github.io/FIRM_demo1/) and [Demo 2](https://mingjingsi.github.io/FIRM_Demo2/)).


# Reference

Jingsi Ming, Zhixiang Lin, Jia Zhao, Xiang Wan, The Tabula Microcebus Consortium, Can Yang, Angela Ruohao Wu, FIRM: Flexible Integration of single-cell RNA-sequencing data for large-scale Multi-tissue cell atlas datasets. Briefings in Bioinformatics.
https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbac167/6585621
