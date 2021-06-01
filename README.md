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

[The 'FIRM' vignette](https://github.com/mingjingsi/FIRM/blob/master/inst/doc/FIRM_package.pdf?raw=true) will provide a good start point for the analysis using FIRM package. Two demos are also provided ([Demo 1](https://drive.google.com/file/d/1nPE2Wjz7Y0BVSFFgIIFyx-u1tBoxzRJu/view?usp=sharing) and [Demo 2](https://drive.google.com/file/d/1tS6UIHB5jhYmOksmS62R29C2s5FgW-lM/view?usp=sharing)).


# Reference

Jingsi Ming, Zhixiang Lin, Xiang Wan, Can Yang, Angela Ruohao Wu, FIRM: Fast Integration of single-cell RNA-sequencing data across Multiple platforms
https://www.biorxiv.org/content/10.1101/2020.06.02.129031v1.full.pdf
