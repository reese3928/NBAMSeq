NBAMSeq
===========

`NBAMSeq` is a Bioconductor package for differential expression analysis based on negative binomial additive model. 

The Bioconductor package can be found [here](https://bioconductor.org/packages/release/bioc/html/NBAMSeq.html).    
The Bioconductor package vignette can be found [here](https://bioconductor.org/packages/release/bioc/vignettes/NBAMSeq/inst/doc/NBAMSeq-vignette.html).   
The `NBAMSeq` paper can be found [here](https://doi.org/10.1186/s12859-020-3506-x).

Package installation
------------
```{r}    
if (!requireNamespace("BiocManager", quietly = TRUE))  install.packages("BiocManager")    
BiocManager::install("NBAMSeq")    
```

Citation
------------
Ren, X., & Kuan, P. F. (2020). Negative binomial additive model for RNA-Seq data analysis. BMC Bioinformatics, 21(1), 1-15.

@article{ren2020negative,    
title={Negative binomial additive model for RNA-Seq data analysis},    
author={Ren, Xu and Kuan, Pei-Fen},    
journal={BMC Bioinformatics},    
volume={21},    
number={1},    
pages={1--15},    
year={2020},    
publisher={BioMed Central}    
}
