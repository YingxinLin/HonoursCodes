# Supplementary Codes

Supplementary codes for my Honours thesis.

## Mixture models for gene distributions in single-cell RNA-seq data

+ Four mixture models: 
  + Gamma-Normal (GN): `gammaNormMix()`
  + Constrained Gamma-Normal (cGN): `gammaNormFix()`
  + Gamma-Gamma (GG): `gammaMix()`
  + Constrained Gamma-Gamma (cGG): `gammaMixFix()`
+ Model selection:
  + Select the best model using BIC: `modelSelect()`

## Runing the functions


Load the R packages

```
require(distr)
require(scales)
require(cvTools)
require(edgeR)
```

Download and read the mixtureModels.R file.

```
source("mixtureModels.R")
```

We start from a count matrix, where each row represents a gene, and each column represents a cell. Here, we simulate a count matrix from a negative binomial distribution with 20000 genes and 250 cells as an example.

```
nGene<-20000
nCell<-250
gene_mean <- 2^runif(nGene, -2, 8)
count <- matrix(rnbinom(nGene*nCell, mu=gene_mean, size=10), nrow=nGene)
```

Then we need to transform the raw count matrix into a log-scale counts per million mapped reads (CPM) matrix. Note that, one can also use log-scale FPKM, RPKM or TPM matrix or other normalisation output as the input of mixture model functions.

```
#Transform the raw count matrix
log2cpm<-log2(cpm(count)+1)
```

### Fitting the GN and GG mixture models for Gene 1

```
GN_gene1<-gammaNormMix(log2cpm[1,])
GG_gene1<-gammaMix(log2cpm[1,])
```

### Fitting the GN and GG for the whole dataset to get global parameters

Note that the fittings for the entire dataset need a few minutes...

```
GN_dataset<-gammaNormMix(log2cpm,plot=FALSE)
GG_dataset<-gammaMix(log2cpm,plot=FALSE)
```

### Fitting the contrained mixture models (cGN, cGG) for Gene 1

```
cGN_gene1<-gammaNormFix(log2cpm[1,],alpha_fix=GN_dataset$alpha,beta_fix=GN_dataset$beta)
cGG_gene1<-gammaMixFix(log2cpm[1,],alpha_fix=GG_dataset$alpha1,beta_fix=GG_dataset$beta1)
```

### Select the best model for Gene 1 using BIC

```

modelSelect(log2cpm[1,],alpha_fix=GG_dataset$alpha1,beta_fix=GG_dataset$beta1,
 	alpha_fix_norm=GN_dataset$alpha,beta_fix_norm =GN_dataset$beta,parameter=TRUE)
  
```


## Author

Yingxin Lin

## Acknowledgement

`gammaNormMix()` function is modified from the codes provided by Shila Ghazanfar.
