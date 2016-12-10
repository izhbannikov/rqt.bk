# rqt: utilities for gene-level meta-analysis

## Installation

```
devtools::install_github("izhbannikov/rqt")

```

## Usage

###Single dataset

```
library(rqt)
# Loading data and constructing the objects #
data <- data.matrix(read.table(system.file("extdata/test.bin1.dat",
    package="rqt"), header=TRUE))
pheno <- data[,1]
geno <- data[, 2:dim(data)[2]]
colnames(geno) <- paste(seq(1, dim(geno)[2]))
geno.obj <- SummarizedExperiment(geno)
obj <- rqtClass(phenotype=pheno, genotype=geno.obj)
# Analysis #
res <- rQTest(obj, method="pca", out.type = "D")
print(res)
```

### Multiple datasets (meta analysis)
```
library(rqt)
data1 <- data.matrix(read.table(system.file("extdata/phengen2.dat",
                                            package="rqt"), skip=1))
pheno <- data1[,1]
geno <- data1[, 2:dim(data1)[2]]
colnames(geno) <- paste(seq(1, dim(geno)[2]))
geno.obj <- SummarizedExperiment(geno)
obj1 <- rqtClass(phenotype=pheno, genotype=geno.obj)

data2 <- data.matrix(read.table(system.file("extdata/phengen3.dat",
                                            package="rqt"), skip=1))
pheno <- data2[,1]
geno <- data2[, 2:dim(data2)[2]]
colnames(geno) <- paste(seq(1, dim(geno)[2]))
geno.obj <- SummarizedExperiment(geno)
obj2 <- rqtClass(phenotype=pheno, genotype=geno.obj)

data3 <- data.matrix(read.table(system.file("extdata/phengen.dat",
                                            package="rqt"), skip=1))
pheno <- data3[,1]
geno <- data3[, 2:dim(data3)[2]]
colnames(geno) <- paste(seq(1, dim(geno)[2]))
geno.obj <- SummarizedExperiment(geno)
obj3 <- rqtClass(phenotype=pheno, genotype=geno.obj)

res.meta <- rQTestMeta(list(obj1, obj2, obj3))
print(res.meta)
```
