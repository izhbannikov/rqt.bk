# rqt: an R-package for gene-level meta analysis

## Installation

```
devtools::install_github("izhbannikov/rqt")

```

## Usage

###Single dataset

```
library(rqt)
data <- read.table(system.file("data/phengen2.dat",package="rqt"), skip=1)
pheno <- data[,1]
geno <- data[, 2:dim(data)[2]]
rqtobj <- rqtClass(phenotype=pheno, genotype=geno)
res <- rQTest(rqtobj)
print(res@results)
```

### Multiple datasets (meta analysis)
```
library(rqt)
data1 <- read.table(system.file("data/phengen2.dat",package="rqt"), skip=1)
obj1 <- rqtClass(phenotype=data1[,1], genotype=data1[, 2:dim(data1)[2]])
data2 <- read.table(system.file("data/phengen3.dat",package="rqt"), skip=1)
obj2 <- rqtClass(phenotype=data2[,1], genotype=data2[, 2:dim(data2)[2]])
data3 <- read.table(system.file("data/phengen.dat",package="rqt"), skip=1)
obj3 <- rqtClass(phenotype=data3[,1], genotype=data3[, 2:dim(data3)[2]])
res <- rQTestMeta(list(obj1, obj2, obj3))
print(res@results)
```
