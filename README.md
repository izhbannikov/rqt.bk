# rqt: an R-package for gene-level meta analysis

## Installation

```devtools::install_github("izhbannikov/rqt")```

## Usage

###Single dataset

```
library(rqt)
data <- read.table(system.file("data/phengen2.dat",package="rqt"), skip=1)
pheno <- data[,1]
geno <- data[, 2:dim(data)[2]]
res <- rQTest(phenotype=pheno, genotype=geno)
```

### Multiple datasets (meta analysis)
```
library(rqt)
data1 <- read.table(system.file("data/phengen2.dat",package="rqt"), skip=1)
data2 <- read.table(system.file("data/phengen3.dat",package="rqt"), skip=1)
data3 <- read.table(system.file("data/phengen.dat",package="rqt"), skip=1)
# Combining datasets:
indata <- list(data1, data2, data3)
# Gene-level meta-analysis:
res <- rQTestMeta(indata)
res
```
