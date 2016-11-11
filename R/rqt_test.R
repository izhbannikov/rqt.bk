
data <- read.table("~/Projects/rqt/data/phengen2.dat", skip=1)
y <- data$V1
newgeno <- data[, 2:dim(data)[2]]
covadat <- NULL
STT <- 0.2
weight <- FALSE
cumvar.threshold <- 90
reg.family <- "binomial"
####
