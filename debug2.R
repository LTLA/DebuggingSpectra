## Generating some data.
#set.seed(10)
#y <- matrix(rnorm(50, sd=2), ncol=10)
#id <- sample(nrow(y), 500, replace=TRUE)
#mat <- y[id,] + matrix(rnorm(5000), ncol=10)
#
#library(FNN)
#res <- FNN::get.knn(mat, k=15)
#library(uwot)
#debug(uwot:::normalized_laplacian_init)
#system.time(ref <- umap(mat, n_epochs=500, n_threads=1))

library(Rcpp)
sourceCpp("debug2.cpp")

library(Matrix)
test <- readRDS("debug2.rds")
run_shift(test@x, test@i, test@p, order=nrow(test), 2)
