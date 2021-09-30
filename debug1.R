#library(snedata)
#
#if (file.exists("mnist_digits.rds")) {
#    mnist <- readRDS("mnist_digits.rds")
#} else {
#    mnist <- download_mnist()
#    saveRDS(mnist, file="mnist_digits.rds")
#}
#
#mat <- as.matrix(mnist[,-ncol(mnist)])
#storage.mode(mat) <- "double"
#
#library(uwot)
#debug(uwot:::spectral_init)
#system.time(ref <- umap(mat, n_epochs=500, n_threads=1))

library(Rcpp)
sourceCpp("debug1.cpp")

library(Matrix)
test <- readRDS("debug1.rds")
run_shift(test@x, test@i, test@p, order=nrow(test), 3)
