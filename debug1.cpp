#include "Spectra/SymEigsShiftSolver.h"
#include "Spectra/MatOp/SparseSymShiftSolve.h"
#include "Spectra/SymEigsSolver.h"
#include "Spectra/MatOp/SparseSymMatProd.h"
#include "Eigen/Sparse"
#include "Rcpp.h"
#include <iostream>

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector run_shift(Rcpp::NumericVector x, Rcpp::IntegerVector i, Rcpp::IntegerVector p, int order, int ndim) {
    std::vector<int> sizes(order);
    for (int c = 0; c < order; ++c) {
        sizes[c] = p[c + 1] - p[c];
    }

    Eigen::SparseMatrix<double> mat(order, order);
    mat.reserve(sizes);

    int counter = 0;
    for (int c = 0; c < order; ++c) {
        for (; counter < p[c + 1]; ++counter) {
            mat.insert(i[counter], c) = x[counter];
        }
    }
    mat.makeCompressed();

    // Finding the smallest eigenvalues & their eigenvectors,
    // using the shift-and-invert mode as recommended.
    const int nobs = mat.rows();
    int nev = std::min(ndim + 1, nobs); // +1 from uwot:::normalized_laplacian_init
    int ncv = std::min(nobs, std::max(2 * nev, 20)); // from RSpectra:::eigs_real_sym. I don't make the rules.

    {
        std::cout << "Setting up the solver..." << std::endl;
        Spectra::SparseSymMatProd<double, Eigen::Upper> op(mat);
        std::cout << "... still setting it up..." << std::endl;
        Spectra::SymEigsSolver<typename std::remove_reference<decltype(op)>::type> eigs(op, nev, ncv); 

        std::cout << "Initializing..." << std::endl;
        eigs.init();
        std::cout << "Computing..." << std::endl;
        eigs.compute(Spectra::SortRule::SmallestMagn, 1000, 1e-20);
        std::cout << "Done!" << std::endl;

        std::cout << eigs.eigenvalues() << std::endl;
    }

    {
        std::cout << "Setting up the solver..." << std::endl;
        Spectra::SparseSymShiftSolve<double, Eigen::Upper> op(mat);
        std::cout << "... still setting it up..." << std::endl;
        Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<double, Eigen::Upper> > eigs(op, nev, ncv, -0.001); // see https://github.com/yixuan/spectra/issues/126

        std::cout << "Initializing..." << std::endl;
        eigs.init();
        std::cout << "Computing..." << std::endl;
        eigs.compute(Spectra::SortRule::LargestMagn);
        std::cout << "Done!" << std::endl;

        std::cout << eigs.eigenvalues() << std::endl;
    }

    return Rcpp::NumericVector(0);
}
