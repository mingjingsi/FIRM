/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2014-2016 Gad Abraham
 * All rights reserved.
 */

#include "randompca.h"
#include "util.h"
#include "svdwide.h"
#include "svdtall.h"

void RandomPCA::pca_fast(MatrixXd& X, unsigned int block_size,
                         unsigned int ndim, unsigned int maxiter,
                         double tol, long seed, bool do_loadings)
{
  unsigned int N, p;
  
  // X_meansd = standardise(X, stand_method_x, verbose);
  N = X.rows();
  p = X.cols();
  
  SVDWide op(X, verbose);
  Spectra::SymEigsSolver<double,
                         Spectra::LARGEST_ALGE, SVDWide> eigs(&op, ndim, ndim * 2 + 1);
  
  eigs.init();
  eigs.compute(maxiter, tol);
  
  double div = 1;
  if(divisor == DIVISOR_N1)
    div = N - 1;
  else if(divisor == DIVISOR_P)
    div = p;
  
  if(eigs.info() == Spectra::SUCCESSFUL)
  {
    U = eigs.eigenvectors();
    // Note: _eigenvalues_, not singular values
    d = eigs.eigenvalues().array() / div;
    if(do_loadings)
    {
      VectorXd s = d.array().sqrt().inverse() / sqrt(div);
      V.noalias() = X.transpose() * U * s.asDiagonal();
    }
    trace = X.array().square().sum() / div;
    pve = d / trace;
    Px = U * d.array().sqrt().matrix().asDiagonal();
    
    verbose && STDOUT << timestamp() << "GRM trace: " << trace << std::endl;
  }
  else
  {
    throw new std::runtime_error(
        std::string("Spectra eigen-decomposition was not successful")
    + ", status: " + std::to_string(eigs.info()));
  }
}