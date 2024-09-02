// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#include "Models/TJModel.h"

#include <armadillo>
#include <iostream>

#include "FermionicBasis.h"

int main() {
  const std::size_t size = 4;
  const std::size_t particles = 3;  // 3/4 filling
  const double t = 1.0;
  const double J = 0.3;

  // Create the model
  TJModel model(t, J, size);

  // Construct a basis
  FermionicBasis basis(size, particles, /*allow_double_occupancy=*/false);

  // Compute matrix elements
  arma::SpMat<std::complex<double>> m(basis.size(), basis.size());
  model.compute_matrix_elements(basis, m);

  // Compute ground state
  arma::cx_vec eigval;
  arma::cx_mat eigvec;
  arma::eigs_gen(eigval, eigvec, m, 2, "sr");

  double gs_energy = std::real(eigval(0));
  std::cout << "Ground state energy: " << gs_energy << std::endl;

  return 0;
}
