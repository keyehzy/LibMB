// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#include "Models/HubbardKagome.h"

#include <armadillo>  //  for eigensolver

#include "BasisFilter.h"
#include "FermionicBasis.h"

template <typename Vec>
std::vector<Term> sorted_terms_from_eigvec(
    const Basis& basis, const Vec& eigvec) {
  std::vector<Term> sorted_terms;
  sorted_terms.reserve(basis.size());

  for (std::size_t i = 0; i < basis.size(); i++) {
    sorted_terms.emplace_back(eigvec(i), basis.element(i));
  }

  std::sort(
      sorted_terms.begin(), sorted_terms.end(),
      [](const auto& a, const auto& b) {
        return std::abs(a.coefficient()) > std::abs(b.coefficient());
      });

  return sorted_terms;
}

int main() {
  const std::size_t size = 12;
  const std::size_t particles = 12;
  const double t = -1.0;
  const double u = 0.5;

  HubbardKagome model(t, u, /*periodic=*/false);
  FermionicBasis basis(size, particles, new TotalSpinFilter(0));

  arma::SpMat<arma::cx_double> m(basis.size(), basis.size());
  model.compute_matrix_elements(basis, m);
  LIBMB_ASSERT(m.is_hermitian());

  std::cout << "This is the compact state" << std::endl;
  std::cout << state_string(basis.element(0), basis.orbitals()) << std::endl;

  arma::cx_vec eigval;
  arma::cx_mat eigvec;
  const std::size_t eigval_count = 2;
  arma::eigs_gen(eigval, eigvec, m, eigval_count, "sr");

  for (std::size_t j = 0; j < eigval_count; j++) {
    const arma::cx_vec& ground_state = eigvec.col(j);
    std::vector<Term> sorted_terms =
        sorted_terms_from_eigvec(basis, ground_state);
    std::cout << "Ground state energy" << eigval(j) << std::endl;
    const std::size_t states_to_print = 10;
    for (std::size_t i = 0; i < states_to_print; i++) {
      std::cout << std::fixed
                << std::real(
                       sorted_terms[i].coefficient() *
                       std::conj(sorted_terms[i].coefficient()))
                << "   "
                << state_string(sorted_terms[i].operators(), basis.orbitals())
                << std::endl;
    }
  }
  return 0;
}
