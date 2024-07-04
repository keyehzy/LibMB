// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#include <armadillo>  //  for eigensolver
#include <iostream>

#include "Assert.h"
#include "BasisFilter.h"
#include "FermionicBasis.h"
#include "Model.h"

using enum Operator::Type;        // for Creation, Annihilation
using enum Operator::Statistics;  // for Fermion
using enum Operator::Spin;        // Up, Down

class HubbardChain : public Model {
 public:
  HubbardChain(double mu, double t, double u, size_t n)
      : m_mu(mu), m_t(t), m_u(u), m_size(n) {}

  ~HubbardChain() override {}

 private:
  Expression hamiltonian() const override {
    Expression result;

    // Chemical potential and Hopping
    for (Operator::Spin spin : {Up, Down}) {
      for (std::size_t i = 0; i < m_size; i++) {
        result += -m_mu * density<Fermion>(spin, i);
        result += -m_t * hopping<Fermion>(spin, i, (i + 1) % m_size);
      }
    }

    // Hubbard U
    for (size_t i = 0; i < m_size; i++) {
      result += m_u * density_density<Fermion>(Up, i, Down, i);
    }

    return result;
  }

  double m_mu;
  double m_t;
  double m_u;
  size_t m_size;
};

int main() {
  const std::size_t size = 2;
  const std::size_t particles = 2;
  const double mu = 0.0;
  const double t = 1.0;
  const double u = 2.0;

  // Create the model
  HubbardChain model(mu, t, u, size);

  // Construct a basis
  FermionicBasis basis(size, particles);

  // Compute matrix elements
  arma::SpMat<arma::cx_double> m(basis.size(), basis.size());
  model.compute_matrix_elements(basis, m);
  LIBMB_ASSERT(m.is_hermitian());

  // Compute ground state using, e.g. Armadillo library
  arma::cx_vec eigval;
  arma::cx_mat eigvec;
  const std::size_t eigval_count = 1;
  arma::eigs_gen(eigval, eigvec, m, eigval_count, "sr");

  // Perform some further analysis here...
}
