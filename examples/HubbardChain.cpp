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
  void hopping_term(Expression& result) const {
    for (Operator::Spin spin : {Up, Down}) {
      // Chemical potential and Hopping
      for (std::size_t i = 0; i < m_size; i++) {
        result += density<Fermion>(-m_mu, spin, i);
        result += hopping<Fermion>(-m_t, spin, i, (i + 1) % m_size);
      }
    }
  }

  void interaction_term(Expression& result) const {
    for (size_t i = 0; i < m_size; i++) {
      // Hubbard U
      result += density_density<Fermion>(m_u, Up, i, Down, i);
    }
  }

  Expression hamiltonian() const override {
    Expression result;
    hopping_term(result);
    interaction_term(result);
    return result;
  }

  double m_mu;
  double m_t;
  double m_u;
  size_t m_size;
};

int main() {
  const std::size_t size = 8;
  const std::size_t particles = 8;
  const double mu = 0.0;
  const double t = 1.0;
  const double u = 2.0;

  HubbardChain model(mu, t, u, size);

  // Construct a basis with total Spin_z equal to zero
  FermionicBasis basis(size, particles, new TotalSpinFilter(0));

  // Compute matrix elements
  arma::SpMat<arma::cx_double> m(basis.size(), basis.size());
  model.compute_matrix_elements(basis, m);
  LIBMB_ASSERT(m.is_hermitian());

  // Compute ground state using, e.g. Armadillo library
  arma::cx_vec eigval;
  arma::cx_mat eigvec;
  const std::size_t eigval_count = 4;
  arma::eigs_gen(eigval, eigvec, m, eigval_count, "sr");

  std::cout << "Eigenvalues:" << std::endl;
  for (std::size_t i = 0; i < eigval.size(); i++) {
    std::cout << eigval(i).real() << std::endl;
  }

  // Perform some further analysis here...
}
