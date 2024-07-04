// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#include <armadillo>

#include "Basis.h"
#include "FermionicBasis.h"
#include "Model.h"

class HeisenbergChain : public Model {
 public:
  HeisenbergChain(std::size_t n, double J, double h)
      : m_size{n}, m_J{J}, m_h{h} {}

  ~HeisenbergChain() override {}

  std::size_t size() const { return m_size; }

 private:
  Expression hamiltonian() const override {
    Expression result;

    for (std::size_t i = 0; i < m_size; i++) {
      result += -m_h * spin_z(i);
    }

    for (std::size_t i = 0; i < m_size; i++) {
      result += -m_J * spin_x(i) * spin_x((i + 1) % m_size);
      result += -m_J * spin_y(i) * spin_y((i + 1) % m_size);
      result += -m_J * spin_z(i) * spin_z((i + 1) % m_size);
    }
    return result;
  }

  std::size_t m_size;
  double m_J;
  double m_h;
};

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

static void analysis(
    const HeisenbergChain& model, const FermionicBasis& basis) {
  arma::SpMat<std::complex<double>> m(basis.size(), basis.size());

  model.compute_matrix_elements(basis, m);
  // std::cout << arma::cx_mat(m) << std::endl;
  LIBMB_ASSERT(m.is_hermitian());

  arma::cx_vec eigval;
  arma::cx_mat eigvec;
  bool ok = arma::eigs_gen(eigval, eigvec, m, 1, "sr");

  if (!ok) {
    std::cerr << "Diagonalization failed" << std::endl;
    exit(1);
  }

  const arma::cx_vec& ground_state = eigvec.col(0);
  std::vector<Term> sorted_terms =
      sorted_terms_from_eigvec(basis, ground_state);

  std::cout << "Ground state: " << eigval(0) << std::endl;

  for (std::size_t i = 0; i < 4; i++) {
    std::cout << std::fixed
              << sorted_terms[i].coefficient() *
                     std::conj(sorted_terms[i].coefficient())
              << "  "
              << state_string(sorted_terms[i].operators(), basis.orbitals())
              << std::endl;
  }
}

int main() {
  const int chain_size = 8;
  const double J = 1.0;
  const double small_h_field = 1e-9;
  FermionicBasis basis(chain_size, chain_size, /*allow_double_occupancy=*/true);

  {
    std::cout << "Ferromagnetic Heisenberg Chain" << std::endl;
    std::cout << "Exact ground state: " << -0.25 * chain_size * J << std::endl;
    HeisenbergChain model(chain_size, J, small_h_field);
    analysis(model, basis);
  }

  {
    std::cout << "Antiferromagnetic Heisenberg Chain" << std::endl;
    std::cout << "Exact ground state for infinite chain: "
              << -chain_size * J * (log(2) - 0.25) << std::endl;
    HeisenbergChain model(chain_size, -J, small_h_field);
    analysis(model, basis);
  }

  return 0;
}
