// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include "Basis.h"
#include "NormalOrderer.h"

class Model {
 public:
  virtual ~Model() = default;

  Model(const Model& other) = delete;
  Model& operator=(const Model& other) = delete;
  Model(Model&& other) = delete;
  Model& operator=(Model&& other) = delete;

  template <typename SpMat>
  void compute_matrix_elements(const Basis& basis, SpMat& mat) const {
    const Expression& hamilt = hamiltonian();
#pragma omp parallel for schedule(dynamic)
    for (const BasisElement& basis_element : basis.elements()) {
      std::size_t basis_index = basis.index(basis_element);
      Expression::ExpressionMap product =
          NormalOrderer(hamilt * basis_element).terms();
      for (const auto& [term, coeff] : product) {
        if (basis.contains(term)) {
          std::size_t term_index = basis.index(term);
#pragma omp critical
          mat(basis_index, term_index) = coeff;
        }
      }
    }
  }

 protected:
  Model() = default;

 private:
  virtual Expression hamiltonian() const = 0;
};
