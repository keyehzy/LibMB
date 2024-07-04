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
    for (const auto& [basis_element, basis_index] : basis.elements()) {
      Expression::ExpressionMap product =
          NormalOrderer(hamilt * basis_element).terms();
      std::erase_if(product, [&](const auto& item) {
        return !basis.contains(item.first);
      });
      for (const auto& [term, coeff] : product) {
        std::size_t term_index = basis.index(term);
#pragma omp critical
        mat(basis_index, term_index) = coeff;
      }
    }
  }

 protected:
  Model() = default;

 private:
  virtual Expression hamiltonian() const = 0;
};
