// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#include "FermionicBasis.h"

void FermionicBasis::generate_combinations(
    BasisElement& current, std::size_t first_orbital, std::size_t depth,
    std::size_t max_depth) {
  if (depth == max_depth) {
    if (m_basis_filter->filter(current)) {
      insert(current);
    }
    return;
  }

  for (std::size_t i = first_orbital; i < m_orbitals; i++) {
    for (int spin_index = 0; spin_index < 2; ++spin_index) {
      Operator::Spin spin = static_cast<Operator::Spin>(spin_index);
      if (current.empty() || current.back().orbital() < i ||
          (m_allow_double_occupancy &&
           (current.back().orbital() == i && spin > current.back().spin()))) {
        current.push_back(Operator(
            Operator::Type::Creation, Operator::Statistics::Fermion, spin, i));
        generate_combinations(current, i, depth + 1, max_depth);
        current.pop_back();
      }
    }
  }
}
