// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#include "BosonicBasis.h"

void BosonicBasis::generate_combinations(
    BasisElement& current, std::size_t first_orbital, std::size_t depth,
    std::size_t max_depth) {
  if (depth == max_depth) {
    if (m_basis_filter->filter(current)) {
      insert(current);
    }
    return;
  }

  for (std::size_t orbital_index = first_orbital; orbital_index < m_orbitals;
       orbital_index++) {
    // TODO: bosonic operator should be integer spin
    Operator::Spin spin = Operator::Spin::Up;
    if (current.empty() || current.back().orbital() <= orbital_index) {
      current.push_back(Operator(
          Operator::Type::Creation, Operator::Statistics::Boson, spin,
          orbital_index));
      generate_combinations(current, orbital_index, depth + 1, max_depth);
      current.pop_back();
    }
  }
}
