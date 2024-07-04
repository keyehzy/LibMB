// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#include "Basis.h"

#include <sstream>

#include "Assert.h"

void prepare_up_and_down_representation(
    const BasisElement& element, std::vector<int>& up, std::vector<int>& down) {
  for (const auto& o : element) {
    (o.type() == Operator::Type::Annihilation
         ? (o.spin() == Operator::Spin::Up ? up[o.orbital()]--
                                           : down[o.orbital()]--)
         : (o.spin() == Operator::Spin::Up ? up[o.orbital()]++
                                           : down[o.orbital()]++));
  }
}

static constexpr std::string_view unicode_empty_cell = "  ";
static constexpr std::string_view unicode_up_arrow_cell = "\u2191 ";
static constexpr std::string_view unicode_down_arrow_cell = " \u2193";
static constexpr std::string_view unicode_double_occ_cell = "\u2191\u2193";

std::string state_string(const BasisElement& element, std::size_t orbitals) {
  std::vector<int> up(orbitals, 0);
  std::vector<int> down(orbitals, 0);
  prepare_up_and_down_representation(element, up, down);

  std::stringstream out;
  out << "|";

  for (std::size_t i = 0; i < orbitals; i++) {
    switch (2 * up[i] + down[i]) {
      case 0:
        out << unicode_empty_cell;
        break;
      case 1:
        out << unicode_down_arrow_cell;
        break;
      case 2:
        out << unicode_up_arrow_cell;
        break;
      case 3:
        out << unicode_double_occ_cell;
        break;
      default:
        LIBMB_UNREACHABLE();
        break;
    }
    if (i < orbitals - 1) {
      out << ",";
    }
  }
  out << ">";
  return out.str();
}

void Basis::generate_combinations(
    BasisElement& current, size_t first_orbital, size_t depth,
    size_t max_depth) {
  if (m_basis_filter->filter(current)) {
    insert(current);
  }

  if (depth == max_depth) {
    return;
  }

  for (size_t orbital_index = first_orbital; orbital_index < m_orbitals;
       orbital_index++) {
    // TODO: bosonic operator should be integer spin
    auto spin = Operator::Spin::Up;
    if (current.empty() || current.back().orbital() <= orbital_index) {
      // In the generic case, we just create a basis with (single-spin) bosons
      // with all possible particles configurations (0..orbitals)
      current.push_back(Operator(
          Operator::Type::Creation, Operator::Statistics::Boson, spin,
          orbital_index));
      generate_combinations(current, orbital_index, depth + 1, max_depth);
      current.pop_back();
    }
  }
}

void Basis::generate_basis() {
  BasisElement current;
  current.reserve(m_particles);
  generate_combinations(current, 0, 0, m_particles);
}
