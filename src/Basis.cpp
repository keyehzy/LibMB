// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#include "Basis.h"

#include <cassert>
#include <sstream>

void prepare_up_and_down_representation(const Basis::BasisElement& element,
                                        std::vector<int>& up,
                                        std::vector<int>& down) {
  for (auto it = element.rbegin(); it != element.rend(); it++) {
    const Operator& o = *it;
    if (o.type() == Operator::Type::ANNIHILATION) {
      if (o.spin() == Operator::Spin::UP) {
        up[o.orbital()] += -1;
      } else {
        down[o.orbital()] += -1;
      }
    } else {
      if (o.spin() == Operator::Spin::UP) {
        up[o.orbital()] += 1;
      } else {
        down[o.orbital()] += 1;
      }
    }
  }
}

std::string Basis::state_string(const Basis::BasisElement& element) const {
  assert(contains(element));
  std::vector<int> up(m_orbitals, 0);
  std::vector<int> down(m_orbitals, 0);
  prepare_up_and_down_representation(element, up, down);
  std::stringstream out;
  out << "|";
  for (std::size_t i = 0; i < m_orbitals; i++) {
    if (up[i] == 1 && down[i] == 1) {
      out << "\u2191\u2193,";
    } else if (up[i] == 1 && down[i] == 0) {
      out << "\u2191 ,";
    } else if (up[i] == 0 && down[i] == 1) {
      out << " \u2193,";
    } else if (up[i] == 0 && down[i] == 0) {
      out << "  ,";
    } else {
      assert(false);
    }
  }
  out << ">";
  return out.str();
}

void Basis::generate_combinations(BasisElement& current, size_t first_orbital,
                                  size_t depth, size_t max_depth) {
  if (depth == max_depth && (!m_filter || m_filter(current))) {
    m_basis_map.insert(current);
    return;
  }

  for (size_t i = first_orbital; i < m_orbitals; i++) {
    for (int spin_index = 0; spin_index < 2; ++spin_index) {
      Operator::Spin spin = static_cast<Operator::Spin>(spin_index);
      if (current.empty() || current.back().orbital() < i ||
          (current.back().orbital() == i && spin > current.back().spin())) {
        current.push_back(Operator::creation(spin, i));
        generate_combinations(current, i, depth + 1, max_depth);
        current.pop_back();
      }
    }
  }
}

void Basis::generate_basis() {
  BasisElement current;
  current.reserve(m_particles);
  generate_combinations(current, 0, 0, m_particles);
}
