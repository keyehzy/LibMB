// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <algorithm>
#include <unordered_map>
#include <vector>

#include "BasisFilter.h"
#include "Operator.h"
#include "Pointers/NonnullOwnPtr.h"

using BasisElement = std::vector<Operator>;

bool operator<(const BasisElement& lhs, const BasisElement& rhs);

class Basis {
 public:
  virtual ~Basis() = default;

  Basis(std::size_t n, std::size_t m)
      : m_orbitals{n}, m_particles{m}, m_basis_filter{make<BasisFilter>()} {}

  Basis(std::size_t n, std::size_t m, BasisFilter* filter)
      : m_orbitals{n}, m_particles{m}, m_basis_filter{adopt_own(filter)} {}

  std::size_t orbitals() const { return m_orbitals; }

  std::size_t particles() const { return m_particles; }

  bool operator==(const Basis& other) const {
    return m_orbitals == other.m_orbitals && m_particles == other.m_particles &&
           m_basis_set == other.m_basis_set;
  }

  bool operator!=(const Basis& other) const { return !(*this == other); }

  const std::vector<BasisElement>& elements() const { return m_basis_set; }

  const BasisElement& element(std::size_t i) const { return m_basis_set[i]; }

  void insert(const BasisElement& value) {
    auto it = std::lower_bound(m_basis_set.begin(), m_basis_set.end(), value);
    m_basis_set.insert(it, value);
  }

  bool contains(const BasisElement& value) const {
    auto it = std::lower_bound(m_basis_set.begin(), m_basis_set.end(), value);
    return it != m_basis_set.end() && *it == value;
  }

  std::size_t index(const BasisElement& value) const {
    auto it = std::lower_bound(m_basis_set.begin(), m_basis_set.end(), value);
    return static_cast<std::size_t>(std::distance(it, m_basis_set.begin()));
  }

  std::size_t size() const { return m_basis_set.size(); }

 protected:
  void generate_basis();
  virtual void generate_combinations(BasisElement&, size_t, size_t, size_t) = 0;

  std::size_t m_orbitals;
  std::size_t m_particles;
  std::vector<BasisElement> m_basis_set;
  NonnullOwnPtr<BasisFilter> m_basis_filter;
};

void prepare_up_and_down_representation(
    const BasisElement& element, std::vector<int>& up, std::vector<int>& down);

std::string state_string(const BasisElement& element, std::size_t orbitals);
