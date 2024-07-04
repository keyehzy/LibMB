// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <algorithm>
#include <unordered_map>
#include <vector>

#include "BasisFilter.h"
#include "IndexedVectorMap.h"
#include "Operator.h"
#include "Pointers/NonnullOwnPtr.h"

using BasisElement = std::vector<Operator>;
using BasisMap = std::unordered_map<BasisElement, std::size_t>;

class Basis {
 public:
  virtual ~Basis() = default;

  Basis(std::size_t n, std::size_t m)
      : m_orbitals{n}, m_particles{m}, m_basis_filter{make<BasisFilter>()} {}

  Basis(std::size_t n, std::size_t m, BasisFilter* filter)
      : m_orbitals{n}, m_particles{m}, m_basis_filter{adopt_own(filter)} {}

  const std::vector<BasisElement>& elements() const {
    return m_basis_map.elements();
  }

  const BasisElement& element(std::size_t i) const { return m_basis_map[i]; }

  std::size_t orbitals() const { return m_orbitals; }

  std::size_t particles() const { return m_particles; }

  bool operator==(const Basis& other) const {
    return m_orbitals == other.m_orbitals && m_particles == other.m_particles &&
           m_basis_map.elements() == other.m_basis_map.elements();
  }

  bool operator!=(const Basis& other) const { return !(*this == other); }

  bool contains(const BasisElement& term) const {
    return m_basis_map.contains(term);
  }

  std::size_t index(const BasisElement& term) const {
    return m_basis_map.index(term);
  }

  std::size_t size() const { return m_basis_map.size(); }

  template <typename CompareFunction>
  void sort(CompareFunction comp) {
    m_basis_map.sort(comp);
  }

 protected:
  void generate_basis();
  virtual void generate_combinations(BasisElement&, size_t, size_t, size_t) = 0;

  std::size_t m_orbitals;
  std::size_t m_particles;
  IndexedVectorMap<BasisElement> m_basis_map;
  NonnullOwnPtr<BasisFilter> m_basis_filter;
};

void prepare_up_and_down_representation(
    const BasisElement& element, std::vector<int>& up, std::vector<int>& down);

std::string state_string(const BasisElement& element, std::size_t orbitals);
