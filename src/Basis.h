// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <algorithm>
#include <boost/container/flat_map.hpp>
#include <unordered_map>
#include <vector>

#include "BasisFilter.h"
#include "Operator.h"
#include "Pointers/NonnullOwnPtr.h"

using BasisElement = std::vector<Operator>;

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
           m_basis_map == other.m_basis_map;
  }

  bool operator!=(const Basis& other) const { return !(*this == other); }

  const boost::container::flat_map<BasisElement, std::size_t>& elements()
      const {
    return m_basis_map;
  }

  const BasisElement& element(std::size_t i) const {
    auto it = m_basis_map.begin();
    std::advance(it, i);
    return it->first;
  }

  void insert(const BasisElement& value) {
    m_basis_map[value] = m_basis_map.size();
  }

  bool contains(const BasisElement& term) const {
    return m_basis_map.find(term) != m_basis_map.end();
  }

  std::size_t index(const BasisElement& term) const {
    return m_basis_map.at(term);
  }

  std::size_t size() const { return m_basis_map.size(); }

  std::vector<BasisElement> elements_as_vector() const {
    std::vector<BasisElement> v;
    for (const auto& [element, index] : m_basis_map) {
      v.push_back(element);
    }
    return v;
  }

 protected:
  void generate_basis();
  virtual void generate_combinations(BasisElement&, size_t, size_t, size_t) = 0;

  std::size_t m_orbitals;
  std::size_t m_particles;
  boost::container::flat_map<BasisElement, std::size_t> m_basis_map;
  NonnullOwnPtr<BasisFilter> m_basis_filter;
};

void prepare_up_and_down_representation(
    const BasisElement& element, std::vector<int>& up, std::vector<int>& down);

std::string state_string(const BasisElement& element, std::size_t orbitals);
