// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <algorithm>
#include <cstddef>
#include <unordered_map>
#include <vector>

template <class T>
class IndexedVectorMap {
 public:
  using element_type = T;
  using map_type = ;

 private:
  std::size_t m_size = 0;
  map_type m_index_map;

 public:
  const map_type& index_map() const { return m_index_map; }

  void insert(const T& value) {
    m_index_map[value] = m_size;
    m_size += 1;
  }

  const T& operator[](std::size_t idx) const {
    auto it = m_index_map.begin();
    std::advance(it, idx);
    return it->first;
  }

  std::size_t index(const T& value) const { return m_index_map.at(value); }

  bool contains(const T& value) const {
    return m_index_map.find(value) != m_index_map.end();
  }

  std::size_t size() const { return m_size; }

  bool operator==(const IndexedVectorMap& other) const {
    return m_size == other.m_size && m_index_map == other.m_index_map;
  }
};
