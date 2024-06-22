// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#include "HubbardChain.h"

void HubbardChain::hopping_term(Expression& result) const {
  for (Operator::Spin spin : {Up, Down}) {
    for (std::size_t i = 0; i < m_size; i++) {
      result += density<Fermion>(-m_mu, spin, i);
      result += hopping<Fermion>(-m_t, spin, i, (i + 1) % m_size);
    }
  }
}

void HubbardChain::interaction_term(Expression& result) const {
  for (size_t i = 0; i < m_size; i++) {
    result += density_density<Fermion>(m_u, Up, i, Down, i);
  }
}
