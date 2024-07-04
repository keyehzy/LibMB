// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#include "HubbardChain.h"

void HubbardChain::hopping_term(Expression& result) const {
  for (Operator::Spin spin : {Up, Down}) {
    for (std::size_t i = 0; i < m_size; i++) {
      result += -m_mu * density<Fermion>(spin, i);
      result += -m_t * hopping<Fermion>(spin, i, (i + 1) % m_size);
    }
  }
}

void HubbardChain::interaction_term(Expression& result) const {
  for (size_t i = 0; i < m_size; i++) {
    result += m_u * density_density<Fermion>(Up, i, Down, i);
  }
}
