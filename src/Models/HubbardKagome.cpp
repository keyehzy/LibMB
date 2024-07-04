// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#include "HubbardKagome.h"

using enum Operator::Statistics;
using enum Operator::Spin;

void HubbardKagome::hopping_term(Expression& result) const {
  Expression e;
  for (Operator::Spin spin : {Up, Down}) {
    // Inner ring
    result += -m_t * hopping<Fermion>(spin, 0, 1);
    result += -m_t * hopping<Fermion>(spin, 1, 2);
    result += -m_t * hopping<Fermion>(spin, 2, 3);
    result += -m_t * hopping<Fermion>(spin, 3, 4);
    result += -m_t * hopping<Fermion>(spin, 4, 5);
    result += -m_t * hopping<Fermion>(spin, 5, 0);

    result += -m_t * hopping<Fermion>(spin, 0, 6);
    result += -m_t * hopping<Fermion>(spin, 1, 7);
    result += -m_t * hopping<Fermion>(spin, 2, 8);
    result += -m_t * hopping<Fermion>(spin, 3, 9);
    result += -m_t * hopping<Fermion>(spin, 4, 10);
    result += -m_t * hopping<Fermion>(spin, 5, 11);

    result += -m_t * hopping<Fermion>(spin, 1, 6);
    result += -m_t * hopping<Fermion>(spin, 2, 7);
    result += -m_t * hopping<Fermion>(spin, 3, 8);
    result += -m_t * hopping<Fermion>(spin, 4, 9);
    result += -m_t * hopping<Fermion>(spin, 5, 10);
    result += -m_t * hopping<Fermion>(spin, 0, 11);

    // PBC
    if (m_periodic) {
      result += -m_t * hopping<Fermion>(spin, 6, 10);
      result += -m_t * hopping<Fermion>(spin, 7, 11);
      result += -m_t * hopping<Fermion>(spin, 8, 6);
      result += -m_t * hopping<Fermion>(spin, 9, 7);
      result += -m_t * hopping<Fermion>(spin, 10, 8);
      result += -m_t * hopping<Fermion>(spin, 11, 9);
    }
  }
}

void HubbardKagome::interaction_term(Expression& result) const {
  for (size_t i1 = 0; i1 < size; i1++) {
    result += m_u * density_density<Fermion>(Up, i1, Down, i1);
  }
}
