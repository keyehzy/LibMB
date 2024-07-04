// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#include "LinearChain.h"

static constexpr auto Fermion = Operator::Statistics::Fermion;

Expression LinearChain::hamiltonian() const {
  Expression result;
  for (Operator::Spin spin : {Operator::Spin::Up, Operator::Spin::Down}) {
    for (std::size_t i = 0; i < m_size - 1; i++) {
      result += -m_t * one_body<Fermion>(spin, i, spin, i + 1);
      result += -m_t * one_body<Fermion>(spin, i, spin, i + 1).adjoint();
    }
    result += -m_t * one_body<Fermion>(spin, m_size - 1, spin, 0);
    result += -m_t * one_body<Fermion>(spin, m_size - 1, spin, 0).adjoint();

    for (std::size_t i = 0; i < m_size; i++) {
      result += -m_u * one_body<Fermion>(spin, i, spin, i);
    }
  }
  return result;
}
