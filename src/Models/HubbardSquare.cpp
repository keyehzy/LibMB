// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#include "HubbardSquare.h"

#include "Operator.h"

static constexpr auto Fermion = Operator::Statistics::Fermion;

void HubbardSquare::hopping_term(Expression& result) const {
  auto index = [&](size_t i, size_t j) { return j * m_nx + i; };
  for (size_t i = 0; i < m_nx; i++) {
    for (size_t j = 0; j < m_ny; j++) {
      for (Operator::Spin spin : {Operator::Spin::Up, Operator::Spin::Down}) {
        size_t dsti = i < m_nx - 1 ? i + 1 : 0;
        size_t dstj = j < m_ny - 1 ? j + 1 : 0;
        result +=
            -m_t * one_body<Fermion>(spin, index(i, j), spin, index(dsti, j));
        result +=
            -m_t * one_body<Fermion>(spin, index(dsti, j), spin, index(i, j));
        result +=
            -m_t * one_body<Fermion>(spin, index(i, j), spin, index(i, dstj));
        result +=
            -m_t * one_body<Fermion>(spin, index(i, dstj), spin, index(i, j));
      }
    }
  }
}

void HubbardSquare::interaction_term(Expression& result) const {
  for (size_t i1 = 0; i1 < m_nx * m_ny; i1++) {
    result += density_density<Fermion>(
        m_u, Operator::Spin::Up, i1, Operator::Spin::Down, i1);
  }
}
