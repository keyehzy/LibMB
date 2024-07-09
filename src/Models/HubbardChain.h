// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include "Model.h"

using enum Operator::Statistics;
using enum Operator::Spin;

class HubbardChain : public Model {
 public:
  HubbardChain(double mu, double t, double u, size_t n)
      : m_mu(mu), m_t(t), m_u(u), m_size(n) {}

  ~HubbardChain() override {}

 private:
  Expression hamiltonian() const override {
    Expression result;

    // Chemical potential and Hopping
    for (Operator::Spin spin : {Up, Down}) {
      for (std::size_t i = 0; i < m_size; i++) {
        result += -m_mu * density<Fermion>(spin, i);
        result += -m_t * hopping<Fermion>(spin, i, (i + 1) % m_size);
      }
    }

    // Hubbard U
    for (size_t i = 0; i < m_size; i++) {
      result += m_u * density_density<Fermion>(Up, i, Down, i);
    }

    return result;
  }

  double m_mu;
  double m_t;
  double m_u;
  size_t m_size;
};
