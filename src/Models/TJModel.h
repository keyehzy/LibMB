// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include "Model.h"

using enum Operator::Statistics;
using enum Operator::Spin;
using enum Operator::Statistics;

class TJModel : public Model {
 public:
  TJModel(double t, double J, std::size_t size, bool periodic = true)
      : m_t(t), m_J(J), m_size(size), m_periodic(periodic) {}

  ~TJModel() override {}

  std::size_t size() const { return m_size; }

 private:
  Expression hamiltonian() const override {
    Expression result;
    hopping_term(result);
    exchange_term(result);
    return result;
  }

  void hopping_term(Expression& result) const {
    for (Operator::Spin spin : {Operator::Spin::Up, Operator::Spin::Down}) {
      for (std::size_t i = 0; i < m_size - 1; ++i) {
        result += -m_t * hopping<Fermion>(spin, i, i + 1);
      }
      if (m_periodic) {
        result += -m_t * hopping<Fermion>(spin, m_size - 1, 0);
      }
    }
  }

  void exchange_term(Expression& result) const {
    for (std::size_t i = 0; i < m_size - 1; ++i) {
      result += exchange_interaction(i, i + 1);
    }
    if (m_periodic) {
      result += exchange_interaction(m_size - 1, 0);
    }
  }

  Expression exchange_interaction(std::size_t i, std::size_t j) const {
    Expression result;
    result += m_J * spin_dot(i, j);
    result +=
        -0.25 * m_J * total_density<Fermion>(i) * total_density<Fermion>(j);
    return result;
  }

  double m_t;
  double m_J;
  std::size_t m_size;
  bool m_periodic;
};
