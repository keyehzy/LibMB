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
  void hopping_term(Expression& result) const;

  void interaction_term(Expression& result) const;

  Expression hamiltonian() const override {
    Expression result;
    hopping_term(result);
    interaction_term(result);
    return result;
  }

  double m_mu;
  double m_t;
  double m_u;
  size_t m_size;
};
