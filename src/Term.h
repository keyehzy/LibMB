// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <algorithm>
#include <complex>
#include <vector>

#include "Operator.h"

class Term {
 public:
  using CoeffType = std::complex<double>;

  Term(CoeffType coefficient, const std::vector<Operator>& operators)
      : m_coefficient{coefficient}, m_operators{operators} {}

  Term() = default;

  CoeffType coefficient() const { return m_coefficient; }

  const std::vector<Operator>& operators() const { return m_operators; }

  std::vector<Operator>& operators() { return m_operators; }

  bool operator==(const Term& other) const {
    return m_coefficient == other.m_coefficient &&
           m_operators == other.m_operators;
  }

  bool operator!=(const Term& other) const { return !(*this == other); }

  Term& operator*=(const Term& other) {
    m_coefficient *= other.coefficient();
    m_operators.insert(
        m_operators.end(), other.operators().begin(), other.operators().end());
    return *this;
  }

  Term& operator*=(const std::vector<Operator>& other) {
    m_operators.insert(m_operators.end(), other.begin(), other.end());
    return *this;
  }

  Term& operator*=(double other) {
    m_coefficient *= other;
    return *this;
  }

  Term product(const Term& other) const {
    std::vector<Operator> new_operators = m_operators;
    new_operators.insert(
        new_operators.end(), other.m_operators.begin(),
        other.m_operators.end());
    return Term(m_coefficient * other.m_coefficient, new_operators);
  }

  Term product(const std::vector<Operator>& operators) const {
    std::vector<Operator> new_operators = m_operators;
    new_operators.insert(
        new_operators.end(), operators.begin(), operators.end());
    return Term(m_coefficient, new_operators);
  }

  friend Term operator*(const Term& lhs, const Term& rhs) {
    return lhs.product(rhs);
  }

  friend Term operator*(const Term& lhs, const std::vector<Operator>& rhs) {
    return lhs.product(rhs);
  }

  friend std::ostream& operator<<(std::ostream& os, const Term& term);

  Term adjoint() const {
    std::vector<Operator> adj_operators;
    for (const auto& op : m_operators) {
      adj_operators.push_back(op.adjoint());
    }
    std::reverse(adj_operators.begin(), adj_operators.end());
    return Term(std::conj(m_coefficient), adj_operators);
  }

  Term negate() const { return Term(-m_coefficient, m_operators); }

 private:
  CoeffType m_coefficient;
  std::vector<Operator> m_operators;
};

template <Operator::Statistics S>
Term one_body(
    Term::CoeffType coefficient, Operator::Spin spin1, std::size_t orbital1,
    Operator::Spin spin2, std::size_t orbital2) {
  return Term(
      coefficient, {Operator::creation<S>(spin1, orbital1),
                    Operator::annihilation<S>(spin2, orbital2)});
}

template <Operator::Statistics S>
Term density(
    Term::CoeffType coefficient, Operator::Spin spin, std::size_t orbital) {
  return Term(
      coefficient, {Operator::creation<S>(spin, orbital),
                    Operator::annihilation<S>(spin, orbital)});
}

template <Operator::Statistics S>
Term two_body(
    Term::CoeffType coefficient, Operator::Spin spin1, std::size_t orbital1,
    Operator::Spin spin2, std::size_t orbital2, Operator::Spin spin3,
    std::size_t orbital3, Operator::Spin spin4, std::size_t orbital4) {
  return Term(
      coefficient, {Operator::creation<S>(spin1, orbital1),
                    Operator::annihilation<S>(spin2, orbital2),
                    Operator::creation<S>(spin3, orbital3),
                    Operator::annihilation<S>(spin4, orbital4)});
}

template <Operator::Statistics S>
Term density_density(
    Term::CoeffType coefficient, Operator::Spin spin1, std::size_t orbital1,
    Operator::Spin spin2, std::size_t orbital2) {
  return Term(
      coefficient, {Operator::creation<S>(spin1, orbital1),
                    Operator::annihilation<S>(spin1, orbital1),
                    Operator::creation<S>(spin2, orbital2),
                    Operator::annihilation<S>(spin2, orbital2)});
}
