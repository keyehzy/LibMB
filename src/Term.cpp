// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#include "Term.h"

Term Term::product(const Term& other) const {
  std::vector<Operator> new_operators(m_operators);
  new_operators.insert(
      new_operators.end(), other.m_operators.begin(), other.m_operators.end());
  return Term(m_coefficient * other.m_coefficient, new_operators);
}

Term Term::product(const std::vector<Operator>& operators) const {
  std::vector<Operator> new_operators(m_operators);
  new_operators.insert(new_operators.end(), operators.begin(), operators.end());
  return Term(m_coefficient, new_operators);
}

Term Term::product(std::complex<double> other) const {
  return Term(other * m_coefficient, m_operators);
}

Term Term::adjoint() const {
  std::vector<Operator> adj_operators;
  adj_operators.reserve(adj_operators.size());
  for (auto it = m_operators.rbegin(); it != m_operators.rend(); it++) {
    adj_operators.push_back(it->adjoint());
  }
  return Term(std::conj(m_coefficient), adj_operators);
}

std::ostream& operator<<(std::ostream& os, const Term& term) {
  os << "Term { Coefficient: " << term.coefficient();
  os << ", Operators: [";
  for (size_t i = 0; i < term.operators().size() - 1; ++i) {
    os << term.operators()[i] << ", ";
  }
  os << term.operators()[term.operators().size() - 1];
  os << "]}";
  return os;
}
