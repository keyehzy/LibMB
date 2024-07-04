// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#include "Expression.h"

#include <iostream>
#include <ostream>

using enum Operator::Statistics;
using enum Operator::Spin;
static constexpr auto ju = std::complex<double>{0, 1};

Expression Expression::add(const Expression& other) const {
  Expression result(*this);
  for (const auto& [operators, coefficient] : other.terms()) {
    result.insert(Term(coefficient, operators));
  }
  return result;
}

Expression Expression::add(const Term& other) const {
  Expression result(*this);
  result.insert(other);
  return result;
}

Expression Expression::product(const Expression& other) const {
  Expression result;
  for (const auto& [operators_a, coefficient_a] : terms()) {
    for (const auto& [operators_b, coefficient_b] : other.terms()) {
      result.insert(
          Term(coefficient_a, operators_a) * Term(coefficient_b, operators_b));
    }
  }
  return result;
}

Expression Expression::product(const Term& other) const {
  Expression result;
  for (const auto& [operators_a, coefficient_a] : terms()) {
    result.insert(Term(coefficient_a, operators_a) * other);
  }
  return result;
}

Expression Expression::product(const std::vector<Operator>& other) const {
  Expression result;
  for (const auto& [operators_a, coefficient_a] : terms()) {
    result.insert(Term(coefficient_a, operators_a) * other);
  }
  return result;
}

Expression Expression::product(std::complex<double> coefficient) const {
  Expression result;
  for (auto& [operators_a, coefficient_a] : terms()) {
    result.insert(Term(coefficient * coefficient_a, operators_a));
  }
  return result;
}

Expression Expression::negate() const {
  Expression result;
  for (const auto& [operators, coefficient] : terms()) {
    result.insert(Term(-coefficient, operators));
  }
  return result;
}

Expression Expression::adjoint() const {
  Expression result;
  for (const auto& [operators, coefficient] : terms()) {
    result.insert(Term(coefficient, operators).adjoint());
  }
  return result;
}

std::ostream& operator<<(std::ostream& os, const Expression& e) {
  for (const auto& [operators, coeff] : e.terms()) {
    os << coeff << "  {";
    for (Operator o : operators) {
      os << o << ", ";
    }
    os << "}\n";
  }
  return os;
}

Expression operator+(const Term& a, const Term& b) {
  Expression result;
  result.insert(a);
  result.insert(b);
  return result;
}

Expression operator-(const Term& a, const Term& b) {
  Expression result;
  result.insert(a);
  result.insert(-b);
  return result;
}

Expression spin_x(std::size_t i) {
  return 0.5 * (spin_flip<Fermion>(i) + spin_flip<Fermion>(i).adjoint());
}

Expression spin_y(std::size_t i) {
  return 0.5 * ju * (-spin_flip<Fermion>(i) + spin_flip<Fermion>(i).adjoint());
}

Expression spin_z(std::size_t i) {
  return 0.5 * (density<Fermion>(Up, i) - density<Fermion>(Down, i));
}
