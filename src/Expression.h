// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <unordered_map>
#include <vector>

#include "Operator.h"
#include "Term.h"

class Expression {
 public:
  using ExpressionMap =
      std::unordered_map<std::vector<Operator>, Term::CoeffType>;

  Expression() = default;

  Expression(const ExpressionMap& terms) : m_terms(terms) {}

  Expression(ExpressionMap&& terms) : m_terms(std::move(terms)) {}

  Expression(const std::vector<Term>& terms) {
    for (const auto& term : terms) {
      m_terms[term.operators()] += term.coefficient();
    }
  }

  void insert(const Term& term) {
    m_terms[term.operators()] += term.coefficient();
  }

  void insert(const Expression& other) {
    for (const auto& [operators, coefficient] : other.terms()) {
      m_terms[operators] += coefficient;
    }
  }

  void insert(double coefficient) { m_terms[{}] += coefficient; }

  Expression& operator+=(const Expression& other) {
    insert(other);
    return *this;
  }

  Expression& operator+=(const Term& other) {
    insert(other);
    return *this;
  }

  Expression& operator+=(double other) {
    insert(other);
    return *this;
  }

  std::size_t size() const { return m_terms.size(); }

  const ExpressionMap& terms() const { return m_terms; }

  ExpressionMap& terms() { return m_terms; }

  bool operator==(const Expression& other) const {
    return m_terms == other.m_terms;
  }

  bool operator!=(const Expression& other) const { return !(*this == other); }

  Expression add(const Expression& other) const {
    Expression result(*this);
    for (const auto& [operators, coefficient] : other.terms()) {
      result.insert(Term(coefficient, operators));
    }
    return result;
  }

  Expression add(const Term& other) const {
    Expression result(*this);
    result.insert(other);
    return result;
  }

  friend Expression operator+(const Expression& lhs, const Expression& rhs) {
    return lhs.add(rhs);
  }

  friend Expression operator+(const Expression& lhs, const Term& rhs) {
    return lhs.add(rhs);
  }

  Expression product(const Expression& other) const {
    Expression result;
    for (const auto& [operators_a, coefficient_a] : terms()) {
      for (const auto& [operators_b, coefficient_b] : other.terms()) {
        result.insert(Term(coefficient_a, operators_a)
                          .product(Term(coefficient_b, operators_b)));
      }
    }
    return result;
  }

  Expression product(const Term& other) const {
    Expression result;
    for (const auto& [operators_a, coefficient_a] : terms()) {
      result.insert(Term(coefficient_a, operators_a).product(other));
    }
    return result;
  }

  Expression product(const std::vector<Operator>& other) const {
    Expression result;
    for (const auto& [operators_a, coefficient_a] : terms()) {
      result.insert(Term(coefficient_a, operators_a).product(other));
    }
    return result;
  }

  Expression product(double coefficient) const {
    Expression result;
    for (const auto& [operators_a, coefficient_a] : terms()) {
      result.insert(Term(coefficient * coefficient_a, operators_a));
    }
    return result;
  }

  friend Expression operator*(const Expression& lhs, const Expression& rhs) {
    return lhs.product(rhs);
  }

  friend Expression operator*(const Expression& lhs, const Term& rhs) {
    return lhs.product(rhs);
  }

  friend Expression operator*(
      const Expression& lhs, const std::vector<Operator>& rhs) {
    return lhs.product(rhs);
  }

  friend Expression operator*(double coefficient, const Expression& other) {
    return other.product(coefficient);
  }

  Expression adjoint() const {
    Expression result;
    for (const auto& [operators, coefficient] : terms()) {
      result.insert(Term(coefficient, operators).adjoint());
    }
    return result;
  }

  Expression negate() const {
    Expression result;
    for (const auto& [operators, coefficient] : terms()) {
      result.insert(Term(-coefficient, operators));
    }
    return result;
  }

  friend std::ostream& operator<<(std::ostream& os, const Expression& e);

 private:
  ExpressionMap m_terms;
};

Expression add(const Term& a, const Term& b);

template <Operator::Statistics S>
Expression hopping(
    double t, Operator::Spin spin, std::size_t i, std::size_t j) {
  return add(
      t * one_body<S>(spin, i, spin, j),
      t * one_body<S>(spin, i, spin, j).adjoint());
}

// NOTE: These are only spin-1/2 operators. We implement them here in terms of
// spin flippings of Fermions. However, for spin-1 particles we would have to
// implement them in terms of Bosons.
Expression spin_x(std::size_t i);

Expression spin_y(std::size_t i);

Expression spin_z(std::size_t i);
