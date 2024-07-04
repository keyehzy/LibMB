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

  Expression operator-() const { return (*this).negate(); }

  Expression operator+(const Expression& rhs) const { return add(rhs); }

  Expression operator-(const Expression& rhs) const {
    return add(rhs.negate());
  }

  Expression operator*(const Expression& rhs) const { return product(rhs); }

  Expression operator*(const std::vector<Operator>& rhs) const {
    return product(rhs);
  }

  Expression operator*(std::complex<double> coefficient) const {
    return product(coefficient);
  }

  Expression adjoint() const;

  friend Expression operator*(
      std::complex<double> coefficient, const Expression& other) {
    return other * coefficient;
  }

  friend std::ostream& operator<<(std::ostream& os, const Expression& e);

 private:
  Expression add(const Expression& other) const;
  Expression add(const Term& other) const;

  Expression product(const Expression& other) const;
  Expression product(const Term& other) const;
  Expression product(const std::vector<Operator>& other) const;
  Expression product(std::complex<double> coefficient) const;

  Expression negate() const;

  ExpressionMap m_terms;
};

Expression operator+(const Term& a, const Term& b);

Expression operator-(const Term& a, const Term& b);

template <Operator::Statistics S>
Expression hopping(Operator::Spin spin, std::size_t i, std::size_t j) {
  return one_body<S>(spin, i, spin, j) +
         one_body<S>(spin, i, spin, j).adjoint();
}

// NOTE: These are only spin-1/2 operators. We implement them here in terms of
// spin flippings of Fermions. However, for spin-1 particles we would have to
// implement them in terms of Bosons.
Expression spin_x(std::size_t i);

Expression spin_y(std::size_t i);

Expression spin_z(std::size_t i);
