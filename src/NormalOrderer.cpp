// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#include "NormalOrderer.h"

#include <deque>
#include <vector>

constexpr Term::CoeffType evaluate_parity(
    Term::CoeffType coefficient, std::size_t phase) {
  return phase % 2 == 0 ? coefficient : -coefficient;
}

NormalOrderer::NormalOrderer(const Term& term) {
  normal_order(term.operators(), term.coefficient());
}

NormalOrderer::NormalOrderer(const std::vector<Term>& terms) {
  for (const Term& term : terms) {
    normal_order(term.operators(), term.coefficient());
  }
}

NormalOrderer::NormalOrderer(const Expression& expression) {
  for (const auto& [operators, coeff] : expression.terms()) {
    normal_order(operators, coeff);
  }
}

NormalOrderer::NormalOrderer(const std::vector<Expression>& expressions) {
  for (const Expression& expression : expressions) {
    for (const auto& [operators, coeff] : expression.terms()) {
      normal_order(operators, coeff);
    }
  }
}

void NormalOrderer::normal_order(
    const std::vector<Operator>& operators, Term::CoeffType coefficient) {
  std::deque<OperatorsPhasePair> queue;
  queue.emplace_back(operators, 0);
  while (!queue.empty()) {
    auto [prev_operators, prev_phase] = std::move(queue.back());
    queue.pop_back();

    if (prev_operators.size() < 2) {
      m_terms_map[prev_operators] += evaluate_parity(coefficient, prev_phase);
      continue;
    }

    auto [new_operators, new_phase] =
        sort_operators(prev_operators, prev_phase, queue);
    m_terms_map[new_operators] += evaluate_parity(coefficient, new_phase);
  }
}

NormalOrderer::OperatorsPhasePair NormalOrderer::sort_operators(
    std::vector<Operator> operators, std::size_t phase,
    std::deque<OperatorsPhasePair>& queue) {
  for (std::size_t i = 1; i < operators.size(); ++i) {
    for (std::size_t j = i; j > 0; --j) {
      Operator& op1 = operators[j - 1];
      Operator& op2 = operators[j];
      if (op1.type() == Operator::Type::Creation &&
          op2.type() == Operator::Type::Creation &&
          op1.identifier() > op2.identifier()) {
        std::swap(op1, op2);
        phase += op1.is_fermion() && op2.is_fermion();
      } else if (
          op1.type() == Operator::Type::Annihilation &&
          op2.type() == Operator::Type::Annihilation &&
          op1.identifier() < op2.identifier()) {
        std::swap(op1, op2);
        phase += op1.is_fermion() && op2.is_fermion();
      } else if (
          op1.type() == Operator::Type::Annihilation &&
          op2.type() == Operator::Type::Creation) {
        if (op1.identifier() == op2.identifier()) {
          std::vector<Operator> elements(operators);
          elements.erase(elements.begin() + j - 1, elements.begin() + j + 1);
          queue.emplace_back(std::move(elements), phase);
        }
        std::swap(op1, op2);
        phase += op1.is_fermion() && op2.is_fermion();
      }
    }
  }
  return OperatorsPhasePair{operators, phase};
}

Expression commute(const Term& term1, const Term& term2) {
  return NormalOrderer({term1.product(term2), term2.product(term1).negate()})
      .expression();
}

Expression commute(
    const Expression& expression1, const Expression& expression2) {
  return NormalOrderer({expression1.product(expression2),
                        expression2.product(expression1).negate()})
      .expression();
}

Expression anticommute(const Term& term1, const Term& term2) {
  return NormalOrderer({term1.product(term2), term2.product(term1)})
      .expression();
}

Expression anticommute(
    const Expression& expression1, const Expression& expression2) {
  return NormalOrderer({expression1.product(expression2),
                        expression2.product(expression1)})
      .expression();
  ;
}
