// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#include "Term.h"

#include <gtest/gtest.h>

using enum Operator::Type;
using enum Operator::Statistics;
using enum Operator::Spin;

TEST(TermTest, ConstructorAndAccessors) {
  std::vector<Operator> operators = {
      Operator::creation<Fermion>(Up, 0),
      Operator::annihilation<Fermion>(Down, 1)};
  Term term(2.5, operators);

  EXPECT_DOUBLE_EQ(term.coefficient().real(), 2.5);
  EXPECT_EQ(term.operators(), operators);
  EXPECT_EQ(term.operators().size(), 2);
  EXPECT_EQ(term.operators().at(0), Operator::creation<Fermion>(Up, 0));
  EXPECT_EQ(term.operators().at(1), Operator::annihilation<Fermion>(Down, 1));
}

TEST(TermTest, EqualityOperator) {
  std::vector<Operator> operators1 = {
      Operator::creation<Fermion>(Up, 0),
      Operator::annihilation<Fermion>(Down, 1)};
  Term term1(2.5, operators1);

  std::vector<Operator> operators2 = {
      Operator::creation<Fermion>(Up, 0),
      Operator::annihilation<Fermion>(Down, 1)};
  Term term2(2.5, operators2);

  std::vector<Operator> operators3 = {
      Operator::annihilation<Fermion>(Down, 1),
      Operator::creation<Fermion>(Up, 0)};
  Term term3(2.5, operators3);

  Term term4(3.0, operators1);

  EXPECT_TRUE(term1 == term2);
  EXPECT_FALSE(term1 == term3);
  EXPECT_FALSE(term1 == term4);
}

TEST(TermTest, InequalityOperator) {
  std::vector<Operator> operators1 = {
      Operator::creation<Fermion>(Up, 0),
      Operator::annihilation<Fermion>(Down, 1)};
  Term term1(2.5, operators1);

  std::vector<Operator> operators2 = {
      Operator::creation<Fermion>(Up, 0),
      Operator::annihilation<Fermion>(Down, 1)};
  Term term2(2.5, operators2);

  std::vector<Operator> operators3 = {
      Operator::annihilation<Fermion>(Down, 1),
      Operator::creation<Fermion>(Up, 0)};
  Term term3(2.5, operators3);

  Term term4(3.0, operators1);

  EXPECT_FALSE(term1 != term2);
  EXPECT_TRUE(term1 != term3);
  EXPECT_TRUE(term1 != term4);
}

TEST(TermTest, OutputOperator) {
  std::vector<Operator> operators = {
      Operator::creation<Fermion>(Up, 0),
      Operator::annihilation<Fermion>(Down, 1)};
  Term term(2.5, operators);

  std::ostringstream os;
  os << term;

  EXPECT_EQ(
      os.str(),
      "Term { Coefficient: (2.5,0), Operators: [Operator { Type: "
      "Creation, Spin: Up, Orbital: 0 }, Operator { Type: Annihilation, "
      "Spin: Down, Orbital: 1 }]}");
}

TEST(TermTest, Product) {
  std::vector<Operator> operators1 = {
      Operator::creation<Fermion>(Up, 0),
      Operator::annihilation<Fermion>(Down, 1)};
  Term term1(2.5, operators1);

  std::vector<Operator> operators2 = {
      Operator::creation<Fermion>(Down, 2),
      Operator::annihilation<Fermion>(Up, 3)};
  Term term2(3.0, operators2);

  Term product = term1.product(term2);

  std::vector<Operator> expected_operators = {
      Operator::creation<Fermion>(Up, 0),
      Operator::annihilation<Fermion>(Down, 1),
      Operator::creation<Fermion>(Down, 2),
      Operator::annihilation<Fermion>(Up, 3)};
  Term expected_product(7.5, expected_operators);

  EXPECT_EQ(product, expected_product);
}

TEST(TermTest, ProductWithVector) {
  std::vector<Operator> operators1 = {
      Operator::creation<Fermion>(Up, 0),
      Operator::annihilation<Fermion>(Down, 1)};
  Term term1(2.5, operators1);

  std::vector<Operator> operators2 = {
      Operator::creation<Fermion>(Down, 2),
      Operator::annihilation<Fermion>(Up, 3)};
  Term term2(3.0, operators2);

  Term product = term1.product(term1).product(term2);

  std::vector<Operator> expected_operators = {
      Operator::creation<Fermion>(Up, 0),
      Operator::annihilation<Fermion>(Down, 1),
      Operator::creation<Fermion>(Up, 0),
      Operator::annihilation<Fermion>(Down, 1),
      Operator::creation<Fermion>(Down, 2),
      Operator::annihilation<Fermion>(Up, 3)};
  Term expected_product(18.75, expected_operators);

  EXPECT_EQ(product, expected_product);
}

TEST(TermTest, AdjointSingleBody) {
  std::vector<Operator> operators = {
      Operator::creation<Fermion>(Up, 0),
      Operator::annihilation<Fermion>(Down, 1)};
  Term term(std::complex<double>(2.5, 1.0), operators);

  Term adjoint = term.adjoint();

  std::vector<Operator> expected_operators = {
      Operator::creation<Fermion>(Down, 1),
      Operator::annihilation<Fermion>(Up, 0)};
  Term expected_adjoint(std::complex<double>(2.5, -1.0), expected_operators);

  EXPECT_EQ(adjoint, expected_adjoint);
}

TEST(TermTest, AdjointManyBody) {
  std::vector<Operator> operators = {
      Operator::creation<Fermion>(Up, 0), Operator::creation<Fermion>(Down, 1),
      Operator::annihilation<Fermion>(Down, 2),
      Operator::annihilation<Fermion>(Up, 3)};
  Term term(std::complex<double>(2.5, 1.0), operators);

  Term adjoint = term.adjoint();

  std::vector<Operator> expected_operators = {
      Operator::creation<Fermion>(Up, 3), Operator::creation<Fermion>(Down, 2),
      Operator::annihilation<Fermion>(Down, 1),
      Operator::annihilation<Fermion>(Up, 0)};
  Term expected_adjoint(std::complex<double>(2.5, -1.0), expected_operators);

  EXPECT_EQ(adjoint, expected_adjoint);
}

TEST(TermTest, OneBodyTerm) {
  Term term = 2.5 * one_body<Fermion>(Up, 0, Down, 1);

  std::vector<Operator> operators = {
      Operator::creation<Fermion>(Up, 0),
      Operator::annihilation<Fermion>(Down, 1)};
  Term expected(2.5, operators);

  EXPECT_EQ(term, expected);
}

TEST(TermTest, ThreeBodyTermUsingNBodyTerm) {
  Term term = Term(
      2.5,
      {Operator::creation<Fermion>(Up, 0), Operator::creation<Fermion>(Down, 1),
       Operator::creation<Fermion>(Up, 2),
       Operator::annihilation<Fermion>(Down, 5),
       Operator::annihilation<Fermion>(Up, 4),
       Operator::annihilation<Fermion>(Down, 3)});

  std::vector<Operator> operators = {
      Operator::creation<Fermion>(Up, 0),
      Operator::creation<Fermion>(Down, 1),
      Operator::creation<Fermion>(Up, 2),
      Operator::annihilation<Fermion>(Down, 5),
      Operator::annihilation<Fermion>(Up, 4),
      Operator::annihilation<Fermion>(Down, 3)};
  Term expected(2.5, operators);

  EXPECT_EQ(term, expected);
}
