// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#include "NormalOrder.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

using testing::IsEmpty;

TEST(NormalOrderTest, NormalOrderTermEqual) {
  {
    Term term(
        1.0,
        {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                  Operator::Spin::UP, 0),
         Operator(Operator::Type::ANNIHILATION, Operator::Statistics::FERMION,
                  Operator::Spin::DOWN, 1)});
    Expression normal_ordered = normal_order(term);

    std::vector<Term> terms = {Term(
        1.0,
        {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                  Operator::Spin::UP, 0),
         Operator(Operator::Type::ANNIHILATION, Operator::Statistics::FERMION,
                  Operator::Spin::DOWN, 1)})};
    Expression expected(terms);

    EXPECT_EQ(normal_ordered, expected);
  }

  {
    Term term(1.0,
              {Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                        Operator::Spin::UP, 0),
               Operator(Operator::Type::ANNIHILATION,
                        Operator::Statistics::BOSON, Operator::Spin::DOWN, 1)});
    Expression normal_ordered = normal_order(term);

    std::vector<Term> terms = {Term(
        1.0, {Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::BOSON, Operator::Spin::DOWN, 1)})};
    Expression expected(terms);

    EXPECT_EQ(normal_ordered, expected);
  }
}

TEST(NormalOrderTest, NormalOrderTermCreationCreation) {
  {
    Term term(1.0,
              {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                        Operator::Spin::UP, 1),
               Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                        Operator::Spin::UP, 0)});
    Expression normal_ordered = normal_order(term);

    std::vector<Term> terms = {Term(
        -1.0, {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                        Operator::Spin::UP, 0),
               Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                        Operator::Spin::UP, 1)})};
    Expression expected(terms);

    EXPECT_EQ(normal_ordered, expected);
  }

  {
    Term term(1.0,
              {Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                        Operator::Spin::UP, 1),
               Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                        Operator::Spin::UP, 0)});
    Expression normal_ordered = normal_order(term);

    std::vector<Term> terms = {Term(
        1.0, {Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 1)})};
    Expression expected(terms);

    EXPECT_EQ(normal_ordered, expected);
  }
}

TEST(NormalOrderTest, NormalOrderTermAnnihilationAnnihilation) {
  {
    Term term(
        1.0,
        {Operator(Operator::Type::ANNIHILATION, Operator::Statistics::FERMION,
                  Operator::Spin::DOWN, 0),
         Operator(Operator::Type::ANNIHILATION, Operator::Statistics::FERMION,
                  Operator::Spin::DOWN, 1)});
    Expression normal_ordered = normal_order(term);

    std::vector<Term> terms = {Term(
        -1.0,
        {Operator(Operator::Type::ANNIHILATION, Operator::Statistics::FERMION,
                  Operator::Spin::DOWN, 1),
         Operator(Operator::Type::ANNIHILATION, Operator::Statistics::FERMION,
                  Operator::Spin::DOWN, 0)})};
    Expression expected(terms);

    EXPECT_EQ(normal_ordered, expected);
  }

  {
    Term term(1.0,
              {Operator(Operator::Type::ANNIHILATION,
                        Operator::Statistics::BOSON, Operator::Spin::DOWN, 0),
               Operator(Operator::Type::ANNIHILATION,
                        Operator::Statistics::BOSON, Operator::Spin::DOWN, 1)});
    Expression normal_ordered = normal_order(term);

    std::vector<Term> terms = {Term(
        1.0, {Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::BOSON, Operator::Spin::DOWN, 1),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::BOSON, Operator::Spin::DOWN, 0)})};
    Expression expected(terms);

    EXPECT_EQ(normal_ordered, expected);
  }
}

TEST(NormalOrderTest, NormalOrderTermCreationAnnihilation) {
  {
    Term term(
        1.0,
        {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                  Operator::Spin::UP, 0),
         Operator(Operator::Type::ANNIHILATION, Operator::Statistics::FERMION,
                  Operator::Spin::DOWN, 0)});
    Expression normal_ordered = normal_order(term);

    std::vector<Term> terms = {Term(
        1.0,
        {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                  Operator::Spin::UP, 0),
         Operator(Operator::Type::ANNIHILATION, Operator::Statistics::FERMION,
                  Operator::Spin::DOWN, 0)})};
    Expression expected(terms);

    EXPECT_EQ(normal_ordered, expected);
  }

  {
    Term term(1.0,
              {Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                        Operator::Spin::UP, 0),
               Operator(Operator::Type::ANNIHILATION,
                        Operator::Statistics::BOSON, Operator::Spin::DOWN, 0)});
    Expression normal_ordered = normal_order(term);

    std::vector<Term> terms = {Term(
        1.0, {Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::BOSON, Operator::Spin::DOWN, 0)})};
    Expression expected(terms);

    EXPECT_EQ(normal_ordered, expected);
  }
}

TEST(NormalOrderTest, NormalOrderTermAnnihilationCreationDifferentSpin) {
  {
    Term term(1.0,
              {Operator(Operator::Type::ANNIHILATION,
                        Operator::Statistics::FERMION, Operator::Spin::DOWN, 0),
               Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                        Operator::Spin::UP, 0)});
    Expression normal_ordered = normal_order(term);

    std::vector<Term> terms = {Term(
        -1.0,
        {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                  Operator::Spin::UP, 0),
         Operator(Operator::Type::ANNIHILATION, Operator::Statistics::FERMION,
                  Operator::Spin::DOWN, 0)})};
    Expression expected(terms);

    EXPECT_EQ(normal_ordered, expected);
  }

  {
    Term term(1.0,
              {Operator(Operator::Type::ANNIHILATION,
                        Operator::Statistics::BOSON, Operator::Spin::DOWN, 0),
               Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                        Operator::Spin::UP, 0)});
    Expression normal_ordered = normal_order(term);

    std::vector<Term> terms = {Term(
        1.0, {Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::BOSON, Operator::Spin::DOWN, 0)})};
    Expression expected(terms);

    EXPECT_EQ(normal_ordered, expected);
  }
}

TEST(NormalOrderTest, NormalOrderTermAnnihilationCreationDifferentOrbital) {
  {
    Term term(1.0,
              {Operator(Operator::Type::ANNIHILATION,
                        Operator::Statistics::FERMION, Operator::Spin::UP, 1),
               Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                        Operator::Spin::UP, 0)});
    Expression normal_ordered = normal_order(term);

    std::vector<Term> terms = {
        Term(-1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::FERMION, Operator::Spin::UP, 1)})};
    Expression expected(terms);

    EXPECT_EQ(normal_ordered, expected);
  }

  {
    Term term(1.0,
              {Operator(Operator::Type::ANNIHILATION,
                        Operator::Statistics::BOSON, Operator::Spin::UP, 1),
               Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                        Operator::Spin::UP, 0)});
    Expression normal_ordered = normal_order(term);

    std::vector<Term> terms = {Term(
        1.0, {Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::BOSON, Operator::Spin::UP, 1)})};
    Expression expected(terms);

    EXPECT_EQ(normal_ordered, expected);
  }
}

TEST(NormalOrderTest, NormalOrderTermAnnihilationCreationSameSpinSameOrbital) {
  {
    Term term(1.0,
              {Operator(Operator::Type::ANNIHILATION,
                        Operator::Statistics::FERMION, Operator::Spin::UP, 0),
               Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                        Operator::Spin::UP, 0)});
    Expression normal_ordered = normal_order(term);

    std::vector<Term> terms = {
        Term(-1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::FERMION, Operator::Spin::UP, 0)}),
        Term(1.0, {})};
    Expression expected(terms);

    EXPECT_EQ(normal_ordered, expected);
  }

  {
    Term term(1.0,
              {Operator(Operator::Type::ANNIHILATION,
                        Operator::Statistics::BOSON, Operator::Spin::UP, 0),
               Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                        Operator::Spin::UP, 0)});
    Expression normal_ordered = normal_order(term);

    std::vector<Term> terms = {
        Term(1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::BOSON, Operator::Spin::UP, 0)}),
        Term(1.0, {})};
    Expression expected(terms);

    EXPECT_EQ(normal_ordered, expected);
  }
}

TEST(NormalOrderTest,
     NormalOrderTermCreationCreationAnnihilationDifferentOrbital) {
  {
    Term term(1.0,
              {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                        Operator::Spin::UP, 1),
               Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                        Operator::Spin::UP, 0),
               Operator(Operator::Type::ANNIHILATION,
                        Operator::Statistics::FERMION, Operator::Spin::UP, 0)});
    Expression normal_ordered = normal_order(term);

    std::vector<Term> terms = {
        Term(-1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                       Operator::Spin::UP, 1),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::FERMION, Operator::Spin::UP, 0)})};
    Expression expected(terms);

    EXPECT_EQ(normal_ordered, expected);
  }

  {
    Term term(1.0,
              {Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                        Operator::Spin::UP, 1),
               Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                        Operator::Spin::UP, 0),
               Operator(Operator::Type::ANNIHILATION,
                        Operator::Statistics::BOSON, Operator::Spin::UP, 0)});
    Expression normal_ordered = normal_order(term);

    std::vector<Term> terms = {Term(
        1.0, {Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 1),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::BOSON, Operator::Spin::UP, 0)})};
    Expression expected(terms);

    EXPECT_EQ(normal_ordered, expected);
  }
}

TEST(NormalOrderTest, NormalOrderTermCreationAnnihilationCreationSameOrbital) {
  {
    Term term(1.0,
              {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                        Operator::Spin::UP, 0),
               Operator(Operator::Type::ANNIHILATION,
                        Operator::Statistics::FERMION, Operator::Spin::UP, 0),
               Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                        Operator::Spin::UP, 0)});
    Expression normal_ordered = normal_order(term);

    std::vector<Term> terms = {
        Term(-1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::FERMION, Operator::Spin::UP, 0)}),
        Term(1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                       Operator::Spin::UP, 0)})};
    Expression expected(terms);

    EXPECT_EQ(normal_ordered, expected);
  }

  {
    Term term(1.0,
              {Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                        Operator::Spin::UP, 0),
               Operator(Operator::Type::ANNIHILATION,
                        Operator::Statistics::BOSON, Operator::Spin::UP, 0),
               Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                        Operator::Spin::UP, 0)});
    Expression normal_ordered = normal_order(term);

    std::vector<Term> terms = {
        Term(1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::BOSON, Operator::Spin::UP, 0)}),
        Term(1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 0)})};
    Expression expected(terms);

    EXPECT_EQ(normal_ordered, expected);
  }
}

TEST(NormalOrderTest,
     NormalOrderTermAnnihilationCreationAnnihilationSameOrbital) {
  {
    Term term(1.0,
              {Operator(Operator::Type::ANNIHILATION,
                        Operator::Statistics::FERMION, Operator::Spin::UP, 0),
               Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                        Operator::Spin::UP, 0),
               Operator(Operator::Type::ANNIHILATION,
                        Operator::Statistics::FERMION, Operator::Spin::UP, 0)});
    Expression normal_ordered = normal_order(term);

    std::vector<Term> terms = {
        Term(-1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::FERMION, Operator::Spin::UP, 0),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::FERMION, Operator::Spin::UP, 0)}),
        Term(1.0,
             {Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::FERMION, Operator::Spin::UP, 0)})};
    Expression expected(terms);

    EXPECT_EQ(normal_ordered, expected);
  }

  {
    Term term(1.0,
              {Operator(Operator::Type::ANNIHILATION,
                        Operator::Statistics::BOSON, Operator::Spin::UP, 0),
               Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                        Operator::Spin::UP, 0),
               Operator(Operator::Type::ANNIHILATION,
                        Operator::Statistics::BOSON, Operator::Spin::UP, 0)});
    Expression normal_ordered = normal_order(term);

    std::vector<Term> terms = {
        Term(1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::BOSON, Operator::Spin::UP, 0),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::BOSON, Operator::Spin::UP, 0)}),
        Term(1.0,
             {Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::BOSON, Operator::Spin::UP, 0)})};
    Expression expected(terms);

    EXPECT_EQ(normal_ordered, expected);
  }
}

TEST(NormalOrderTest,
     NormalOrderTermAnnihilationAnnihilationCreationSameOrbital) {
  {
    Term term(1.0,
              {Operator(Operator::Type::ANNIHILATION,
                        Operator::Statistics::FERMION, Operator::Spin::UP, 0),
               Operator(Operator::Type::ANNIHILATION,
                        Operator::Statistics::FERMION, Operator::Spin::UP, 0),
               Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                        Operator::Spin::UP, 0)});
    Expression normal_ordered = normal_order(term);

    std::vector<Term> terms = {
        Term(1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::FERMION, Operator::Spin::UP, 0),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::FERMION, Operator::Spin::UP, 0)}),
        Term(0.0,
             {Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::FERMION, Operator::Spin::UP, 0)})};
    Expression expected(terms);

    EXPECT_EQ(normal_ordered, expected);
  }

  {
    Term term(1.0,
              {Operator(Operator::Type::ANNIHILATION,
                        Operator::Statistics::BOSON, Operator::Spin::UP, 0),
               Operator(Operator::Type::ANNIHILATION,
                        Operator::Statistics::BOSON, Operator::Spin::UP, 0),
               Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                        Operator::Spin::UP, 0)});
    Expression normal_ordered = normal_order(term);

    std::vector<Term> terms = {
        Term(1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::BOSON, Operator::Spin::UP, 0),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::BOSON, Operator::Spin::UP, 0)}),
        Term(2.0,
             {Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::BOSON, Operator::Spin::UP, 0)})};
    Expression expected(terms);

    EXPECT_EQ(normal_ordered, expected);
  }
}

TEST(NormalOrderTest, NormalOrderExpression) {
  {
    std::vector<Term> terms = {
        Term(1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::FERMION, Operator::Spin::UP, 1)}),
        Term(1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                       Operator::Spin::UP, 1),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::FERMION, Operator::Spin::UP, 0)})};
    Expression normal_ordered = normal_order(terms);

    std::vector<Term> normal_terms = {
        Term(1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::FERMION, Operator::Spin::UP, 1)}),
        Term(1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                       Operator::Spin::UP, 1),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::FERMION, Operator::Spin::UP, 0)})};
    Expression expected(normal_terms);

    EXPECT_EQ(normal_ordered, expected);
  }

  {
    std::vector<Term> terms = {
        Term(1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::BOSON, Operator::Spin::UP, 1)}),
        Term(1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 1),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::BOSON, Operator::Spin::UP, 0)})};
    Expression normal_ordered = normal_order(terms);

    std::vector<Term> normal_terms = {
        Term(1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::BOSON, Operator::Spin::UP, 1)}),
        Term(1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 1),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::BOSON, Operator::Spin::UP, 0)})};
    Expression expected(normal_terms);

    EXPECT_EQ(normal_ordered, expected);
  }
}

TEST(NormalOrderTest, NormalOrderExpressionWrongOrder) {
  {
    std::vector<Term> terms = {
        Term(1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                       Operator::Spin::UP, 1),
              Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                       Operator::Spin::UP, 0)}),
        Term(1.0,
             {Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::FERMION, Operator::Spin::UP, 0),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::FERMION, Operator::Spin::UP, 1)})};
    Expression normal_ordered = normal_order(terms);

    std::vector<Term> normal_terms = {
        Term(-1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                       Operator::Spin::UP, 1)}),
        Term(-1.0,
             {Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::FERMION, Operator::Spin::UP, 1),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::FERMION, Operator::Spin::UP, 0)})};
    Expression expected(normal_terms);

    EXPECT_EQ(normal_ordered, expected);
  }

  {
    std::vector<Term> terms = {
        Term(1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 1),
              Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 0)}),
        Term(1.0,
             {Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::BOSON, Operator::Spin::UP, 0),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::BOSON, Operator::Spin::UP, 1)})};
    Expression normal_ordered = normal_order(terms);

    std::vector<Term> normal_terms = {
        Term(1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 1)}),
        Term(1.0,
             {Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::BOSON, Operator::Spin::UP, 1),
              Operator(Operator::Type::ANNIHILATION,
                       Operator::Statistics::BOSON, Operator::Spin::UP, 0)})};
    Expression expected(normal_terms);

    EXPECT_EQ(normal_ordered, expected);
  }
}

TEST(NormalOrderTest, NormalOrderExpressionResultingInZero) {
  {
    std::vector<Term> terms = {
        Term(1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                       Operator::Spin::UP, 1),
              Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                       Operator::Spin::UP, 0)}),
        Term(1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                       Operator::Spin::UP, 1)})};
    Expression normal_ordered = normal_order(terms);

    std::vector<Term> normal_terms = {Term(
        0.0, {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                       Operator::Spin::UP, 1)})};
    Expression expected(normal_terms);

    EXPECT_EQ(normal_ordered, expected);
  }

  {
    std::vector<Term> terms = {
        Term(1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 1),
              Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 0)}),
        Term(1.0,
             {Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 1)})};
    Expression normal_ordered = normal_order(terms);

    std::vector<Term> normal_terms = {Term(
        2.0, {Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 0),
              Operator(Operator::Type::CREATION, Operator::Statistics::BOSON,
                       Operator::Spin::UP, 1)})};
    Expression expected(normal_terms);

    EXPECT_EQ(normal_ordered, expected);
  }
}

TEST(NormalOrderTest,
     NormalOrderTermAnnihilationAnnihilationCreationSameOrbitalAfterClean) {
  Term term(1.0,
            {Operator(Operator::Type::ANNIHILATION,
                      Operator::Statistics::FERMION, Operator::Spin::UP, 0),
             Operator(Operator::Type::ANNIHILATION,
                      Operator::Statistics::FERMION, Operator::Spin::UP, 0),
             Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                      Operator::Spin::UP, 0)});
  Expression normal_ordered = normal_order(term);
  std::erase_if(normal_ordered.terms(),
                [](const auto &term) { return std::abs(term.second) < 1e-10; });

  std::vector<Term> terms = {Term(
      1.0, {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                     Operator::Spin::UP, 0),
            Operator(Operator::Type::ANNIHILATION,
                     Operator::Statistics::FERMION, Operator::Spin::UP, 0),
            Operator(Operator::Type::ANNIHILATION,
                     Operator::Statistics::FERMION, Operator::Spin::UP, 0)})};
  Expression expected(terms);

  EXPECT_EQ(normal_ordered, expected);
}

TEST(NormalOrderTest, NormalOrderExpressionResultingInZeroAfterClean) {
  std::vector<Term> terms = {
      Term(1.0,
           {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                     Operator::Spin::UP, 1),
            Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                     Operator::Spin::UP, 0)}),
      Term(1.0,
           {Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                     Operator::Spin::UP, 0),
            Operator(Operator::Type::CREATION, Operator::Statistics::FERMION,
                     Operator::Spin::UP, 1)})};
  Expression normal_ordered = normal_order(terms);
  std::erase_if(normal_ordered.terms(),
                [](const auto &term) { return std::abs(term.second) < 1e-10; });
  EXPECT_THAT(normal_ordered.terms(), IsEmpty());
}
