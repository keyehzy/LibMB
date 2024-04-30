#include "Expression.h"

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using testing::IsEmpty;

TEST(ExpressionTest, ConstructorAndAccessors) {
  std::vector<Term> terms = {
      Term(2.5, {Operator(OperatorType::CREATION, Spin::UP, 0),
                 Operator(OperatorType::ANNIHILATION, Spin::DOWN, 1)}),
      Term(3.0, {Operator(OperatorType::CREATION, Spin::DOWN, 0),
                 Operator(OperatorType::ANNIHILATION, Spin::UP, 1)})};
  Expression expression(terms);

  EXPECT_EQ(expression.terms().size(), 2);
  EXPECT_DOUBLE_EQ(expression.terms().at(
                       {Operator(OperatorType::CREATION, Spin::UP, 0),
                        Operator(OperatorType::ANNIHILATION, Spin::DOWN, 1)}),
                   2.5);
  EXPECT_DOUBLE_EQ(expression.terms().at(
                       {Operator(OperatorType::CREATION, Spin::DOWN, 0),
                        Operator(OperatorType::ANNIHILATION, Spin::UP, 1)}),
                   3.0);
}

TEST(ExpressionTest, EqualityOperator) {
  std::vector<Term> terms1 = {
      Term(2.5, {Operator(OperatorType::CREATION, Spin::UP, 0),
                 Operator(OperatorType::ANNIHILATION, Spin::DOWN, 1)}),
      Term(3.0, {Operator(OperatorType::CREATION, Spin::DOWN, 0),
                 Operator(OperatorType::ANNIHILATION, Spin::UP, 1)})};
  Expression expression1(terms1);

  std::vector<Term> terms2 = {
      Term(2.5, {Operator(OperatorType::CREATION, Spin::UP, 0),
                 Operator(OperatorType::ANNIHILATION, Spin::DOWN, 1)}),
      Term(3.0, {Operator(OperatorType::CREATION, Spin::DOWN, 0),
                 Operator(OperatorType::ANNIHILATION, Spin::UP, 1)})};
  Expression expression2(terms2);

  std::vector<Term> terms3 = {
      Term(3.0, {Operator(OperatorType::CREATION, Spin::DOWN, 0),
                 Operator(OperatorType::ANNIHILATION, Spin::UP, 1)}),
      Term(2.5, {Operator(OperatorType::CREATION, Spin::UP, 0),
                 Operator(OperatorType::ANNIHILATION, Spin::DOWN, 1)})};
  Expression expression3(terms3);

  std::vector<Term> terms4 = {
      Term(2.5, {Operator(OperatorType::CREATION, Spin::UP, 0),
                 Operator(OperatorType::ANNIHILATION, Spin::DOWN, 1)}),
      Term(3.0, {Operator(OperatorType::CREATION, Spin::DOWN, 0),
                 Operator(OperatorType::ANNIHILATION, Spin::UP, 1)}),
      Term(1.0, {Operator(OperatorType::CREATION, Spin::UP, 0),
                 Operator(OperatorType::ANNIHILATION, Spin::DOWN, 1)})};
  Expression expression4(terms4);

  EXPECT_TRUE(expression1 == expression2);
  EXPECT_TRUE(expression1 == expression3);
  EXPECT_FALSE(expression1 == expression4);
}

TEST(ExpressionTest, InequalityOperator) {
  std::vector<Term> terms1 = {
      Term(2.5, {Operator(OperatorType::CREATION, Spin::UP, 0),
                 Operator(OperatorType::ANNIHILATION, Spin::DOWN, 1)}),
      Term(3.0, {Operator(OperatorType::CREATION, Spin::DOWN, 0),
                 Operator(OperatorType::ANNIHILATION, Spin::UP, 1)})};
  Expression expression1(terms1);

  std::vector<Term> terms2 = {
      Term(2.5, {Operator(OperatorType::CREATION, Spin::UP, 0),
                 Operator(OperatorType::ANNIHILATION, Spin::DOWN, 1)}),
      Term(3.0, {Operator(OperatorType::CREATION, Spin::DOWN, 0),
                 Operator(OperatorType::ANNIHILATION, Spin::UP, 1)})};
  Expression expression2(terms2);

  std::vector<Term> terms3 = {
      Term(3.0, {Operator(OperatorType::CREATION, Spin::DOWN, 0),
                 Operator(OperatorType::ANNIHILATION, Spin::UP, 1)}),
      Term(2.5, {Operator(OperatorType::CREATION, Spin::UP, 0),
                 Operator(OperatorType::ANNIHILATION, Spin::DOWN, 1)})};
  Expression expression3(terms3);

  std::vector<Term> terms4 = {
      Term(2.5, {Operator(OperatorType::CREATION, Spin::UP, 0),
                 Operator(OperatorType::ANNIHILATION, Spin::DOWN, 1)}),
      Term(3.0, {Operator(OperatorType::CREATION, Spin::DOWN, 0),
                 Operator(OperatorType::ANNIHILATION, Spin::UP, 1)}),
      Term(1.0, {Operator(OperatorType::CREATION, Spin::UP, 0),
                 Operator(OperatorType::ANNIHILATION, Spin::DOWN, 1)})};
  Expression expression4(terms4);

  EXPECT_FALSE(expression1 != expression2);
  EXPECT_FALSE(expression1 != expression3);
  EXPECT_TRUE(expression1 != expression4);
}

TEST(NormalOrderTest, ExpressionResultingInZero) {
  std::vector<Term> terms = {
      Term(1.0, {Operator(OperatorType::CREATION, Spin::UP, 0),
                 Operator(OperatorType::ANNIHILATION, Spin::UP, 1)}),
      Term(-1.0, {Operator(OperatorType::CREATION, Spin::UP, 0),
                  Operator(OperatorType::ANNIHILATION, Spin::UP, 1)})};
  Expression expression(terms);

  std::vector<Term> normal_terms = {
      Term(0.0, {Operator(OperatorType::CREATION, Spin::UP, 0),
                 Operator(OperatorType::ANNIHILATION, Spin::UP, 1)})};
  Expression expected(normal_terms);

  EXPECT_EQ(expression, expected);
}

TEST(NormalOrderTest, ExpressionResultingInZeroAfterClean) {
  std::vector<Term> terms = {
      Term(1.0, {Operator(OperatorType::CREATION, Spin::UP, 0),
                 Operator(OperatorType::ANNIHILATION, Spin::UP, 1)}),
      Term(-1.0, {Operator(OperatorType::CREATION, Spin::UP, 0),
                  Operator(OperatorType::ANNIHILATION, Spin::UP, 1)})};
  Expression expression(terms);
  expression.clean();

  EXPECT_THAT(expression.terms(), IsEmpty());
}