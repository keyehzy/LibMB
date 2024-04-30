#include "Basis.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <unordered_set>

using testing::UnorderedElementsAre;

constexpr int binomial(int n, int k) {
  if (k == 0 || k == n) {
    return 1;
  }
  return binomial(n - 1, k - 1) + binomial(n - 1, k);
}

TEST(BasisTest, ConstructorAndAttributes) {
  Basis basis(2, 2);

  EXPECT_EQ(basis.particles(), 2);
  EXPECT_EQ(basis.orbitals(), 2);
  EXPECT_EQ(basis.elements().size(), binomial(4, 2));
}

TEST(BasisTest, BasisGeneration) {
  Basis basis(2, 2);
  EXPECT_EQ(basis.elements().size(), 6);
  EXPECT_THAT(basis.elements(),
              UnorderedElementsAre(
                  std::vector<Operator>{Operator::creation(Spin::DOWN, 0),
                                        Operator::creation(Spin::DOWN, 1)},
                  std::vector<Operator>{Operator::creation(Spin::UP, 1),
                                        Operator::creation(Spin::DOWN, 1)},
                  std::vector<Operator>{Operator::creation(Spin::DOWN, 0),
                                        Operator::creation(Spin::UP, 1)},
                  std::vector<Operator>{Operator::creation(Spin::UP, 0),
                                        Operator::creation(Spin::DOWN, 1)},
                  std::vector<Operator>{Operator::creation(Spin::UP, 0),
                                        Operator::creation(Spin::DOWN, 0)},
                  std::vector<Operator>{Operator::creation(Spin::UP, 0),
                                        Operator::creation(Spin::UP, 1)}));
}

TEST(BasisTest, IndexingUnique) {
  Basis basis(2, 2);

  std::unordered_set<std::size_t> indices;
  for (const auto& term : basis.elements()) {
    indices.insert(basis.index(term));
  }
  EXPECT_EQ(indices.size(), basis.elements().size());

  for (const auto& term : basis.elements()) {
    EXPECT_GE(basis.index(term), 0);
    EXPECT_LT(basis.index(term), basis.elements().size());
  }
}

TEST(BasisTest, IndexingInsideBasis) {
  Basis basis(2, 2);

  std::vector<Operator> term = {Operator::creation(Spin::UP, 0),
                                Operator::creation(Spin::DOWN, 1)};
  std::size_t index = basis.index(term);
  EXPECT_TRUE(basis.contains(term));
}

TEST(BasisTest, IndexingOutsideBasis) {
  Basis basis(2, 2);

  std::vector<Operator> term = {Operator::creation(Spin::UP, 0),
                                Operator::creation(Spin::UP, 2)};
  EXPECT_FALSE(basis.contains(term));
}

TEST(BasisTest, EqualityOperator) {
  Basis basis1(2, 2);
  Basis basis2(2, 2);
  Basis basis3(2, 3);
  Basis basis4(3, 2);

  EXPECT_EQ(basis1, basis2);
  EXPECT_NE(basis1, basis3);
  EXPECT_NE(basis1, basis4);
}

TEST(BasisTest, IndexingEmptyTerm) {
  Basis basis(2, 2);

  std::vector<Operator> term = {};
  EXPECT_FALSE(basis.contains(term));
}

TEST(BasisTest, IndexingSingleTermSingleBodyBasis) {
  Basis basis(1, 1);

  std::vector<Operator> term = {Operator::creation(Spin::UP, 0)};
  EXPECT_TRUE(basis.contains(term));
}

TEST(BasisTest, IndexingSingleTermManyBodyBasis) {
  Basis basis(2, 2);

  std::vector<Operator> term = {Operator::creation(Spin::UP, 0)};
  EXPECT_FALSE(basis.contains(term));
}